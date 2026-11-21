//! xpoa-cli:按分组生成 POA consensus 的命令行工具。
//!
//! 输入 FASTA 中每条序列的 qname 以 `|group=GROUP_ID` 结尾;相同 GROUP_ID 的
//! 序列必须连续出现(否则 panic)。每组经 POA 生成 consensus,以 GROUP_ID 为
//! qname 输出到 FASTA。
//!
//! 流水线:reader 线程解析 + 分组 → rayon 并行计算(在飞组数有界)→ 主线程
//! 按输入组顺序写出;进度条输出到 stderr。

use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::thread;
use std::time::Duration;

use clap::{Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::ThreadPoolBuilder;
use xpoa::error::XpoaResult;
use xpoa::xpoa::get_default_poa_setting;
use xpoa::{BACKEND_CPP, BACKEND_RUST, PoaSetting, poa_consensus};

/// reader → 主线程的消息;Eof 用于区分正常读完与 reader 异常退出。
enum Msg {
    Group(Group),
    Eof,
}

/// 一个连续分组的全部 reads。
struct Group {
    id: String,
    reads: Vec<String>,
    /// 组内首条 header 的行号,报错定位用。
    first_line: u64,
}

/// 单组计算结果,回传主线程按序写出:(组 id, consensus, reads 数, 组首行号)。
type ComputeResult = (String, XpoaResult<Option<(String, usize)>>, usize, u64);

/// reader 组通道容量:限制同时缓存的组数,内存有界。
const GROUP_CHANNEL_CAP: usize = 16;

#[derive(Debug, Parser)]
#[command(
    name = "xpoa-cli",
    version,
    about = "POA consensus per group from a grouped FASTA"
)]
struct Cli {
    /// 输入 FASTA:每个 qname 必须以 "|group=GROUP_ID" 结尾,且同组序列必须连续
    input: PathBuf,

    /// 输出 FASTA(缺省 stdout)
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// consensus 工作线程数(0 = 全核)
    #[arg(short, long, default_value_t = 0, value_name = "N")]
    threads: usize,

    /// POA 后端(缺省 = 编译期默认;纯 Rust 构建下选 cpp 会报错)
    #[arg(long, value_enum, value_name = "BACKEND")]
    backend: Option<BackendArg>,

    /// POA 算法版本(1 或 2)
    #[arg(long, default_value_t = 1, value_name = "1|2")]
    poa_version: i32,

    /// 输出 FASTA 折行宽度(0 = 单行)
    #[arg(short = 'w', long, default_value_t = 60, value_name = "N")]
    line_width: usize,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum BackendArg {
    Cpp,
    Rust,
}

impl BackendArg {
    fn backend(self) -> i32 {
        match self {
            Self::Cpp => BACKEND_CPP,
            Self::Rust => BACKEND_RUST,
        }
    }
}

#[derive(Default)]
struct Stats {
    groups: u64,
    reads: u64,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> Result<(), String> {
    if cli.poa_version != 1 && cli.poa_version != 2 {
        return Err(format!(
            "invalid --poa-version {} (expected 1 or 2)",
            cli.poa_version
        ));
    }
    if cli.threads > 0 {
        ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .map_err(|e| format!("build rayon thread pool: {e}"))?;
    }
    let mut setting = get_default_poa_setting();
    if let Some(backend) = cli.backend {
        setting.backend = backend.backend();
    }
    setting.version = cli.poa_version;

    let pb = new_spinner();
    match pipeline(cli, &setting, &pb) {
        Ok(stats) => {
            pb.finish_with_message(format!(
                "done: {} groups, {} reads",
                stats.groups, stats.reads
            ));
            Ok(())
        }
        Err(e) => {
            pb.abandon();
            Err(e)
        }
    }
}

/// 读 ‖ 算 ‖ 写流水线:收组 → 派发 rayon 计算 → 超出在飞上限时按序写出。
fn pipeline(cli: &Cli, setting: &PoaSetting, pb: &ProgressBar) -> Result<Stats, String> {
    let (tx, rx) = sync_channel::<Msg>(GROUP_CHANNEL_CAP);
    let input = cli.input.clone();
    let reader = thread::Builder::new()
        .name("fasta-reader".into())
        .spawn(move || read_fasta_groups(&input, tx))
        .map_err(|e| format!("spawn reader thread: {e}"))?;

    let mut out: Box<dyn Write> = match &cli.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).map_err(|e| format!("create {}: {e}", path.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };

    // 有界在飞:队列超过线程数时先按序写出队首;内存 ≈ 线程数 + 组通道容量。
    let max_inflight = rayon::current_num_threads().max(1);
    let mut queue: VecDeque<Receiver<ComputeResult>> = VecDeque::new();
    let mut stats = Stats::default();

    loop {
        match rx.recv() {
            Ok(Msg::Group(group)) => {
                let (rtx, rrx) = channel();
                let setting = *setting;
                rayon::spawn(move || {
                    let res = poa_consensus(&group.reads, &setting);
                    let _ = rtx.send((group.id, res, group.reads.len(), group.first_line));
                });
                queue.push_back(rrx);
                if queue.len() > max_inflight {
                    write_result(
                        out.as_mut(),
                        queue.pop_front().expect("queue non-empty"),
                        pb,
                        cli.line_width,
                        &mut stats,
                    )?;
                }
            }
            Ok(Msg::Eof) => break,
            Err(_) => {
                // Eof 前通道关闭:reader 提前返回或 panic
                join_reader(reader, pb)?;
                return Err("reader thread ended without Eof".into());
            }
        }
    }
    join_reader(reader, pb)?;

    for rrx in queue.drain(..) {
        write_result(out.as_mut(), rrx, pb, cli.line_width, &mut stats)?;
    }
    out.flush().map_err(|e| format!("flush output: {e}"))?;
    Ok(stats)
}

/// reader 线程:解析 FASTA,按连续分组切批后经通道发出;结束时发 Eof。
/// 同组不连续 → panic(主线程 join 后 resume_unwind 原样传播)。
fn read_fasta_groups(path: &Path, tx: SyncSender<Msg>) -> Result<(), String> {
    let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let mut seen: HashSet<String> = HashSet::new();
    let mut cur_id: Option<String> = None;
    let mut cur_group_line = 0u64;
    let mut cur_qname = String::new();
    let mut cur_record_line = 0u64;
    let mut cur_seq = String::new();
    let mut reads: Vec<String> = Vec::new();
    let mut lineno = 0u64;

    let send_group = |group| {
        tx.send(Msg::Group(group))
            .map_err(|_| "output closed; reader aborting".to_string())
    };

    for line in BufReader::new(file).lines() {
        lineno += 1;
        let line = line.map_err(|e| format!("read {}: line {lineno}: {e}", path.display()))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            // 收尾上一条 record:必须有序列
            if cur_id.is_some() {
                if cur_seq.is_empty() {
                    return Err(format!(
                        "line {lineno}: record \"{cur_qname}\" (header at line {cur_record_line}) has no sequence"
                    ));
                }
                reads.push(std::mem::take(&mut cur_seq));
            }
            let qname = header.split_whitespace().next().unwrap_or("");
            let gid = parse_group(qname, lineno)?;
            if cur_id.as_deref() != Some(gid) {
                // 组切换:该组若之前出现过 → 非连续,panic
                if let Some(prev) = &cur_id {
                    if seen.contains(gid) {
                        panic!(
                            "group \"{gid}\" (record \"{qname}\" at line {lineno}) reappears after \
                             group \"{prev}\"; same-group records must be contiguous in the input FASTA"
                        );
                    }
                    let finished = Group {
                        id: prev.clone(),
                        reads: std::mem::take(&mut reads),
                        first_line: cur_group_line,
                    };
                    send_group(finished)?;
                }
                seen.insert(gid.to_string());
                cur_id = Some(gid.to_string());
                cur_group_line = lineno;
            }
            cur_qname = qname.to_string();
            cur_record_line = lineno;
        } else {
            if cur_id.is_none() {
                return Err(format!(
                    "line {lineno}: sequence data before the first FASTA header"
                ));
            }
            cur_seq.push_str(&line.to_ascii_uppercase());
        }
    }

    // EOF:收尾最后一条 record 与最后一个组
    if cur_id.is_some() {
        if cur_seq.is_empty() {
            return Err(format!(
                "record \"{cur_qname}\" (header at line {cur_record_line}) has no sequence"
            ));
        }
        reads.push(cur_seq);
    }
    if let Some(id) = cur_id {
        send_group(Group {
            id,
            reads,
            first_line: cur_group_line,
        })?;
    }
    tx.send(Msg::Eof).map_err(|_| "output closed".to_string())?;
    Ok(())
}

/// qname → GROUP_ID:取最后一个 "|group=" 之后的内容(即以该后缀结尾)。
fn parse_group(qname: &str, lineno: u64) -> Result<&str, String> {
    match qname.rsplit_once("|group=") {
        Some((_, gid)) if !gid.is_empty() => Ok(gid),
        Some(_) => Err(format!(
            "line {lineno}: empty GROUP_ID in qname \"{qname}\""
        )),
        None => Err(format!(
            "line {lineno}: qname \"{qname}\" does not end with \"|group=GROUP_ID\""
        )),
    }
}

/// join reader 线程:Err 返回值原样传播;panic 经 resume_unwind 原样重抛。
fn join_reader(
    reader: thread::JoinHandle<Result<(), String>>,
    pb: &ProgressBar,
) -> Result<(), String> {
    match reader.join() {
        Ok(Ok(())) => Ok(()),
        Ok(Err(msg)) => Err(msg),
        Err(payload) => {
            pb.abandon();
            std::panic::resume_unwind(payload);
        }
    }
}

/// 按序取出一组计算结果并写出;更新进度与统计。
fn write_result(
    out: &mut dyn Write,
    rx: Receiver<ComputeResult>,
    pb: &ProgressBar,
    line_width: usize,
    stats: &mut Stats,
) -> Result<(), String> {
    let (gid, res, n_reads, first_line) = rx
        .recv()
        .map_err(|_| "consensus worker died without returning a result".to_string())?;
    match res {
        Ok(Some((seq, _n_passes))) => {
            write_record(out, &gid, &seq, line_width)
                .map_err(|e| format!("write group \"{gid}\": {e}"))?;
        }
        Ok(None) => {
            return Err(format!(
                "group \"{gid}\" (line {first_line}): produced no consensus"
            ));
        }
        Err(e) => {
            return Err(format!(
                "group \"{gid}\" (line {first_line}): poa_consensus failed: {e}"
            ));
        }
    }
    stats.groups += 1;
    stats.reads += n_reads as u64;
    pb.inc(1);
    pb.set_message(format!("groups: {}  reads: {}", stats.groups, stats.reads));
    Ok(())
}

/// 写一条 FASTA 记录:>{id}\n + 按 line_width 折行(0 = 单行)。
fn write_record(
    out: &mut dyn Write,
    id: &str,
    seq: &str,
    line_width: usize,
) -> std::io::Result<()> {
    out.write_all(b">")?;
    out.write_all(id.as_bytes())?;
    out.write_all(b"\n")?;
    if line_width == 0 {
        out.write_all(seq.as_bytes())?;
        out.write_all(b"\n")?;
    } else {
        for chunk in seq.as_bytes().chunks(line_width) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn new_spinner() -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} {msg} [{elapsed_precise}]")
            .expect("spinner template is valid"),
    );
    pb.enable_steady_tick(Duration::from_millis(200));
    pb
}
