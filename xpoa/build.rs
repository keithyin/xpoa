use std::{
    env,
    path::{Path, PathBuf},
    process::{Command, ExitStatus},
};
fn check_command_status(status: ExitStatus, msg: &str) {
    if !status.success() {
        panic!("check command status fail. {}", msg);
    }
}

fn main() {
    // 纯 Rust 构建(default-features=false):跳过 C++,无需 cmake/Boost/SeqAn
    if std::env::var("CARGO_FEATURE_CPP").is_err() {
        return;
    }

    let build_option = "Release";

    let cpp_code_dir_name = "cpp_code";
    let out_path_str = env::var("OUT_DIR").unwrap();
    let build_out_path = Path::new(&out_path_str);
    // let c_cur_dir = env::current_dir().unwrap();
    // let build_out_path = Path::new(c_cur_dir.to_str().unwrap());

    Command::new("sh")
        .arg("-c")
        .arg(format!(
            "cp -r {} {}",
            cpp_code_dir_name,
            build_out_path.to_str().unwrap()
        ))
        .output()
        .expect("cp cpp_code error");

    let c_src_file_dir = build_out_path.join(cpp_code_dir_name);
    let c_src_file_dir = &c_src_file_dir;
    let current_dir = env::current_dir().unwrap();

    // third_party 头文件(优先于系统路径);见 third_party/README.md。
    // 查找顺序:环境变量 XPOA_BOOST_INCLUDE / XPOA_SEQAN_INCLUDE(供
    // `cargo package`/`cargo publish` 的 verify 构建使用 —— 其工作目录是
    // 解包目录,cwd 相对路径会失效)> cwd 相对的 third_party/。
    let mut third_party_includes = String::new();
    let boost_dir = env::var("XPOA_BOOST_INCLUDE")
        .ok()
        .map(PathBuf::from)
        .unwrap_or_else(|| current_dir.join("third_party/boost/include"));
    let seqan_dir = env::var("XPOA_SEQAN_INCLUDE")
        .ok()
        .map(PathBuf::from)
        .unwrap_or_else(|| current_dir.join("third_party/seqan/include"));
    for dir in [boost_dir, seqan_dir] {
        if dir.is_dir() {
            third_party_includes.push_str(&format!(" -I{}", dir.to_str().unwrap()));
        }
    }

    Command::new("sh")
        .arg("-c")
        .arg("rm -rf build")
        .current_dir(format!("{}/{}", out_path_str, cpp_code_dir_name))
        .status()
        .expect("remove build error");

    check_command_status(
        Command::new("sh")
            .arg("-c")
            .arg(format!(
                "/bin/cmake -DCMAKE_BUILD_TYPE:STRING={} -S./ -B./build -G 'Unix Makefiles' -DCMAKE_CXX_FLAGS='{}'",
                build_option, third_party_includes
            ))
            .current_dir(format!("{}/{}", out_path_str, cpp_code_dir_name))
            .status()
            .expect("cmake config error"),
        "cmake config",
    );

    check_command_status(
        Command::new("sh")
            .arg("-c")
            .arg(format!(
                "/bin/cmake --build build/ --config {} --target all -j{} --",
                build_option,
                num_cpus::get_physical().max(4)
            ))
            .current_dir(format!("{}/{}", out_path_str, cpp_code_dir_name))
            .status()
            .expect("first build"),
        "build",
    );

    let obj_filedir = c_src_file_dir.join("build");
    let obj_filedir = &obj_filedir;

    println!("cargo:rerun-if-changed={}", current_dir.to_str().unwrap());
    println!("cargo:rerun-if-changed={}", current_dir.to_str().unwrap());
    println!(
        "cargo:rustc-link-search=native={}",
        obj_filedir.to_str().unwrap()
    );
    println!("cargo:rustc-link-lib=static=xpoa");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=dylib=stdc++");

    // println!("cargo:rustc-link-lib=pthread");
}
