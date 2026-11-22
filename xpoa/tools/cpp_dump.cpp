// 调试工具:用 C++ 库重建指定 case 的 PoaDraft 流程,输出每步的图 CSV
// 与最终共识/摘要,供与 Rust native 实现差分对比。
// 独立编译(不属于 CMake 目标):
//   g++ -std=c++17 -O2 -I cpp_code -I cpp_code/edlib \
//     -I third_party/boost/include -I third_party/seqan/include \
//     tools/cpp_dump.cpp -o /tmp/cpp_dump \
//     target/release/build/xpoa-*/out/cpp_code/build/libxpoa.a -lz -lm

#include <Settings.h>
#include <SparsePoa.h>
#include <Sequence.h>
#include <EdlibGlobalAlign.h>
#include <PoaDraft.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using SparsePoa = PacBio::Poa::SparsePoa;
using Summary = PacBio::Poa::PoaAlignmentSummary;

static const int ADAPTER_BEFORE = 1;
static const int ADAPTER_AFTER = 2;

int main(int argc, char** argv)
{
    // case_004 的输入(grid_l50_e0.05_n3_v1_u0)
    std::vector<std::string> reads = {
        "CCCAAGAGGGACTTCAGACGGACGCGGGGCGGGTAGTTGTATCTAACCAGG",
        "CCCAAGAGGGACTTCAGAACTGACGCGGGGCGGGTAGTTGTATCTAACCACGG",
        "CCCAAGAGGGACTTCAAACGACGCGGGGCGTGGTAGGTTGTATCTAACCAGG",
    };
    std::vector<int> flags = {3, 3, 3};

    PoaSetting setting;
    setting.version = 1;
    setting.ed_unify_strand = 0;

    if (argc > 1) setting.version = std::atoi(argv[1]);

    SparsePoa poa(setting.version);
    std::vector<SparsePoa::ReadKey> keys;
    size_t cov = 0;

    for (size_t i = 0; i < reads.size(); i++) {
        const std::string& seq = reads[i];
        SparsePoa::ReadKey key;
        if (setting.ed_unify_strand) {
            if (i == 0) {
                key = poa.AddRead(seq, setting, true);
            } else {
                std::string rc = PacBio::Data::ReverseComplement(seq);
                int d1 = Geneus::Align::EdlibGlobalAlignDistance(seq, reads[0]);
                int d2 = Geneus::Align::EdlibGlobalAlignDistance(rc, reads[0]);
                key = (d1 > d2) ? poa.AddRead(rc, setting, true) : poa.AddRead(seq, setting, false);
            }
        } else {
            key = poa.OrientAndAddRead(seq, setting);
        }
        keys.push_back(key);
        if (key >= 0) cov++;

        char name[128];
        std::snprintf(name, sizeof(name), "cpp_graph_after_%zu.csv", i);
        poa.WriteGraphCsvFile(name);
        std::cerr << "read" << i << " key=" << key << std::endl;
    }

    const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
    std::vector<Summary> summaries;
    auto pc = poa.FindConsensus(static_cast<int>(minCov), setting, &summaries);
    poa.WriteGraphCsvFile("cpp_graph_final.csv");

    std::cout << "consensus: " << pc->Sequence << std::endl;
    for (size_t i = 0; i < summaries.size(); i++) {
        const auto& s = summaries[i];
        std::cout << "read" << i << " identity=" << s.AlignmentIdentity
                  << " rc=" << s.ReverseComplementedRead << " read=[" << s.ExtentOnRead.Left()
                  << "," << s.ExtentOnRead.Right() << ") css=[" << s.ExtentOnConsensus.Left() << ","
                  << s.ExtentOnConsensus.Right() << ")" << std::endl;
    }
    size_t n_passes = 0;
    for (size_t i = 0; i < reads.size(); i++) {
        if (keys[i] < 0) continue;
        if (summaries[keys[i]].AlignmentIdentity < setting.min_identity) continue;
        if (!(flags[i] & ADAPTER_BEFORE && flags[i] & ADAPTER_AFTER)) continue;
        n_passes++;
    }
    std::cout << "n_passes: " << n_passes << std::endl;
    return 0;
}
