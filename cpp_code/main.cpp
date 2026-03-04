#include <iostream>
#include "PoaConsensus.h"
#include "SparsePoa.h"
#include "Sequence.h"
#include "EdlibGlobalAlign.h"
// using namespace PacBio;
using SparsePoa = PacBio::Poa::SparsePoa;
using PoaAlignmentSummary = PacBio::Poa::PoaAlignmentSummary;

struct Subread
{
    char *seq;
    int flags;
};

std::string PoaDraftGenCore(const Subread *reads,
                            size_t num_sbr,
                            std::vector<SparsePoa::ReadKey> *readKeys,
                            std::vector<PoaAlignmentSummary> *summaries,
                            const size_t maxPoaCov,
                            const PoaSetting *settings)
{
    SparsePoa poa;
    size_t cov = 0;

    readKeys->clear();

    // readKeys->resize(sorted.size());
    std::string reads_0(reads->seq);

    for (int i = 0; i < num_sbr; ++i)
    {

        SparsePoa::ReadKey key = -1;
        std::string seq((reads + i)->seq);
        // std::cerr << "processing: " << seq << std::endl;
        if (settings->ed_unify_strand)
        {
            if (i == 0)
            {
                key = poa.AddRead(seq, *settings, true);
            }
            else
            {
                std::string rcSeq = PacBio::Data::ReverseComplement(seq);
                int distance1 =
                    Geneus::Align::EdlibGlobalAlignDistance(seq, reads_0);

                int distance2 = Geneus::Align::EdlibGlobalAlignDistance(rcSeq, reads_0);
                if (distance1 > distance2)
                {
                    key = poa.AddRead(rcSeq, *settings, true);
                }
                else
                {
                    key = poa.AddRead(seq, *settings, false);
                }
            }
        }
        else
        {
            key = poa.OrientAndAddRead(seq, *settings);
        }

        readKeys->emplace_back(key);
        if (key >= 0)
        {
            if ((++cov) >= maxPoaCov)
                break;
        }
    }

    // at least 50% of the reads should cover
    // TODO(lhepler) revisit this minimum coverage equation
    const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
    return poa.FindConsensus(minCov, *settings, &(*summaries))->Sequence;
}

int main()
{
    std::vector<SparsePoa::ReadKey> readKeys;
    std::vector<PoaAlignmentSummary> summaries;

    // const Subread *reads;
    Subread reads[] = {
        Subread{"ACGTACGTACGTACGT", 3},
        Subread{"ACGTACGTACGTACGT", 3},
        Subread{"ACGTACGTACGTACGT", 3},
    };

    PoaSetting setting;

    std::string result = PoaDraftGenCore(reads, 3, &readKeys, &summaries, 1000000, &setting);
    std::cout << "ConsensusSeq = " << result << std::endl;
    // SparsePoa poa;
    // size_t cov = 0;
}