// Author: Lance Hepler

#include <Settings.h>
#include <PoaDraft.h>
#include <Sequence.h>
#include <PoaConsensus.h>
#include <SparsePoa.h>
#include <EdlibGlobalAlign.h>

#include <algorithm>
#include <boost/optional.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using SparsePoa = PacBio::Poa::SparsePoa;
using PoaAlignmentSummary = PacBio::Poa::PoaAlignmentSummary;
// using ConsensusSettings = PacBio::CCS::ConsensusSettings;

enum LocalContextFlags : uint8_t
{
    NO_LOCAL_CONTEXT = 0, ///< No context information available
    ADAPTER_BEFORE = 1,   ///< Adapter precedes subread
    ADAPTER_AFTER = 2,    ///< Adapter follows subread
    BARCODE_BEFORE = 4,   ///< Barcode precedes subread
    BARCODE_AFTER = 8,    ///< Barcode follows subread
    FORWARD_PASS = 16,    ///< Subread's orientation is 'forward pass'
    REVERSE_PASS = 32     ///< Subread's orientation is 'reverse pass'
};

extern "C"
{

    size_t _ComputeNPassesForNoPolish(const Subread *reads,
                                      size_t num_sbr,
                                      std::vector<SparsePoa::ReadKey> &readKeys,
                                      std::vector<PoaAlignmentSummary> &summaries, float minIdentity)
    {

        const size_t nReads = num_sbr;
        size_t nPasses = 0;
        for (int i = 0; i < nReads; i++)
        {
            if (readKeys[i] < 0)
                continue;
            if (summaries[readKeys[i]].AlignmentIdentity < minIdentity)
            {
                continue;
            }

            int flags = (reads + i)->flags;
            if (!(flags & ADAPTER_BEFORE && flags & ADAPTER_AFTER))
                continue;
            nPasses += 1;
        }
        return nPasses;
    }

    std::string _PoaDraftGenCore(const Subread *reads,
                                 size_t num_sbr,
                                 std::vector<SparsePoa::ReadKey> *readKeys,
                                 std::vector<PoaAlignmentSummary> *summaries,
                                 const size_t maxPoaCov,
                                 const PoaSetting *settings)
    {
        SparsePoa poa(settings->version);
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

    std::string _PoaDraftGenCoreWithAllFwdStrandCore(const Subread *reads,
                                                     size_t num_sbr,
                                                     std::vector<SparsePoa::ReadKey> *readKeys,
                                                     std::vector<PoaAlignmentSummary> *summaries,
                                                     const size_t maxPoaCov,
                                                     const PoaSetting *settings)
    {
        SparsePoa poa;
        size_t cov = 0;

        readKeys->clear();

        std::string reads_0(reads->seq);

        for (int i = 0; i < num_sbr; ++i)
        {
            std::string seq((reads + i)->seq);
            SparsePoa::ReadKey key = poa.AddRead(seq, *settings, false);
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

    void PoaDraftFreeResult(Result result)
    {
        free(result.seq);
    }

    Result PoaDraftGen(const Subread *reads, size_t num_sbr, const PoaSetting *setting)
    {
        std::vector<SparsePoa::ReadKey> readKeys;
        std::vector<PoaAlignmentSummary> summaries;
        std::string consensus_seq = _PoaDraftGenCore(reads, num_sbr, &readKeys, &summaries, 1000000, setting);
        // std::cerr << "res_seq: " << consensus_seq << std::endl;
        size_t n_passes = _ComputeNPassesForNoPolish(reads, num_sbr, readKeys, summaries, setting->min_identity);

        // char *consensus = new char[consensus_seq.size() + 1];
        char *consensus = (char *)malloc(sizeof(char) * (consensus_seq.size() + 1));
        std::strcpy(consensus, consensus_seq.c_str());
        Result result;
        result.n_passes = n_passes;
        result.seq = consensus;

        return result;
    }

    Result PoaDraftGenWithAllFwdStrand(const Subread *reads, size_t num_sbr, const PoaSetting *setting)
    {
        std::vector<SparsePoa::ReadKey> readKeys;
        std::vector<PoaAlignmentSummary> summaries;
        std::string consensus_seq = _PoaDraftGenCoreWithAllFwdStrandCore(reads, num_sbr, &readKeys, &summaries, 1000000, setting);
        // std::cerr << "res_seq: " << consensus_seq << std::endl;
        size_t n_passes = _ComputeNPassesForNoPolish(reads, num_sbr, readKeys, summaries, setting->min_identity);

        // char *consensus = new char[consensus_seq.size() + 1];
        char *consensus = (char *)malloc(sizeof(char) * (consensus_seq.size() + 1));
        std::strcpy(consensus, consensus_seq.c_str());
        Result result;
        result.n_passes = n_passes;
        result.seq = consensus;

        return result;
    }
}