// Authors: David Alexander, Lance Hepler

#pragma once

#include <Interval.h>
#include <RangeFinder.h>
#include <Settings.h>

#include <memory>
#include <string>
#include <vector>

namespace PacBio
{

    // fwd decls
    namespace Consensus
    {
        class ScoredMutation;
    } // namespace Consensus

    namespace Poa
    {
        // fwd decls
        class PoaGraph;

        using Interval = PacBio::Data::Interval;
        using PoaGraph = PacBio::Poa::PoaGraph;

                //
        // Details of how a read aligns to the POA graph
        //
        struct PoaAlignmentSummary
        {
            bool ReverseComplementedRead;
            Interval ExtentOnRead;
            Interval ExtentOnConsensus;
            float AlignmentScore;
            float AlignmentIdentity;

            PoaAlignmentSummary()
                : ReverseComplementedRead{false}, ExtentOnRead(0, 0), ExtentOnConsensus(0, 0), AlignmentScore{0}, AlignmentIdentity{0}
            {
            }
        };

        struct PoaAlignmentOptions
        {
            bool ClipBegin;
            bool ClipEnd;
        };

        //
        // Partial order aligner with parsimonious memory usage
        //
        class SparsePoa
        {
        public:
            // ReadKey: Integer key representing
            //  >= 0: a read in the POA graph
            //    -1: read could not be inserted into POA graph
            using ReadKey = int;

        public:
            SparsePoa();
            SparsePoa(int version);
            ~SparsePoa();

            //
            // Add read, which must already have been oriented to be in the
            // "forward" convention
            //
            ReadKey AddRead(const std::string &readSequence, const PoaSetting &settings,
                            bool rc,
                            const PoaAlignmentOptions &alnOptions = PoaAlignmentOptions(),
                            float minScoreToAdd = 0);

            //
            // Find better orientation, (fwd or RC) and add as such
            //
            ReadKey OrientAndAddRead(const std::string &readSequence,
                                     const PoaSetting &settings,
                                     const PoaAlignmentOptions &alnOptions = PoaAlignmentOptions(),
                                     float minScoreToAdd = 0);

            //
            // Walk the POA and get the optimal consensus path
            //
            std::shared_ptr<const PacBio::Poa::PoaConsensus> FindConsensus(
                int minCoverage, const PoaSetting &settings,
                std::vector<PoaAlignmentSummary> *summaries = NULL) const;

            //
            // Serialize the POA graph to std::string
            //
            std::string ToGraphViz(int flags = 0, const PacBio::Poa::PoaConsensus *pc = nullptr) const;

            //
            // Serialize the POA graph to a file
            //
            void WriteGraphVizFile(const std::string &filename, int flags = 0,
                                   const PacBio::Poa::PoaConsensus *pc = nullptr) const;

            void WriteGraphCsvFile(const std::string &filename) const;

            //
            // Clean up the POA graph, pruning minority paths, to speed up
            // successive AddRead operations.
            //
            void PruneGraph(int minCoverage);

        private:
            void repCheck();

        private:
            using Path = std::vector<PoaGraph::Vertex>;

            PoaGraph *graph_;
            std::vector<Path> readPaths_;
            std::vector<bool> reverseComplemented_;
            PacBio::Poa::detail::SdpRangeFinder *rangeFinder_;
        };

    } // namespace Poa
} // namespace PacBio
