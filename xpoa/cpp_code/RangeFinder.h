// Author: Armin Töpfer

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <PoaGraph.h>
#include <algorithm>

namespace PacBio
{
    namespace Poa
    {
        namespace detail
        {

            // an Anchor represents a point (cssPos, readPos)
            typedef std::pair<size_t, size_t> SdpAnchor;
            typedef std::vector<SdpAnchor> SdpAnchorVector;

            class PoaGraphImpl;

            // using PacBio::Consensus::PoaGraph;

            //
            // SdpRangeFinder objects are responsible for identifying the range
            // of read positions that we should seek to align to a POA vertex;
            // this implementation uses SDP to identify fairly narrow bands,
            // enabling sparse memory usage.
            //
            // This is an abstract class that will be inherited in a client
            // library that has access to an SDP method.
            //
            // RangeFinder state goes away on next call to InitRangeFinder.  We could
            // have dealt with this using a factory pattern but bleh.
            //
            class SdpRangeFinder
            {
            protected:
                std::map<PoaGraph::Vertex, std::pair<int, int>> alignableReadIntervalByVertex_;

            public:
                virtual ~SdpRangeFinder();

                virtual void InitRangeFinder(const PoaGraphImpl &poaGraph,
                                             const std::vector<PoaGraph::Vertex> &consensusPath,
                                             const std::string &consensusSequence, const std::string &readSequence);

                // TODO: write contract
                virtual std::pair<int, int> FindAlignableRange(PoaGraph::Vertex v);

            protected:
                // TODO: write contract
                virtual SdpAnchorVector FindAnchors(const std::string &consensusSequence,
                                                    const std::string &readSequence) const = 0;
            };

            static inline bool compareAnchorsOnCssPos(const SdpAnchor &a1, const SdpAnchor &a2)
            {
                return a1.first < a2.first;
            }

            static inline const SdpAnchor *binarySearchAnchors(const SdpAnchorVector &anchors, size_t cssPosition)
            {
                auto found = std::lower_bound(anchors.begin(), anchors.end(), std::make_pair(cssPosition, -1),
                                              compareAnchorsOnCssPos);
                if (found != anchors.end() && (*found).first == cssPosition)
                {
                    return &(*found);
                }
                else
                {
                    return nullptr;
                }
            }

            using std::max;
            using std::min;

            using Interval = std::pair<int, int>;

            inline std::string formatIntervalEndpoint(int i)
            {
                if (i == INT_MAX / 2)
                {
                    return "inf";
                }
                else if (i == -INT_MAX / 2)
                {
                    return "-inf";
                }
                else
                {
                    return std::to_string(i);
                }
            }

            inline std::string formatInterval(const Interval &ival)
            {
                return (std::string("[") + formatIntervalEndpoint(ival.first) + std::string(", ") +
                        formatIntervalEndpoint(ival.second) + std::string(")"));
            }

            // Canonical empty interval....
            const Interval emptyInterval = Interval(INT_MAX / 2, -INT_MAX / 2);

            inline Interval RangeUnion(const Interval &range1, const Interval &range2)
            {
                return Interval(min(range1.first, range2.first), max(range1.second, range2.second));
            }

            inline Interval RangeUnion(const std::vector<Interval> &ranges)
            {
                Interval result = emptyInterval;
                for (const Interval &r : ranges)
                {
                    result = RangeUnion(result, r);
                }
                return result;
            }

            inline Interval next(const Interval &v, int upperBound)
            {
                if (v == emptyInterval)
                    return emptyInterval;
                else
                    return Interval(min(v.first + 1, upperBound), min(v.second + 1, upperBound));
            }

            inline Interval prev(const Interval &v, int lowerBound = 0)
            {
                if (v == emptyInterval)
                    return emptyInterval;
                else
                    return Interval(max(v.first - 1, lowerBound), max(v.second - 1, lowerBound));
            }

        } // namespace detail
    } // namespace Poa
} // namespace PacBio
