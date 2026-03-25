// Author: Armin Töpfer

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "PoaGraph.h"
#include "RangeFinder.h"

namespace PacBio
{
    namespace Poa
    {
        namespace detail
        {

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
            class SdpRangeFinderV1 : public SdpRangeFinder
            {
            public:
                virtual ~SdpRangeFinderV1();

                SdpAnchorVector FindAnchors(const std::string &consensusSequence,
                                            const std::string &readSequence) const override;
            };

        } // namespace detail
    } // namespace Poa
} // namespace PacBio
