// Author: Lance Hepler

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "SparseAlignment.h"
#include "RangeFinder.h"
#include "RangeFinderV1.h"

#include "PoaGraphImpl.h"

#define WIDTH 30
#define DEBUG_RANGE_FINDER 0
#define KMER 6 

#if DEBUG_RANGE_FINDER
#include <iostream>
using std::cout;
using std::endl;
#endif // DEBUG_RANGE_FINDER

namespace PacBio
{
    namespace Poa
    {
        namespace detail
        {

            SdpRangeFinderV1::~SdpRangeFinderV1() = default;
          

            SdpAnchorVector SdpRangeFinderV1::FindAnchors(const std::string &consensusSequence,
                                                          const std::string &readSequence) const
            {
                return CCS::SparseAlign(KMER, consensusSequence, readSequence);
            }

        } // namespace detail
    } // namespace Poa
} // namespace PacBio
