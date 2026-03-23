#pragma once

// Author: rocyin

#include "Settings.h"
#include "Sequence.h"
#include "PoaConsensus.h"
#include "SparsePoa.h"
#include "edlib/EdlibGlobalAlign.h"

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

extern "C"
{
    struct Subread
    {
        char *seq;
        int flags;
    };

    struct Result
    {
        char *seq;
        size_t n_passes;
    };

    void PoaDraftFreeResult(Result result);

    Result PoaDraftGen(const Subread *reads, size_t num_sbr, const PoaSetting *setting);

    Result PoaDraftGenWithAllFwdStrand(const Subread *reads, size_t num_sbr, const PoaSetting *setting);
}