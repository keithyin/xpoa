#pragma once

extern "C"
{
    struct PoaSetting
    {
        float min_identity = 0.82;
        int match_score = 3;
        int mismatch_score = -5;
        int insertion_score = -2;
        int deletion_score = -2;
        int ed_unify_strand = 1;
    };
}