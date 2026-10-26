#pragma once

#include "minimizer.hpp"
#include "fenwick.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace chain {

struct Anchor {
    int64_t ref_pos = 0;
    int64_t query_pos = 0;
    uint32_t span = 0;
    bool is_rev = false;
};

struct Chain {
    int64_t ref_start = 0;
    int64_t ref_end = 0;
    int64_t query_start = 0;
    int64_t query_end = 0;
    int64_t score = 0;
    uint32_t num_anchors = 0;
    bool is_rev = false;
    bool is_circular_wrap = false;
    std::vector<Anchor> anchors;
};

struct ChainParams {
    uint32_t k = 15;
    uint32_t w = 10;
    int64_t max_gap = 5000;
    int32_t bandwidth = 64;
    double gap_cost_diff = 0.05;
    uint32_t max_occ = 200;
    uint32_t min_anchors = 2;
    bool circular = false;
    bool use_fenwick = false;
};

inline int64_t compute_gap_penalty(int64_t dx, int64_t dy, double gap_cost_diff) {
    int64_t diff = std::abs(dx - dy);
    int64_t min_gap = std::min(dx, dy);
    // Gap cost: non-linear or affine-like
    return static_cast<int64_t>(std::round(gap_cost_diff * static_cast<double>(diff) + 0.01 * static_cast<double>(min_gap)));
}

class Chainer {
public:
    explicit Chainer(ChainParams params) : params_(params) {}

    std::vector<Chain> align(std::string_view ref_seq,
                            std::string_view query_seq) {
        int64_t ref_len = static_cast<int64_t>(ref_seq.size());
        int64_t query_len = static_cast<int64_t>(query_seq.size());

        auto ref_min = extract_minimizers(ref_seq, params_.k, params_.w);
        auto query_min = extract_minimizers(query_seq, params_.k, params_.w);

        // Index reference minimizers by hash
        struct RefSeed {
            int64_t pos;
            uint32_t span;
            bool is_rev;
        };
        std::unordered_map<uint64_t, std::vector<RefSeed>> ref_index;
        ref_index.reserve(ref_min.size() * (params_.circular ? 2 : 1));

        for (const auto& m : ref_min) {
            ref_index[m.hash].push_back({static_cast<int64_t>(m.pos), m.span, m.is_rev});
            if (params_.circular) {
                // Duplicate reference minimizer into [L_ref, 2*L_ref)
                ref_index[m.hash].push_back({static_cast<int64_t>(m.pos) + ref_len, m.span, m.is_rev});
            }
        }

        // Collect anchors separating forward and reverse strand
        std::vector<Anchor> fwd_anchors;
        std::vector<Anchor> rev_anchors;

        for (const auto& qm : query_min) {
            auto it = ref_index.find(qm.hash);
            if (it == ref_index.end() || it->second.size() > params_.max_occ) continue;

            for (const auto& rm : it->second) {
                bool is_rev = (qm.is_rev != rm.is_rev);
                int64_t q_pos = static_cast<int64_t>(qm.pos);
                if (is_rev) {
                    // For reverse strand, reflect query pos: q' = query_len - 1 - q_pos
                    int64_t mapped_q = query_len - static_cast<int64_t>(qm.span) - q_pos;
                    rev_anchors.push_back({rm.pos, mapped_q, qm.span, true});
                } else {
                    fwd_anchors.push_back({rm.pos, q_pos, qm.span, false});
                }
            }
        }

        std::vector<Chain> all_chains;

        auto chain_strand = [&](std::vector<Anchor>& anchors, bool is_rev) {
            if (anchors.empty()) return;

            std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
                if (a.ref_pos != b.ref_pos) return a.ref_pos < b.ref_pos;
                return a.query_pos < b.query_pos;
            });

            size_t n = anchors.size();
            std::vector<int64_t> f(n, 0);
            std::vector<int32_t> parent(n, -1);

            for (size_t i = 0; i < n; ++i) {
                f[i] = static_cast<int64_t>(anchors[i].span);
            }

            if (!params_.use_fenwick) {
                // Banded Dynamic Programming: O(N * B)
                for (size_t i = 0; i < n; ++i) {
                    size_t start_j = (i > static_cast<size_t>(params_.bandwidth)) ? (i - params_.bandwidth) : 0;
                    for (size_t j = start_j; j < i; ++j) {
                        int64_t dx = anchors[i].ref_pos - anchors[j].ref_pos;
                        int64_t dy = anchors[i].query_pos - anchors[j].query_pos;

                        if (dx > 0 && dy > 0 && dx <= params_.max_gap && dy <= params_.max_gap) {
                            int64_t penalty = compute_gap_penalty(dx, dy, params_.gap_cost_diff);
                            int64_t match_gain = std::min(static_cast<int64_t>(anchors[i].span), dy);
                            int64_t score = f[j] + match_gain - penalty;

                            if (score > f[i]) {
                                f[i] = score;
                                parent[i] = static_cast<int32_t>(j);
                            }
                        }
                    }
                }
            } else {
                // Coordinate-compressed Fenwick RMQ Tree
                std::vector<int64_t> q_coords;
                q_coords.reserve(n);
                for (const auto& a : anchors) q_coords.push_back(a.query_pos);
                std::sort(q_coords.begin(), q_coords.end());
                q_coords.erase(std::unique(q_coords.begin(), q_coords.end()), q_coords.end());

                FenwickRMQ fenwick(q_coords.size());

                for (size_t i = 0; i < n; ++i) {
                    auto it = std::lower_bound(q_coords.begin(), q_coords.end(), anchors[i].query_pos);
                    size_t rank = static_cast<size_t>(std::distance(q_coords.begin(), it));

                    if (rank > 0) {
                        auto best = fenwick.query(rank);
                        if (best.anchor_idx >= 0) {
                            int64_t dx = anchors[i].ref_pos - anchors[best.anchor_idx].ref_pos;
                            int64_t dy = anchors[i].query_pos - anchors[best.anchor_idx].query_pos;
                            if (dx <= params_.max_gap && dy <= params_.max_gap) {
                                int64_t penalty = compute_gap_penalty(dx, dy, params_.gap_cost_diff);
                                int64_t match_gain = std::min(static_cast<int64_t>(anchors[i].span), dy);
                                int64_t score = best.score + match_gain - penalty;
                                if (score > f[i]) {
                                    f[i] = score;
                                    parent[i] = best.anchor_idx;
                                }
                            }
                        }
                    }
                    fenwick.update(rank + 1, f[i], static_cast<int32_t>(i));
                }
            }

            // Extract disjoint chains by finding highest scoring unvisited end
            std::vector<bool> used(n, false);
            std::vector<size_t> order(n);
            for (size_t i = 0; i < n; ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return f[a] > f[b];
            });

            for (size_t idx : order) {
                if (used[idx] || f[idx] < static_cast<int64_t>(params_.k)) continue;

                Chain chain;
                chain.score = f[idx];
                chain.is_rev = is_rev;

                int32_t curr = static_cast<int32_t>(idx);
                std::vector<Anchor> chain_anchors;
                while (curr >= 0 && !used[curr]) {
                    used[curr] = true;
                    chain_anchors.push_back(anchors[curr]);
                    curr = parent[curr];
                }

                if (chain_anchors.empty()) continue;
                std::reverse(chain_anchors.begin(), chain_anchors.end());

                chain.num_anchors = static_cast<uint32_t>(chain_anchors.size());
                if (chain.num_anchors < params_.min_anchors) continue;
                chain.ref_start = chain_anchors.front().ref_pos;
                chain.ref_end = chain_anchors.back().ref_pos + chain_anchors.back().span;

                if (!is_rev) {
                    chain.query_start = chain_anchors.front().query_pos;
                    chain.query_end = chain_anchors.back().query_pos + chain_anchors.back().span;
                } else {
                    // Map back to original query coordinates
                    int64_t orig_start = query_len - (chain_anchors.back().query_pos + chain_anchors.back().span);
                    int64_t orig_end = query_len - chain_anchors.front().query_pos;
                    chain.query_start = orig_start;
                    chain.query_end = orig_end;
                }

                // Check circular wrap
                if (params_.circular && chain.ref_end >= ref_len) {
                    chain.is_circular_wrap = true;
                    // Normalize reference coordinates modulo ref_len
                    if (chain.ref_start >= ref_len) {
                        chain.ref_start -= ref_len;
                        chain.ref_end -= ref_len;
                    }
                }

                chain.anchors = std::move(chain_anchors);
                all_chains.push_back(std::move(chain));
            }
        };

        chain_strand(fwd_anchors, false);
        chain_strand(rev_anchors, true);

        std::sort(all_chains.begin(), all_chains.end(), [](const Chain& a, const Chain& b) {
            return a.score > b.score;
        });

        return all_chains;
    }

private:
    ChainParams params_;
};

} // namespace chain
