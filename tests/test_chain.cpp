#include "chain/minimizer.hpp"
#include "chain/fenwick.hpp"
#include "chain/chainer.hpp"
#include "chain/paf.hpp"
#include <iostream>
#include <cassert>
#include <string>
#include <random>

void test_minimizer_extraction() {
    std::string seq = "ATGCGATCGATCGATCGATCGATCGATCGAACGTACGTTAGCTAGCT";
    uint32_t k = 7;
    uint32_t w = 5;

    auto mins = chain::extract_minimizers(seq, k, w);
    assert(!mins.empty());

    // Minimizers must have valid positions within bounds
    for (const auto& m : mins) {
        assert(m.pos + k <= seq.size());
        assert(m.span == k);
    }

    std::cout << "[PASS] test_minimizer_extraction (" << mins.size() << " minimizers extracted)\n";
}

void test_fenwick_tree() {
    chain::FenwickRMQ rmq(10);
    rmq.update(3, 100, 1);
    rmq.update(7, 250, 2);
    rmq.update(5, 180, 3);

    assert(rmq.query(2).anchor_idx == -1);
    assert(rmq.query(4).score == 100 && rmq.query(4).anchor_idx == 1);
    assert(rmq.query(6).score == 180 && rmq.query(6).anchor_idx == 3);
    assert(rmq.query(8).score == 250 && rmq.query(8).anchor_idx == 2);

    std::cout << "[PASS] test_fenwick_tree\n";
}

void test_forward_and_reverse_chaining() {
    std::string ref = "ACGTGTCGATCGATAGCTACGATCGATCGTACGATCGATCGAATCGATCGATCGATCG";
    std::string query_fwd = ref.substr(5, 45);
    std::string query_rev = chain::reverse_complement(query_fwd);

    chain::ChainParams params;
    params.k = 8;
    params.w = 4;
    chain::Chainer chainer(params);

    // Forward alignment
    auto fwd_chains = chainer.align(ref, query_fwd);
    assert(!fwd_chains.empty());
    assert(!fwd_chains[0].is_rev);
    assert(fwd_chains[0].num_anchors >= 3);

    // Reverse alignment
    auto rev_chains = chainer.align(ref, query_rev);
    assert(!rev_chains.empty());
    assert(rev_chains[0].is_rev);
    assert(rev_chains[0].num_anchors >= 3);

    std::cout << "[PASS] test_forward_and_reverse_chaining\n";
}

void test_circular_origin_wrapping() {
    // 500 bp circular genome using Mersenne Twister
    std::mt19937 rng(12345);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string ref(500, 'A');
    for (size_t i = 0; i < ref.size(); ++i) {
        ref[i] = bases[rng() % 4];
    }

    // Query spans the circular origin: last 80 bp of ref + first 80 bp of ref
    std::string query = ref.substr(420, 80) + ref.substr(0, 80);

    chain::ChainParams params;
    params.k = 11;
    params.w = 5;
    params.circular = true;

    chain::Chainer chainer(params);
    auto chains = chainer.align(ref, query);
    assert(!chains.empty());
    // Best chain must wrap around circular origin
    bool found_wrap = false;
    for (const auto& c : chains) {
        if (c.is_circular_wrap) {
            found_wrap = true;
            assert(c.num_anchors >= 5);
            break;
        }
    }
    assert(found_wrap);

    std::cout << "[PASS] test_circular_origin_wrapping\n";
}

void test_fenwick_vs_exact_dp() {
    std::mt19937 rng(42);
    const char bases[] = {'A', 'C', 'G', 'T'};

    std::string ref(1000, 'A');
    for (size_t i = 0; i < ref.size(); ++i) {
        ref[i] = bases[rng() % 4];
    }
    std::string query = ref.substr(50, 400);

    chain::ChainParams p_dp;
    p_dp.k = 11;
    p_dp.w = 5;
    p_dp.use_fenwick = false;

    chain::ChainParams p_fen;
    p_fen.k = 11;
    p_fen.w = 5;
    p_fen.use_fenwick = true;

    chain::Chainer c_dp(p_dp);
    chain::Chainer c_fen(p_fen);

    auto chains_dp = c_dp.align(ref, query);
    auto chains_fen = c_fen.align(ref, query);

    assert(!chains_dp.empty() && !chains_fen.empty());
    assert(chains_dp[0].score > 0);
    assert(chains_fen[0].score > 0);
    assert(chains_dp[0].score == chains_fen[0].score);

    std::cout << "[PASS] test_fenwick_vs_exact_dp (Score: " << chains_dp[0].score << ")\n";
}

int main() {
    std::cout << "Running chain unit tests...\n";
    test_minimizer_extraction();
    test_fenwick_tree();
    test_forward_and_reverse_chaining();
    test_circular_origin_wrapping();
    test_fenwick_vs_exact_dp();
    std::cout << "All chain unit tests passed successfully!\n";
    return 0;
}
