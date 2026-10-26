#pragma once

#include "dna.hpp"
#include <cstdint>
#include <vector>
#include <deque>
#include <bit>
#include <optional>
#include <array>
#include <algorithm>

namespace chain {

inline constexpr std::array<uint64_t, 4> NTHASH_SEEDS = {
    0x3c8bf092713723bfULL,
    0x91849a6e2345ef01ULL,
    0xa4b1928424a1b029ULL,
    0x573b2a89364a2731ULL
};

inline constexpr uint64_t get_seed(uint8_t c) noexcept {
    return (c < 4) ? NTHASH_SEEDS[c] : 0ULL;
}

inline constexpr uint64_t get_comp_seed(uint8_t c) noexcept {
    return (c < 4) ? NTHASH_SEEDS[3 - c] : 0ULL;
}

struct Minimizer {
    uint64_t hash = 0;
    uint32_t pos = 0;      // 0-indexed start position in sequence
    uint32_t span = 0;     // k
    bool is_rev = false;   // true if reverse complement was chosen
};

inline std::vector<Minimizer> extract_minimizers(std::string_view seq, uint32_t k, uint32_t w) {
    std::vector<Minimizer> minimizers;
    if (seq.size() < k) return minimizers;

    // Collect all valid k-mers with their canonical hash and strand
    struct KmerInfo {
        uint64_t hash;
        uint32_t pos;
        bool is_rev;
    };
    std::vector<KmerInfo> kmers;
    kmers.reserve(seq.size() - k + 1);

    uint64_t hf = 0;
    uint64_t hr = 0;
    uint64_t hr_init = 0;
    std::vector<uint8_t> window;
    window.reserve(k);
    size_t window_head = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t code = base_to_2bit(seq[i]);
        if (code >= 4) {
            hf = 0; hr = 0; hr_init = 0;
            window.clear();
            window_head = 0;
            continue;
        }

        if (window.size() < k) {
            window.push_back(code);
            hf = std::rotl(hf, 1) ^ get_seed(code);
            hr_init ^= std::rotl(get_comp_seed(code), static_cast<int>(window.size() - 1));

            if (window.size() == k) {
                hr = hr_init;
                bool is_rev = (hr < hf);
                kmers.push_back({is_rev ? hr : hf, static_cast<uint32_t>(i + 1 - k), is_rev});
            }
        } else {
            uint8_t out_code = window[window_head];
            window[window_head] = code;
            window_head = (window_head + 1) % k;

            hf = std::rotl(hf, 1) ^ std::rotl(get_seed(out_code), static_cast<int>(k)) ^ get_seed(code);
            hr = std::rotr(hr, 1) ^ std::rotr(get_comp_seed(out_code), 1) ^ std::rotl(get_comp_seed(code), static_cast<int>(k - 1));

            bool is_rev = (hr < hf);
            kmers.push_back({is_rev ? hr : hf, static_cast<uint32_t>(i + 1 - k), is_rev});
        }
    }

    if (kmers.empty()) return minimizers;

    // Sliding window of size w over kmers using deque for monotonic min-queue
    std::deque<size_t> dq;
    uint32_t last_pos = UINT32_MAX;
    uint64_t last_hash = UINT64_MAX;

    for (size_t i = 0; i < kmers.size(); ++i) {
        while (!dq.empty() && kmers[dq.back()].hash >= kmers[i].hash) {
            dq.pop_back();
        }
        dq.push_back(i);

        // Remove elements outside sliding window of size w
        if (i >= w && dq.front() <= i - w) {
            dq.pop_front();
        }

        // When window has at least w elements (or at end of sequence), record minimum
        if (i + 1 >= w || i == kmers.size() - 1) {
            size_t min_idx = dq.front();
            const auto& m = kmers[min_idx];
            if (m.pos != last_pos || m.hash != last_hash) {
                minimizers.push_back({m.hash, m.pos, k, m.is_rev});
                last_pos = m.pos;
                last_hash = m.hash;
            }
        }
    }

    return minimizers;
}

} // namespace chain
