#pragma once

#include "chainer.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

namespace chain {

inline std::string format_paf_record(const std::string& query_name,
                                    size_t query_len,
                                    const std::string& ref_name,
                                    size_t ref_len,
                                    const Chain& chain,
                                    bool is_primary = true) {
    std::ostringstream ss;
    int64_t q_span = std::abs(chain.query_end - chain.query_start);
    int64_t r_span = std::abs(chain.ref_end - chain.ref_start);
    int64_t block_len = std::max(q_span, r_span);
    int64_t matches = std::min(block_len, static_cast<int64_t>(chain.num_anchors * chain.anchors.front().span));

    uint8_t mapq = is_primary ? 60 : 0;

    ss << query_name << "\t"
       << query_len << "\t"
       << chain.query_start << "\t"
       << chain.query_end << "\t"
       << (chain.is_rev ? '-' : '+') << "\t"
       << ref_name << "\t"
       << ref_len << "\t"
       << chain.ref_start << "\t"
       << chain.ref_end << "\t"
       << matches << "\t"
       << block_len << "\t"
       << static_cast<int>(mapq) << "\t"
       << "tp:A:" << (is_primary ? "P" : "S") << "\t"
       << "cm:i:" << chain.num_anchors << "\t"
       << "s1:i:" << chain.score << "\t"
       << "circ:Z:" << (chain.is_circular_wrap ? "true" : "false");

    return ss.str();
}

} // namespace chain
