#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>
#include <limits>

namespace chain {

struct MaxEntry {
    int64_t score = std::numeric_limits<int64_t>::min();
    int32_t anchor_idx = -1;
};

class FenwickRMQ {
public:
    explicit FenwickRMQ(size_t size) : size_(size), tree_(size + 1) {}

    void update(size_t idx, int64_t score, int32_t anchor_idx) {
        for (; idx <= size_; idx += idx & -idx) {
            if (score > tree_[idx].score) {
                tree_[idx].score = score;
                tree_[idx].anchor_idx = anchor_idx;
            }
        }
    }

    MaxEntry query(size_t idx) const {
        MaxEntry best;
        for (; idx > 0; idx -= idx & -idx) {
            if (tree_[idx].score > best.score) {
                best = tree_[idx];
            }
        }
        return best;
    }

    void clear() {
        std::fill(tree_.begin(), tree_.end(), MaxEntry{});
    }

private:
    size_t size_;
    std::vector<MaxEntry> tree_;
};

} // namespace chain
