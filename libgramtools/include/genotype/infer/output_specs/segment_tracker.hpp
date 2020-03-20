#ifndef GT_OUT_COMMON_HPP
#define GT_OUT_COMMON_HPP

#include <iostream>
#include <cassert>
#include <string>
#include <vector>

namespace gram::genotype {
    struct segment {
        std::string ID;
        std::size_t size;
    };

/**
 * Supports ID queries based on a position.
 * Model:
 *  - Stored segments are contiguous, with increasing positions
 *  - You can only query positions within a segment or in subsequent (=increasing) segments
 */
    class segmentTracker {
    private:
        std::vector<segment> segments;
        std::size_t min, max, global_max;
        std::size_t cur_idx;
    public:
        segmentTracker() = default;
        segmentTracker(std::istream &coords_file) : cur_idx(0), min(0), global_max(0) {
            segment next_segment{
                    "gramtools_prg",
                    std::numeric_limits<std::size_t>::max()
            };
            while (coords_file >> next_segment.ID >> next_segment.size) {
                global_max += next_segment.size;
                segments.push_back(next_segment);
            }
            // Case: pushes a single, maximally large segment
            if (segments.size() == 0) {
                segments.push_back(next_segment);
                global_max = std::numeric_limits<std::size_t>::max();
            }
            max = segments.at(0).size - 1;
        }

        std::string const &get_ID(std::size_t pos) {
            assert(pos >= min and pos < global_max);
            while (pos > max) {
                cur_idx++;
                min = max + 1;
                max += segments.at(cur_idx).size;
            }
            return segments.at(cur_idx).ID;
        }

        auto const& edge() const{ return max; }
        auto global_edge() const{ return global_max - 1; }
        void reset() {
            min = 0;
            cur_idx = 0;
            max = segments.at(0).size - 1;
        }
        auto const& get_segments() const {return segments;}
        auto num_segments() const {return segments.size();}
    };
}
#endif //GT_OUT_COMMON_HPP
