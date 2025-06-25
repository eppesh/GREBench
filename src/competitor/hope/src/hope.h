#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "greedy_spline_corridor.h"

namespace hopens {

template <typename KeyType, typename ValueType>
class Hope {
   public:
    struct Segment {
        KeyType start;
        KeyType end;
        bool is_line;
        std::vector<std::pair<KeyType, ValueType>> data;
        KeyType step;  // if it's line
        double slope = 0.0;
        double intercept = 0.0;

        // Tree structure fields
        std::unique_ptr<std::vector<Segment>> children;
        size_t num_keys_covered =
            0;  // Number of original keys this segment covers

        // Density tracking for retraining
        size_t insert_count =
            0;  // Number of insertions in this segment since last retrain

        Segment() : start(0), end(0), is_line(false), step(0) {}

        Segment(KeyType s, KeyType e, bool line = false)
            : start(s), end(e), is_line(line), step(0) {}

        // Move constructor
        Segment(Segment&& other) noexcept
            : start(other.start),
              end(other.end),
              is_line(other.is_line),
              data(std::move(other.data)),
              step(other.step),
              slope(other.slope),
              intercept(other.intercept),
              children(std::move(other.children)),
              num_keys_covered(other.num_keys_covered),
              insert_count(other.insert_count) {}

        // Move assignment operator
        Segment& operator=(Segment&& other) noexcept {
            if (this != &other) {
                start = other.start;
                end = other.end;
                is_line = other.is_line;
                data = std::move(other.data);
                step = other.step;
                slope = other.slope;
                intercept = other.intercept;
                children = std::move(other.children);
                num_keys_covered = other.num_keys_covered;
                insert_count = other.insert_count;
            }
            return *this;
        }

        // Copy constructor (deep copy)
        Segment(const Segment& other)
            : start(other.start),
              end(other.end),
              is_line(other.is_line),
              data(other.data),
              step(other.step),
              slope(other.slope),
              intercept(other.intercept),
              num_keys_covered(other.num_keys_covered),
              insert_count(other.insert_count) {
            if (other.children) {
                children =
                    std::make_unique<std::vector<Segment>>(*other.children);
            }
        }

        // Copy assignment operator
        Segment& operator=(const Segment& other) {
            if (this != &other) {
                start = other.start;
                end = other.end;
                is_line = other.is_line;
                data = other.data;
                step = other.step;
                slope = other.slope;
                intercept = other.intercept;
                num_keys_covered = other.num_keys_covered;
                insert_count = other.insert_count;

                if (other.children) {
                    children =
                        std::make_unique<std::vector<Segment>>(*other.children);
                } else {
                    children.reset();
                }
            }
            return *this;
        }
    };

   private:
    static constexpr size_t kMinLineLength = 5;
    static constexpr double kDensityFactorHigh = 2.0;
    static constexpr double kDensityFactorLow = 0.5;
    static constexpr size_t kMaxError = 32;

    std::vector<Segment> root_segments_;
    size_t node_capacity_;
    size_t total_num_segments_;  // Number of segments in the whole index
    double top_k_percentage_;

   public:
    Hope(size_t node_capacity = 100, double top_k = 0.05)
        : node_capacity_(node_capacity),
          total_num_segments_(0),
          top_k_percentage_(top_k) {}

    // Bulk loading
    void BulkLoad(const std::pair<KeyType, ValueType>* key_value, size_t num) {
        if (num == 0 || key_value == nullptr) return;

        // Create vector from array without copying - just wrap the data
        std::vector<std::pair<KeyType, ValueType>> data(key_value,
                                                        key_value + num);

        // Step 0: Sort and deduplicate
        SortAndDeduplicate(data);

        // Step 1: Find lines (arithmetic sequences)
        std::vector<Segment> segments = FindLinesAndCreateSegments(data);

        // Step 2: Build tree structure
        total_num_segments_ = segments.size();
        root_segments_ = BuildTreeLevel(std::move(segments), node_capacity_);
    }

    // Insert a single key-value pair
    bool Insert(KeyType key, const ValueType& value) {
        return InsertIntoSegments(root_segments_, key, value);
    }

    // Lookup a key and return its value
    bool Lookup(KeyType key, ValueType& value) const {
        return LookupInSegments(root_segments_, key, value);
    }

    // Utility functions
    size_t GetRootSegmentCount() const { return root_segments_.size(); }
    size_t GetTotalSegmentCount() const { return total_num_segments_; }

    void PrintTree(int depth = 0) const {
        PrintSegments(root_segments_, depth);
    }

   private:
    void SortAndDeduplicate(std::vector<std::pair<KeyType, ValueType>>& data) {
        // Sort by key
        std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        // Remove duplicates - keep first occurrence
        auto last = std::unique(
            data.begin(), data.end(),
            [](const auto& a, const auto& b) { return a.first == b.first; });
        data.erase(last, data.end());
    }
    std::vector<Segment> FindLinesAndCreateSegments(
        const std::vector<std::pair<KeyType, ValueType>>& data) {
        std::vector<Segment> segments;
        std::vector<bool> used(data.size(), false);

        // Step 1: find all lines
        FindLines(data, used, segments);

        // Step 2: create segments for non-line keys
        CreateNonLineSegments(data, used, segments);

        // Sort segments by start key
        std::sort(segments.begin(), segments.end(),
                  [](const Segment& a, const Segment& b) {
                      return a.start < b.start;
                  });
        return segments;
    }

    void FindLines(const std::vector<std::pair<KeyType, ValueType>>& data,
                   std::vector<bool>& used, std::vector<Segment>& segments) {
        for (size_t i = 0; i < data.size(); ++i) {
            if (used[i]) {
                continue;
            }

            // Try to find a line starting from position i
            size_t line_end = FindLongestLine(data, i);
            if (line_end - i + 1 >= kMinLineLength) {
                Segment line_segment;
                line_segment.start = data[i].first;
                line_segment.end = data[line_end].first;
                line_segment.is_line = true;
                line_segment.num_keys_covered = line_end - i + 1;

                // calculate step (common difference)
                if (line_end > i) {
                    line_segment.step = data[i + 1].first - data[i].first;
                }

                // Calculate slope and intercept for the line
                CalculateLineParameters(data, i, line_end, line_segment);

                // Store data keys
                for (size_t j = i; j <= line_end; ++j) {
                    line_segment.data.push_back(data[j]);
                    used[j] = true;
                }

                segments.push_back(std::move(line_segment));
            }
        }
    }

    size_t FindLongestLine(
        const std::vector<std::pair<KeyType, ValueType>>& data, size_t start) {
        if (start + 1 >= data.size()) {
            return start;
        }

        KeyType diff = data[start + 1].first - data[start].first;
        size_t end = start + 1;

        for (size_t i = start + 2; i < data.size(); ++i) {
            if (data[i].first - data[i - 1].first == diff) {
                end = i;
            } else {
                break;
            }
        }
        return end;
    }
    void CalculateLineParameters(
        const std::vector<std::pair<KeyType, ValueType>>& data, size_t start,
        size_t end, Segment& segment) {
        if (end <= start) {
            return;
        }

        // For lines (arithmetic sequences), we treat index as y-coordinate
        // x = key, y = position in the sequence (index)
        // For example: keys 2,4,6,8,10,12 -> points (2,0), (4,1), (6,2), (8,3),
        // (10,4), (12,5)

        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        size_t n = end - start + 1;

        for (size_t i = start; i <= end; ++i) {
            double x = static_cast<double>(data[i].first);  // key
            double y = static_cast<double>(
                i - start);  // index in sequence (0, 1, 2, ...)
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }

        // calculate the denominator
        double denom = n * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) > 1e-10) {
            segment.slope = (n * sum_xy - sum_x * sum_y) / denom;
            segment.intercept = (sum_y - segment.slope * sum_x) / n;
        }
    }

    void CreateNonLineSegments(
        const std::vector<std::pair<KeyType, ValueType>>& data,
        const std::vector<bool>& used, std::vector<Segment>& segments) {
        // Collect consecutive non-line keys
        for (size_t i = 0; i < data.size(); ++i) {
            if (used[i]) {
                continue;
            }

            size_t start = i;
            while (i < data.size() && !used[i]) {
                ++i;
            }
            --i;  // Adjust for the last increment

            // Collect data for this non-line range
            std::vector<std::pair<KeyType, ValueType>> non_line_data;
            for (size_t j = start; j <= i; ++j) {
                non_line_data.push_back(data[j]);
            }

            // Use GreedySplineCorridor to create optimized spline segments
            auto spline_segments =
                spline_utils::CreateSplineSegments<KeyType, ValueType>(
                    non_line_data,
                    kMaxError);  // max_error = 32, similar to RadixSpline
                                 // default

            // Convert spline_utils::SplineSegment to our Segment format
            for (const auto& spline_seg : spline_segments) {
                Segment segment;
                segment.start = spline_seg.start;
                segment.end = spline_seg.end;
                segment.is_line = spline_seg.is_line;
                segment.data = spline_seg.data;
                segment.step = spline_seg.step;
                segment.slope = spline_seg.slope;
                segment.intercept = spline_seg.intercept;
                segment.num_keys_covered = spline_seg.num_keys_covered;

                segments.push_back(std::move(segment));
            }
        }
    }

    std::vector<Segment> BuildTreeLevel(std::vector<Segment> segments,
                                        size_t current_node_capacity) {
        if (segments.size() <= current_node_capacity) {
            return segments;
        }

        // Step 1: select top-k segments to keep at this level
        size_t num_top_segments = static_cast<size_t>(
            std::max(1.0, top_k_percentage_ * current_node_capacity));

        std::vector<Segment> top_segments =
            SelectTopSegments(segments, num_top_segments);
        std::vector<bool> is_top_segment(segments.size(), false);

        // Mark top segments
        for (const auto& top_seg : top_segments) {
            for (size_t i = 0; i < segments.size(); ++i) {
                if (segments[i].start == top_seg.start &&
                    segments[i].end == top_seg.end) {
                    is_top_segment[i] = true;
                    break;
                }
            }
        }

        // Step 2: Group remaining segments and create representations
        std::vector<Segment> result_segments = std::move(top_segments);
        // We select a representation every 'segments_per_group' segments
        size_t segments_per_group = segments.size() / current_node_capacity;
        if (segments_per_group < 1) {
            segments_per_group = 1;
        }

        CreateGroupRepresentations(segments, is_top_segment, segments_per_group,
                                   current_node_capacity, result_segments);

        return result_segments;
    }

    std::vector<Segment> SelectTopSegments(const std::vector<Segment>& segments,
                                           size_t num_top) {
        std::vector<std::pair<size_t, size_t>>
            segment_coverage;  // (coverage, index)
        segment_coverage.reserve(segments.size());

        for (size_t i = 0; i < segments.size(); ++i) {
            segment_coverage.emplace_back(segments[i].num_keys_covered, i);
        }

        // Sort by coverage (descending)
        std::partial_sort(
            segment_coverage.begin(),
            segment_coverage.begin() +
                std::min(num_top, segment_coverage.size()),
            segment_coverage.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

        std::vector<Segment> top_segments;
        top_segments.reserve(std::min(num_top, segment_coverage.size()));
        for (size_t i = 0; i < std::min(num_top, segment_coverage.size());
             ++i) {
            top_segments.emplace_back(segments[segment_coverage[i].second]);
        }

        // Sort by start key
        std::sort(top_segments.begin(), top_segments.end(),
                  [](const Segment& a, const Segment& b) {
                      return a.start < b.start;
                  });

        return top_segments;
    }

    void CreateGroupRepresentations(const std::vector<Segment>& segments,
                                    const std::vector<bool>& is_top_segment,
                                    size_t segments_per_group,
                                    size_t current_node_capacity,
                                    std::vector<Segment>& result_segments) {
        if (segments.empty()) return;

        // Create groups while respecting top-k segment boundaries
        std::vector<std::vector<size_t>> groups;
        CreateGroupsWithTopSegmentBoundaries(segments, is_top_segment,
                                             segments_per_group, groups);

        // Create representation segments for each group
        for (const auto& group_indices : groups) {
            if (group_indices.empty()) continue;

            // Create representation segment
            Segment repr_segment;
            repr_segment.start = segments[group_indices.front()].start;
            repr_segment.end = segments[group_indices.back()].end;
            repr_segment.is_line =
                false;  // Representation segments are not lines

            // Calculate total coverage and collect child segments
            std::vector<Segment> child_segments;
            child_segments.reserve(group_indices.size());
            for (size_t idx : group_indices) {
                repr_segment.num_keys_covered += segments[idx].num_keys_covered;
                child_segments.emplace_back(segments[idx]);
            }

            // Recursively build child level if needed
            if (child_segments.size() > current_node_capacity) {
                child_segments = BuildTreeLevel(std::move(child_segments),
                                                current_node_capacity);
            }

            repr_segment.children = std::make_unique<std::vector<Segment>>(
                std::move(child_segments));

            result_segments.push_back(std::move(repr_segment));
        }

        // Sort final result by start key
        std::sort(result_segments.begin(), result_segments.end(),
                  [](const Segment& a, const Segment& b) {
                      return a.start < b.start;
                  });
    }

    // Mark the segments that should be grouped; top-seg may divide the groups
    // to avoid overlap
    void CreateGroupsWithTopSegmentBoundaries(
        const std::vector<Segment>& segments,
        const std::vector<bool>& is_top_segment, size_t segments_per_group,
        std::vector<std::vector<size_t>>& groups) {
        std::vector<size_t> current_group;
        size_t target_group_size = segments_per_group;

        for (size_t i = 0; i < segments.size(); ++i) {
            if (is_top_segment[i]) {
                // Found a top segment - finalize current group if not empty
                if (!current_group.empty()) {
                    groups.push_back(std::move(current_group));
                    current_group.clear();
                }
                // Top segments are already added to result_segments, skip them
                // here
                continue;
            }

            // Add non-top segment to current group
            current_group.push_back(i);

            // Check if we should finalize this group
            bool should_finalize = false;

            // finalize with three cases:
            // Case 1: Group has reached target size
            // Case 2: Next segment is a top segment (avoid splitting around it)
            // Case 3: This is the last segment
            should_finalize =
                (current_group.size() >= target_group_size) ||
                (i + 1 < segments.size() && is_top_segment[i + 1]) ||
                (i == segments.size() - 1);

            if (should_finalize && !current_group.empty()) {
                groups.push_back(std::move(current_group));
                current_group.clear();
            }
        }

        // Handle any remaining segments
        if (!current_group.empty()) {
            groups.push_back(std::move(current_group));
        }
    }

    bool InsertIntoSegments(std::vector<Segment>& segments, KeyType key,
                            const ValueType& value) {
        return true;
    }

   private:
    // Binary search for segment selection
    // Find the last segment whose start <= key; check if the key is within that
    // segment's [start,end] range
    typename std::vector<Segment>::const_iterator FindSegment(
        const std::vector<Segment>& segments, KeyType key) const {
        if (segments.empty()) return segments.end();

        size_t left = 0;
        size_t right = segments.size();

        // Find the rightmost segment where start <= key
        while (left < right) {
            size_t mid = left + (right - left) / 2;
            if (segments[mid].start <= key) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // Check if the previous segment contains the key
        if (left > 0) {
            size_t idx = left - 1;
            if (key >= segments[idx].start && key <= segments[idx].end) {
                return segments.begin() + idx;
            }
        }

        return segments.end();
    }

    typename std::vector<Segment>::iterator FindSegmentMutable(
        std::vector<Segment>& segments, KeyType key) {
        if (segments.empty()) return segments.end();

        size_t left = 0;
        size_t right = segments.size();

        // Find the rightmost segment where start <= key
        while (left < right) {
            size_t mid = left + (right - left) / 2;
            if (segments[mid].start <= key) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // Check if the previous segment contains the key
        if (left > 0) {
            size_t idx = left - 1;
            if (key >= segments[idx].start && key <= segments[idx].end) {
                return segments.begin() + idx;
            }
        }

        return segments.end();
    }

    // Find the position of the first key >= target_key using binary search
    size_t SearchInData(const std::vector<std::pair<KeyType, ValueType>>& data,
                        KeyType key) const {
        if (data.empty()) return data.size();

        size_t left = 0;
        size_t right = data.size();

        while (left < right) {
            size_t mid = left + (right - left) / 2;
            if (data[mid].first < key) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        // At this point, left == right, and points to the first index >= key
        return left;
    }

    // Custom binary search in range [begin, end)
    // Return the first index i in that range such that data[i].first >= key.
    size_t SearchInRange(const std::vector<std::pair<KeyType, ValueType>>& data,
                         KeyType key, size_t begin, size_t end) const {
        while (begin < end) {
            size_t mid = begin + (end - begin) / 2;
            if (data[mid].first < key) {
                begin = mid + 1;
            } else {
                end = mid;
            }
        }
        return begin;  // position of first element >= key
    }

    bool LookupInSegments(const std::vector<Segment>& segments, KeyType key,
                          ValueType& value) const {
        auto segment_it = FindSegment(segments, key);
        if (segment_it != segments.end()) {
            if (segment_it->children) {
                // This is a representation segment, search in children
                return LookupInSegments(*segment_it->children, key, value);
            } else {
                // This is a leaf segment, search in data
                return LookupInSegmentData(*segment_it, key, value);
            }
        }
        return false;  // Key not found in any segment
    }

    bool LookupInSegmentData(const Segment& segment, KeyType key,
                             ValueType& value) const {
        if (segment.is_line) {
            return LookupInLineSegment(segment, key, value);
        } else {
            return LookupInSplineSegment(segment, key, value);
        }
    }

    bool LookupInLineSegment(const Segment& segment, KeyType key,
                             ValueType& value) const {
        // For line segments, we can calculate the position directly
        if (segment.step > 0 && ((key - segment.start) % segment.step == 0)) {
            size_t position = (key - segment.start) / segment.step;
            if (position < segment.data.size()) {
                value = segment.data[position].second;
                return true;
            }
        }

        // Fallback: use spline interpolation + binary search
        std::cout << "[LookupInLineSegment] Warning! key=" << key
                  << ", segment.start=" << segment.start
                  << ", segment.end=" << segment.end << std::endl;
        return LookupInSplineSegment(segment, key, value);
    }

    bool LookupInSplineSegment(const Segment& segment, KeyType key,
                               ValueType& value) const {
        if (segment.data.empty()) {
            return false;
        }

        // Use spline interpolation to estimate position
        double estimated_pos = EstimatePosition(segment, key);

        if (estimated_pos >= 0) {
            // Local search with error bounds using binary search
            return LocalSearchAroundEstimate(segment, key, value,
                                             estimated_pos);
        } else {
            // Fallback to binary search
            return SearchInSegment(segment, key, value);
        }
    }

    bool LocalSearchAroundEstimate(const Segment& segment, KeyType key,
                                   ValueType& value,
                                   double estimated_pos) const {
        size_t estimate = static_cast<size_t>(std::round(estimated_pos));

        // Calculate search bounds
        size_t begin = (estimate < kMaxError) ? 0 : (estimate - kMaxError);
        size_t end = std::min(estimate + kMaxError + 1, segment.data.size());

        // First check the estimated position
        if (estimate < segment.data.size() &&
            segment.data[estimate].first == key) {
            value = segment.data[estimate].second;
            return true;
        }

        // binary search in the error window
        size_t pos = SearchInRange(segment.data, key, begin, end);

        if (pos < end && segment.data[pos].first == key) {
            value = segment.data[pos].second;
            return true;
        }

        return false;
    }

    bool SearchInSegment(const Segment& segment, KeyType key,
                         ValueType& value) const {
        size_t pos = SearchInData(segment.data, key);

        if (pos < segment.data.size() && segment.data[pos].first == key) {
            value = segment.data[pos].second;
            return true;
        }

        return false;
    }

    double EstimatePosition(const Segment& segment, KeyType key) const {
        // Use slope and intercept: y = ax + b
        if (std::abs(segment.slope) > 1e-10) {
            double pos = segment.slope * key + segment.intercept;
            // Clamp to valid range
            return std::max(
                0.0,
                std::min(pos, static_cast<double>(segment.data.size() - 1)));
        }
        return -1.0;  // Invalid estimate
    }

    void PrintSegments(const std::vector<Segment>& segments, int depth) const {
        std::string indent(depth * 2, ' ');

        for (const auto& seg : segments) {
            std::cout << indent;

            if (seg.is_line) {
                std::cout << "Line [" << static_cast<int>(seg.start) << ", "
                          << static_cast<int>(seg.end)
                          << "], step: " << static_cast<int>(seg.step)
                          << ", slope: " << std::fixed << std::setprecision(3)
                          << seg.slope << ", intercept: " << std::fixed
                          << std::setprecision(3) << seg.intercept
                          << ", keys: " << seg.num_keys_covered
                          << ", inserts: " << seg.insert_count << "\n";
            } else {
                std::cout << "Spline [" << static_cast<int>(seg.start) << ", "
                          << static_cast<int>(seg.end)
                          << "], keys: " << seg.num_keys_covered
                          << ", inserts: " << seg.insert_count << "\n";
            }

            if (seg.children) {
                PrintSegments(*seg.children, depth + 1);
            }
        }
    }

};  // Class Hope

}  // namespace hopens