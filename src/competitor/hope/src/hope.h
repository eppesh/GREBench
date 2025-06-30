#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>  // std::setprecision(n), std::fixed
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

// #include "greedy_spline_corridor.h"
#include "hope_util.h"

namespace hopens {

template <typename KeyType, typename ValueType>
class Hope {
   public:
    struct alignas(64) Segment {
        KeyType start;
        KeyType end;
        double slope = 0.0;  // Used for quickly locating key in data
        double intercept = 0.0;
        size_t offset = 0;  // Used for concecutive spline segments
        std::unique_ptr<std::vector<Segment>> children;

        bool is_line;
        bool is_representation = false;  // distinguish representation segments
        KeyType step;                    // if it's line
        // Number of original keys this segment covers
        size_t num_keys_covered = 0;
        size_t segment_count = 1;  // total segments including children
        std::vector<std::pair<KeyType, ValueType>> data;

        // Hierarchical indexing for children
        std::vector<int> child_seg_index;
        // child_slope - used for quickly locating segment in children
        double child_slope = 0.0;
        double child_intercept = 0.0;

        Segment() : start(0), end(0), is_line(false), step(0) {}

        Segment(KeyType s, KeyType e, bool line = false)
            : start(s), end(e), is_line(line), step(0) {}

        // Move constructor
        Segment(Segment&& other) noexcept
            : start(other.start),
              end(other.end),
              is_line(other.is_line),
              is_representation(other.is_representation),
              data(std::move(other.data)),
              step(other.step),
              slope(other.slope),
              offset(other.offset),
              intercept(other.intercept),
              children(std::move(other.children)),
              num_keys_covered(other.num_keys_covered),
              segment_count(other.segment_count),
              child_seg_index(std::move(other.child_seg_index)),
              child_slope(other.child_slope),
              child_intercept(other.child_intercept) {}

        // Move assignment operator
        Segment& operator=(Segment&& other) noexcept {
            if (this != &other) {
                start = other.start;
                end = other.end;
                is_line = other.is_line;
                is_representation = other.is_representation;
                data = std::move(other.data);
                step = other.step;
                slope = other.slope;
                intercept = other.intercept;
                offset = other.offset;
                children = std::move(other.children);
                num_keys_covered = other.num_keys_covered;
                segment_count = other.segment_count;
                child_seg_index = std::move(other.child_seg_index);
                child_slope = other.child_slope;
                child_intercept = other.child_intercept;
            }
            return *this;
        }

        // Copy constructor (deep copy)
        Segment(const Segment& other)
            : start(other.start),
              end(other.end),
              is_line(other.is_line),
              is_representation(other.is_representation),
              data(other.data),
              step(other.step),
              slope(other.slope),
              intercept(other.intercept),
              offset(other.offset),
              num_keys_covered(other.num_keys_covered),
              segment_count(other.segment_count),
              child_seg_index(other.child_seg_index),
              child_slope(other.child_slope),
              child_intercept(other.child_intercept) {
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
                is_representation = other.is_representation;
                data = other.data;
                step = other.step;
                slope = other.slope;
                intercept = other.intercept;
                offset = other.offset;
                num_keys_covered = other.num_keys_covered;
                segment_count = other.segment_count;
                child_seg_index = other.child_seg_index;
                child_slope = other.child_slope;
                child_intercept = other.child_intercept;

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
   static constexpr size_t kRadixBits = 8;  // For radix table
   static constexpr size_t kRadixSize = 1 << kRadixBits;
   
   std::vector<Segment> root_segments_;
   size_t total_num_segments_;  // Number of segments in the whole index
   // Parameters
   size_t node_capacity_ = 100;
   double top_k_percentage_;
   size_t max_error_ = 32;
   size_t min_line_length_ = 10;
   size_t default_temp_node_capacity_ = 5;  // After retrain, 3->5
   double density_factor_ = 2.0;

    // Root-level indexing
    double root_slope_ = 0.0;
    double root_intercept_ = 0.0;
    std::vector<int> root_seg_index_;  // index for root segments
    // Alternative: Radix table for root indexing
    std::vector<uint32_t> root_radix_table_;
    bool use_radix_table_ = false;  // whether use radix table for root index

    // Key range for bounds checking
    KeyType min_key_;
    KeyType max_key_;
    bool key_range_initialized_ = false;

   public:
    // ====================================================================
    // PUBLIC API METHODS
    // ====================================================================

    Hope(size_t node_capacity = 100, double top_k = 0.05, size_t max_error = 32,
         size_t min_line_len = 10, bool use_radix = false)
        : node_capacity_(node_capacity),
          total_num_segments_(0),
          top_k_percentage_(top_k),
          max_error_(max_error),
          min_line_length_(min_line_len),
          use_radix_table_(use_radix) {}

    void SetParameters(size_t node_capacity, double top_k, size_t max_error,
                       size_t min_line_length, bool use_radix,
                       size_t temp_node_capacity) {
        node_capacity_ = node_capacity;
        top_k_percentage_ = top_k;
        max_error_ = max_error;
        min_line_length_ = min_line_length;
        use_radix_table_ = use_radix;
        default_temp_node_capacity_ = temp_node_capacity;
    }
    // Bulk loading
    void BulkLoad(const std::pair<KeyType, ValueType>* key_value, size_t num) {
        if (num == 0 || key_value == nullptr) return;

        // Create vector from array without copying - just wrap the data
        std::vector<std::pair<KeyType, ValueType>> data(key_value,
                                                        key_value + num);

        // Step 0: Sort and deduplicate
        SortAndDeduplicate(data);

        // Initialize key range
        if (!data.empty()) {
            min_key_ = data.front().first;
            max_key_ = data.back().first;
            key_range_initialized_ = true;
        }

        // Step 1: Find lines (arithmetic sequences)
        std::vector<Segment> segments = FindLinesAndCreateSegments(data);

        // Step 2: Build tree structure
        total_num_segments_ = segments.size();
        root_segments_ = BuildTreeLevel(std::move(segments), node_capacity_);

        // Update segment counts
        UpdateSegmentCounts(root_segments_);

        // Step 3: Build hierarchical indexing
        BuildHierarchicalIndex();

        // Debug Info
        std::cout << "[BulkLoad] min_key=" << min_key_
                  << ", max_key=" << max_key_
                  << "; total segments after bulk load:" << total_num_segments_
                  << "; root segments:" << root_segments_.size() << std::endl;
        PrintTree();
    }

    // Insert a single key-value pair
    bool Insert(KeyType key, const ValueType& value) {
        // Create point segment
        Segment point_seg;
        point_seg.start = key;
        point_seg.end = key;
        point_seg.is_line = false;
        point_seg.is_representation = false;
        point_seg.data.emplace_back(key, value);
        point_seg.num_keys_covered = 1;
        point_seg.segment_count = 1;

        // Insert and get the root segment that contains it
        int root_idx = -1;
        bool result =
            InsertPointSegment(root_segments_, point_seg, nullptr, root_idx);

        if (result && root_idx != -1) {
            // Check if retrain is needed
            CheckAndRetrain(root_idx);
        }

        return result;
    }

    // Lookup a key and return its value
    bool Lookup(KeyType key, ValueType& value) const {
        return LookupInSegments(root_segments_, key, value);
    }

    // Utility functions
    size_t GetRootSegmentCount() const { return root_segments_.size(); }
    size_t GetTotalSegmentCount() const { return total_num_segments_; }

    void PrintTree(int depth = 0) const {
        std::ofstream ofs("hope_segments_info.log");
        PrintSegments(root_segments_, depth, ofs);
        ofs.close();
    }

   private:
    // ====================================================================
    // LOOKUP METHODS
    // ====================================================================

    bool LookupInSegments(const std::vector<Segment>& segments, KeyType key,
                          ValueType& value) const {
        int seg_index = FindSegmentForLookup(segments, key);
        if (seg_index == -1) return false;

        const Segment& seg = segments[seg_index];

        // For representation segments, go directly to children
        if (seg.is_representation) {
            if (seg.children) {
                return LookupInSegments(*seg.children, key, value);
            }
            return false;
        }

        // For normal segments, check data first
        if (!seg.data.empty()) {
            if (LookupInSegmentData(seg, key, value)) {
                return true;
            }
        }

        // Then check children if any
        if (seg.children) {
            return LookupInSegments(*seg.children, key, value);
        }

        return false;
    }

    int FindSegmentForLookup(const std::vector<Segment>& segments,
                             KeyType key) const {
        if (&segments == &root_segments_) {
            if (use_radix_table_ && !root_radix_table_.empty()) {
                return SearchInRadixTable(key);
            } else if (!root_seg_index_.empty()) {
                return SearchInRootSegmentIndex(key);
            }
        }
        return BinarySearchSegment(segments, key);
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
        // For splines with only 1 key
        if (segment.data.size() == 1 && segment.data[0].first == key) {
            value = segment.data[0].second;
            return true;
        }

        // Use spline interpolation to estimate position
        double estimated_pos = EstimatePositionInSegment(segment, key);

        if (estimated_pos >= 0) {
            // Local search with error bounds using binary search
            return LocalSearchAroundEstimate(segment, key, value,
                                             estimated_pos);
        } else {
            // Fallback to binary search
            std::cout << "[LookupInSplineSegment] warning! key=" << key
                      << std::endl;
            return SearchInSegment(segment, key, value);
        }
    }

    bool LocalSearchAroundEstimate(const Segment& segment, KeyType key,
                                   ValueType& value,
                                   double estimated_pos) const {
        size_t estimate = static_cast<size_t>(std::round(estimated_pos));

        // Calculate search bounds
        size_t begin = (estimate < max_error_) ? 0 : (estimate - max_error_);
        size_t end = std::min(estimate + max_error_ + 1, segment.data.size());

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

    // ====================================================================
    // INSERT METHODS
    // ====================================================================

    bool InsertPointSegment(std::vector<Segment>& segments,
                            const Segment& point_seg, Segment* parent_segment,
                            int& root_segment_idx) {
        if (segments.empty()) {
            segments.push_back(point_seg);
            total_num_segments_++;
            if (parent_segment) {
                parent_segment->segment_count++;
            }
            // Rebuild index if this is root level
            if (&segments == &root_segments_) {
                root_segment_idx = 0;
                RebuildRootIndex();
            }
            return true;
        }

        KeyType key = point_seg.start;

        // Find insertion position
        int insert_pos = FindInsertPosition(segments, key);

        // Check if it falls on an existing segment
        if (insert_pos < segments.size() && key >= segments[insert_pos].start &&
            key <= segments[insert_pos].end) {
            // Insert as child of this segment
            if (&segments == &root_segments_) {
                root_segment_idx = insert_pos;
            }
            return InsertIntoSegment(segments[insert_pos], point_seg,
                                     root_segment_idx);
        }

        // Falls between segments
        bool is_full = (segments.size() >= node_capacity_);

        if (!is_full) {
            // Node not full, insert directly
            segments.insert(segments.begin() + insert_pos, point_seg);
            total_num_segments_++;
            if (parent_segment) {
                parent_segment->segment_count++;
            }

            if (&segments == &root_segments_) {
                root_segment_idx = insert_pos;
                RebuildRootIndex();
            }
            return true;
        } else {
            // Node is full, insert into previous segment
            int target_idx = (insert_pos > 0) ? insert_pos - 1 : 0;

            if (&segments == &root_segments_) {
                root_segment_idx = target_idx;
            }

            return InsertIntoSegment(segments[target_idx], point_seg,
                                     root_segment_idx);
        }
    }

    bool InsertIntoSegment(Segment& segment, const Segment& point_seg,
                           int& root_idx) {
        KeyType key = point_seg.start;

        // For non-representation segments, check if key exists in data
        if (!segment.is_representation && !segment.data.empty()) {
            size_t pos = SearchInData(segment.data, key);
            if (pos < segment.data.size() && segment.data[pos].first == key) {
                // Update existing key
                segment.data[pos].second = point_seg.data[0].second;
                return true;
            }
        }

        // Insert as child
        if (!segment.children) {
            segment.children = std::make_unique<std::vector<Segment>>();
        }

        return InsertPointSegment(*segment.children, point_seg, &segment,
                                  root_idx);
    }

    int FindInsertPosition(const std::vector<Segment>& segments,
                           KeyType key) const {
        // For non-root nodes, use binary search directly
        if (&segments != &root_segments_) {
            return BinarySearchInsertPosition(segments, key);
        }

        // For root node, use index if available
        if (use_radix_table_ && !root_radix_table_.empty()) {
            return SearchInRadixTable(key);
        } else if (!root_seg_index_.empty()) {
            return SearchInRootSegmentIndex(key);
        }

        return BinarySearchInsertPosition(segments, key);
    }

    int BinarySearchInsertPosition(const std::vector<Segment>& segments,
                                   KeyType key) const {
        int left = 0;
        int right = segments.size();

        while (left < right) {
            int mid = left + (right - left) / 2;
            if (segments[mid].start <= key) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // Adjust if key falls within previous segment
        if (left > 0 && key <= segments[left - 1].end) {
            return left - 1;
        }

        return left;
    }

    // ====================================================================
    // SEGMENT SEARCH AND INDEXING METHODS
    // ====================================================================

    int SearchInSegmentIndex(const std::vector<Segment>& segments,
                             const std::vector<int>& seg_index, double slope,
                             double intercept, KeyType key) const {
        if (segments.empty()) {
            return -1;
        }
        // Handle boundary cases
        if (key <= segments.front().start) {
            return 0;
        } else if (key >= segments.back().start) {
            return static_cast<int>(segments.size() - 1);
        }

        // Calculate estimated position with bounds checking
        int64_t position = static_cast<int64_t>(slope * key + intercept);
        if (position < 0 || position >= seg_index.size()) {
            // Fallback to binary search
            return BinarySearchSegment(segments, key);
        }

        int estimated_index = seg_index[position];

        // Bounds check for estimated_index
        if (estimated_index < 0 || estimated_index >= segments.size()) {
            return BinarySearchSegment(segments, key);
        }

        // Check immediate neighbors first (most common case)
        if (IsKeyInSegmentRange(segments, estimated_index, key)) {
            return estimated_index;
        }
        // Check adjacent segments
        if (estimated_index > 0 &&
            IsKeyInSegmentRange(segments, estimated_index - 1, key)) {
            return estimated_index - 1;
        }

        if (estimated_index + 1 < segments.size() &&
            IsKeyInSegmentRange(segments, estimated_index + 1, key)) {
            return estimated_index + 1;
        }

        // If not found in immediate vicinity, use exponential search + binary
        // search
        return ExponentialBinarySearch(segments, estimated_index, key);
    }

    // Check if key falls within a segment's range
    bool IsKeyInSegmentRange(const std::vector<Segment>& segments, int index,
                             KeyType key) const {
        if (index < 0 || index >= segments.size()) return false;

        KeyType seg_start = segments[index].start;
        KeyType seg_end = (index + 1 < segments.size())
                              ? segments[index + 1].start - 1
                              : segments[index].end;

        return key >= seg_start && key <= seg_end;
    }

    // Expoential search followed by binary search
    int ExponentialBinarySearch(const std::vector<Segment>& segments,
                                int estimated_index, KeyType key) const {
        int n = static_cast<int>(segments.size());

        if (segments[estimated_index].start < key) {
            // Search right
            int low = estimated_index + 2;
            int high = low;
            int step = 1;

            while (high < n && segments[high].start <= key) {
                low = high;
                step *= 2;
                high = std::min(low + step, n - 1);
            }
            return BinarySearchInRange(segments, low, high, key);
        } else {
            // Search left
            int high = estimated_index - 2;
            int low = high;
            int step = 1;

            while (low > 0 && segments[low].start >= key) {
                high = low;
                step *= 2;
                low = std::max(high - step, 0);
            }

            return BinarySearchInRange(segments, low, high, key);
        }
        return -1;
    }

    // Binary search in a specific range
    int BinarySearchInRange(const std::vector<Segment>& segments, int low,
                            int high, KeyType key) const {
        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (IsKeyInSegmentRange(segments, mid, key)) {
                return mid;
            }

            if (segments[mid].start <= key) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return -1;
    }

    // Fallback binary search for entire segments array
    int BinarySearchSegment(const std::vector<Segment>& segments,
                            KeyType key) const {
        int left = 0;
        int right = static_cast<int>(segments.size()) - 1;

        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (IsKeyInSegmentRange(segments, mid, key)) {
                return mid;
            }

            if (segments[mid].start <= key) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        return -1;
    }

    // Root-level segment search
    int SearchInRootSegmentIndex(KeyType key) const {
        return SearchInSegmentIndex(root_segments_, root_seg_index_,
                                    root_slope_, root_intercept_, key);
    }

    // ====================================================================
    // SEGMENT BUILDING AND TREE CONSTRUCTION METHODS
    // ====================================================================

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
            if (line_end - i + 1 >= min_line_length_) {
                Segment line_segment;
                line_segment.start = data[i].first;
                line_segment.end = data[line_end].first;
                line_segment.is_line = true;
                line_segment.num_keys_covered = line_end - i + 1;
                line_segment.segment_count = 1;

                // calculate step (common difference)
                if (line_end > i) {
                    line_segment.step = data[i + 1].first - data[i].first;
                }

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
                    max_error_);  // max_error = 32, similar to RadixSpline
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
                segment.offset = spline_seg.offset;
                segment.num_keys_covered = spline_seg.num_keys_covered;
                segment.segment_count = 1;

                /* CalculateLineParameters(spline_seg.data, 0,
                                        spline_seg.data.size() - 1, segment); */

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
        size_t segments_per_group =
            std::ceil(segments.size() / current_node_capacity);
        if (segments_per_group < 1) {
            segments_per_group = 1;
        }

        CreateGroupRepresentations(segments, is_top_segment, segments_per_group,
                                   current_node_capacity, result_segments);

        return result_segments;
    }

    std::vector<Segment> SelectTopSegments(const std::vector<Segment>& segments,
                                           size_t num_top) {
        // Sorted by how many keys the segment covers
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
            // Representation segments are not lines
            repr_segment.is_line = false;
            repr_segment.is_representation = true;  // Mark as representation
            repr_segment.segment_count = 0;

            // Calculate total coverage and collect child segments
            std::vector<Segment> child_segments;
            child_segments.reserve(group_indices.size());
            for (size_t idx : group_indices) {
                repr_segment.num_keys_covered += segments[idx].num_keys_covered;
                repr_segment.segment_count += segments[idx].segment_count;
                child_segments.emplace_back(segments[idx]);
            }

            // Recursively build child level if needed
            if (child_segments.size() > current_node_capacity) {
                child_segments = BuildTreeLevel(std::move(child_segments),
                                                current_node_capacity);
            }

            repr_segment.children = std::make_unique<std::vector<Segment>>(
                std::move(child_segments));
            repr_segment.segment_count++;  // Add self

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

    // ====================================================================
    // Root INDEXING METHODS
    // ====================================================================

    void BuildRadixTable() {
        if (root_segments_.empty()) return;
        size_t num_shift_bits = 18;
        uint32_t max_prefix = (max_key_ - min_key_) >> num_shift_bits;
        root_radix_table_.resize(max_prefix + 2, 0);
        size_t prev_prefix = 0;

        for (size_t i = 0; i < root_segments_.size(); ++i) {
            KeyType key = root_segments_[i].start;
            KeyType curr_prefix = (key - min_key_) >> num_shift_bits;

            if (curr_prefix != prev_prefix) {
                for (KeyType prefix = prev_prefix + 1; prefix <= curr_prefix;
                     ++prefix) {
                    root_radix_table_[prefix] = i;
                }
                prev_prefix = curr_prefix;
            }
        }

        // Finalize radix table
        ++prev_prefix;
        for (; prev_prefix < root_radix_table_.size(); ++prev_prefix) {
            root_radix_table_[prev_prefix] = root_segments_.size();
        }
    }

    int SearchInRadixTable(KeyType key) const {
        if (key < min_key_) return 0;
        if (key > max_key_) return root_segments_.size() - 1;

        KeyType prefix = (key - min_key_) >> 18;
        assert(prefix + 1 < root_radix_table_.size());
        uint32_t begin = root_radix_table_[prefix];
        uint32_t end = root_radix_table_[prefix + 1];
        if (end - begin < 32) {
            // Do linear search over narrowed range.
            uint32_t current = begin;
            while (root_segments_[current].start < key) ++current;
            return current;
        }

        // Do binary search over narrowed range.
        int pos = BinarySearchInRange(root_segments_, begin, end, key);
        return pos;
    }

    void RebuildRootIndex() {
        if (use_radix_table_) {
            BuildRadixTable();
        } else {
            BuildSegmentIndexForLevel(root_segments_, root_seg_index_,
                                      root_slope_, root_intercept_);
        }
    }
    // Build hierarchical indexing for all levels
    void BuildHierarchicalIndex() {
        /* BuildSegmentIndexForLevel(root_segments_, root_seg_index_,
           root_slope_, root_intercept_); */
        RebuildRootIndex();
        // BuildChildIndicesRecursively(root_segments_);
    }

    // Build segment index for a specific level
    void BuildSegmentIndexForLevel(const std::vector<Segment>& segments,
                                   std::vector<int>& seg_index, double& slope,
                                   double& intercept) {
        if (segments.empty()) {
            return;
        }
        size_t redundant_size = segments.size() * 90;  // 90x redundancy
        seg_index.resize(redundant_size, -1);

        KeyType start_key = segments.front().start;
        KeyType end_key = segments.back().start;

        if (end_key > start_key) {
            slope =
                static_cast<double>(redundant_size - 1) / (end_key - start_key);
            intercept = -slope * start_key;
        } else {
            slope = 0.0;
            intercept = 0.0;
        }

        // Fill the index array
        for (size_t i = 0; i < segments.size(); ++i) {
            int64_t position =
                static_cast<int64_t>(slope * (segments[i].start) + intercept);
            if (position >= 0 && position < redundant_size) {
                seg_index[position] = static_cast<int>(i);
            }
        }

        // Fill gaps with nearest valid indices
        int last_valid_index = 0;
        for (size_t i = 0; i < redundant_size; ++i) {
            if (seg_index[i] == -1) {
                seg_index[i] = last_valid_index;
            } else {
                last_valid_index = seg_index[i];
            }
        }
    }

    // Recursively build indices for all child segments
    /* void BuildChildIndicesRecursively(std::vector<Segment>& segments) {
        for (auto& segment : segments) {
            if (segment.children && !segment.children->empty()) {
                // For non-root nodes, we use binary search, so no index needed
                // But keep this for potential future optimization
                BuildSegmentIndexForLevel(
                    *segment.children, segment.child_seg_index,
                    segment.child_slope, segment.child_intercept);
                BuildChildIndicesRecursively(*segment.children);
            }
        }
    } */

    // ====================================================================
    // RETRAIN AND REBALANCE METHODS
    // ====================================================================

    void CheckAndRetrain(int root_idx) {
        if (root_idx < 0 || root_idx >= root_segments_.size()) return;

        const Segment& segment = root_segments_[root_idx];
        double density = static_cast<double>(segment.segment_count) *
                         node_capacity_ / total_num_segments_;

        if (density > density_factor_) {
            RetrainWithNeighbors(root_idx);
        }
    }

    void RetrainWithNeighbors(int root_idx) {
        // Select target segment and neighbors
        int start_idx = std::max(0, root_idx - 1);
        int end_idx =
            std::min(static_cast<int>(root_segments_.size() - 1), root_idx + 1);

        // Collect data from segments
        std::vector<std::pair<KeyType, ValueType>> retrain_data;
        CollectDataFromSegmentRange(root_segments_, start_idx, end_idx,
                                    retrain_data);

        // Determine temp node capacity (how many segments after retrain)
        size_t temp_node_capacity =
            default_temp_node_capacity_;  // 3->5 by default

        SortAndDeduplicate(retrain_data);

        // Retrain
        std::vector<Segment> new_segments =
            FindLinesAndCreateSegments(retrain_data);
        if (new_segments.size() > temp_node_capacity) {
            new_segments =
                BuildTreeLevel(std::move(new_segments), temp_node_capacity);
        }

        // Update segment counts and reset insert counts
        for (auto& seg : new_segments) {
            UpdateSegmentCount(seg);
        }

        // Calculate how many extra slots needed
        int old_slots = end_idx - start_idx + 1;
        int new_slots = new_segments.size();
        int extra_slots_needed = new_slots - old_slots;

        if (extra_slots_needed > 0) {
            // Need to compress other segments to make room
            CompressSegmentsForSpace(extra_slots_needed, start_idx, end_idx);
        }

        // Replace segments
        ReplaceSegmentRange(root_segments_, start_idx, end_idx,
                            std::move(new_segments));

        // Rebuild root index
        RebuildRootIndex();
    }

    void CompressSegmentsForSpace(int slots_needed, int exclude_start,
                                  int exclude_end) {
        // Find consecutive segments with lowest density to compress
        int best_start = -1;
        double min_density = std::numeric_limits<double>::max();
        int required_segments =
            slots_needed + 1;  // Need to compress N+1 to 1 to get N slots

        for (int i = 0; i <= root_segments_.size() - required_segments; ++i) {
            // Skip if overlaps with excluded range
            if ((i >= exclude_start && i <= exclude_end) ||
                (i + required_segments - 1 >= exclude_start &&
                 i + required_segments - 1 <= exclude_end)) {
                continue;
            }

            double total_density = 0;
            for (int j = i; j < i + required_segments; ++j) {
                total_density +=
                    static_cast<double>(root_segments_[j].segment_count) *
                    node_capacity_ / total_num_segments_;
            }

            if (total_density < min_density) {
                min_density = total_density;
                best_start = i;
            }
        }

        if (best_start != -1) {
            // Compress these segments into one
            std::vector<std::pair<KeyType, ValueType>> compress_data;
            CollectDataFromSegmentRange(root_segments_, best_start,
                                        best_start + required_segments - 1,
                                        compress_data);

            // Create single representation segment
            Segment repr_seg;
            repr_seg.start = root_segments_[best_start].start;
            repr_seg.end =
                root_segments_[best_start + required_segments - 1].end;
            repr_seg.is_representation = true;
            repr_seg.num_keys_covered = compress_data.size();

            SortAndDeduplicate(compress_data);

            // Build children
            std::vector<Segment> child_segments =
                FindLinesAndCreateSegments(compress_data);
            repr_seg.children = std::make_unique<std::vector<Segment>>(
                std::move(child_segments));
            UpdateSegmentCount(repr_seg);

            // Replace
            ReplaceSegmentRange(root_segments_, best_start,
                                best_start + required_segments - 1,
                                {std::move(repr_seg)});
        }
    }

    void CollectDataFromSegmentRange(
        const std::vector<Segment>& segments, size_t start_idx, size_t end_idx,
        std::vector<std::pair<KeyType, ValueType>>& data) {
        for (size_t i = start_idx; i <= end_idx; ++i) {
            CollectDataFromSegment(segments[i], data);
        }
    }

    void CollectDataFromSegment(
        const Segment& segment,
        std::vector<std::pair<KeyType, ValueType>>& data) {
        if (segment.children) {
            for (const auto& child : *(segment.children)) {
                CollectDataFromSegment(child, data);
            }
        } else {
            data.insert(data.end(), segment.data.begin(), segment.data.end());
        }
    }

    void ReplaceSegmentRange(std::vector<Segment>& segments, size_t start_idx,
                             size_t end_idx,
                             std::vector<Segment> new_segments) {
        size_t old_count = end_idx - start_idx + 1;
        size_t new_count = new_segments.size();

        // Update total segment count
        size_t old_segment_count = 0;
        for (size_t i = start_idx; i <= end_idx; ++i) {
            old_segment_count += segments[i].segment_count;
        }

        size_t new_segment_count = 0;
        for (const auto& seg : new_segments) {
            new_segment_count += seg.segment_count;
        }

        total_num_segments_ =
            total_num_segments_ - old_segment_count + new_segment_count;

        auto start_it = segments.begin() + start_idx;
        auto end_it = segments.begin() + end_idx + 1;
        segments.erase(start_it, end_it);
        segments.insert(segments.begin() + start_idx,
                        std::make_move_iterator(new_segments.begin()),
                        std::make_move_iterator(new_segments.end()));
    }

    // ====================================================================
    // UTILITY AND HELPER METHODS
    // ====================================================================

    void UpdateSegmentCounts(std::vector<Segment>& segments) {
        for (auto& seg : segments) {
            UpdateSegmentCount(seg);
        }
    }

    void UpdateSegmentCount(Segment& segment) {
        segment.segment_count = 1;  // Self
        if (segment.children) {
            UpdateSegmentCounts(*segment.children);
            for (const auto& child : *segment.children) {
                segment.segment_count += child.segment_count;
            }
        }
    }

    // Update key range when new keys are inserted
    /*     void UpdateKeyRange(KeyType key) {
            if (!key_range_initialized_) {
                min_key_ = key;
                max_key_ = key;
                key_range_initialized_ = true;
            } else {
                if (key < min_key_) {
                    min_key_ = key;
                }
                if (key > max_key_) {
                    max_key_ = key;
                }
            }
        } */

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

    double EstimatePositionInSegment(const Segment& segment,
                                     KeyType key) const {
        // Use slope and intercept: y = ax + b
        if (std::abs(segment.slope) > 1e-10) {
            double pos = segment.slope * key + segment.intercept;
            pos = pos - segment.offset;
            // Clamp to valid range
            return std::max(
                0.0,
                std::min(pos, static_cast<double>(segment.data.size() - 1)));
        }
        return -1.0;  // Invalid estimate
    }

    void PrintSegments(const std::vector<Segment>& segments, int depth,
                       std::ostream& out) const {
        std::string indent(depth * 2, ' ');

        for (const auto& seg : segments) {
            out << indent;

            if (seg.is_representation) {
                out << "Repr ";
            } else if (seg.is_line) {
                out << "Line ";
            } else {
                out << "Spline ";
            }

            out << "[" << static_cast<int>(seg.start) << ", "
                << static_cast<int>(seg.end) << "]";

            if (seg.is_line) {
                out << ", step: " << static_cast<int>(seg.step);
            }

            out << ", keys: " << seg.num_keys_covered
                << ", segments: " << seg.segment_count
                << ", slope: " << seg.slope << ", intercept: " << seg.intercept
                << ", offset: " << seg.offset;

            if (seg.children) {
                out << ", children: " << seg.children->size();
            }

            out << "\n";

            if (seg.children) {
                PrintSegments(*seg.children, depth + 1, out);
            }
        }
    }

};  // Class Hope

}  // namespace hopens