#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace spline_utils {

// A CDF coordinate similar to RS
template <class KeyType>
struct Coord {
    KeyType x;
    double y;
};

enum Orientation { Collinear, CW, CCW };
static constexpr double precision = std::numeric_limits<double>::epsilon();

// Compute orientation of three points (from RS)
static Orientation ComputeOrientation(const double dx1, const double dy1,
                                      const double dx2, const double dy2) {
    const double expr = std::fma(dy1, dx2, -std::fma(dy2, dx1, 0));
    if (expr > precision)
        return Orientation::CW;
    else if (expr < -precision)
        return Orientation::CCW;
    return Orientation::Collinear;
}

// Generic segment structure for the utility
template <typename KeyType, typename ValueType>
struct alignas(64) SplineSegment {
    KeyType start;
    KeyType end;
    double slope = 0.0;
    double intercept = 0.0;
    size_t offset = 0;
    bool is_line;
    KeyType step;
    size_t num_keys_covered = 0;
    std::vector<std::pair<KeyType, ValueType>> data;
};

// GreedySplineCorridor implementation adapted from RS
template <typename KeyType, typename ValueType>
class GreedySplineCorridorBuilder {
   public:
    struct SplinePoint {
        KeyType key;
        double position;
        ValueType value;
    };

    GreedySplineCorridorBuilder(size_t max_error = 32)
        : max_error_(max_error) {}

    std::vector<SplinePoint> BuildSplinePoints(
        const std::vector<std::pair<KeyType, ValueType>>& data) {
        if (data.empty()) return {};

        spline_points_.clear();
        curr_num_distinct_keys_ = 0;
        data_ref_ = &data;

        for (size_t i = 0; i < data.size(); ++i) {
            PossiblyAddPointToSpline(data[i].first, static_cast<double>(i),
                                     data[i].second);
        }

        // Ensure the last point is added
        if (!data.empty() && (spline_points_.empty() ||
                              spline_points_.back().key != data.back().first)) {
            AddPointToSpline(data.back().first,
                             static_cast<double>(data.size() - 1),
                             data.back().second);
        }

        return spline_points_;
    }

   private:
    void AddPointToSpline(KeyType key, double position,
                          const ValueType& value) {
        spline_points_.push_back({key, position, value});
    }

    void SetUpperLimit(KeyType key, double position) {
        upper_limit_ = {key, position};
    }

    void SetLowerLimit(KeyType key, double position) {
        lower_limit_ = {key, position};
    }

    void RememberPreviousCDFPoint(KeyType key, double position) {
        prev_point_ = {key, position};
    }

    void PossiblyAddPointToSpline(KeyType key, double position,
                                  const ValueType& value) {
        if (curr_num_distinct_keys_ == 0) {
            AddPointToSpline(key, position, value);
            ++curr_num_distinct_keys_;
            RememberPreviousCDFPoint(key, position);
            return;
        }

        if (curr_num_distinct_keys_ > 0 && key == spline_points_.back().key) {
            return;
        }

        ++curr_num_distinct_keys_;

        if (curr_num_distinct_keys_ == 2) {
            SetUpperLimit(key, position + max_error_);
            SetLowerLimit(key,
                          (position < max_error_) ? 0 : position - max_error_);
            RememberPreviousCDFPoint(key, position);
            return;
        }

        const auto& last = spline_points_.back();
        const double upper_y = position + max_error_;
        const double lower_y =
            (position < max_error_) ? 0 : position - max_error_;

        assert(upper_limit_.x >= last.key);
        assert(lower_limit_.x >= last.key);
        assert(key >= last.key);
        const double upper_limit_x_diff = upper_limit_.x - last.key;
        const double lower_limit_x_diff = lower_limit_.x - last.key;
        const double x_diff = key - last.key;

        assert(upper_limit_.y >= last.position);
        assert(position >= last.position);
        const double upper_limit_y_diff = upper_limit_.y - last.position;
        const double lower_limit_y_diff = lower_limit_.y - last.position;
        const double y_diff = position - last.position;

        assert(prev_point_.x != last.key);

        if ((ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                                y_diff) != Orientation::CW) ||
            (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                                y_diff) != Orientation::CCW)) {
            ValueType prev_value{};
            for (const auto& pair : *data_ref_) {
                if (pair.first == prev_point_.x) {
                    prev_value = pair.second;
                    break;
                }
            }

            AddPointToSpline(prev_point_.x, prev_point_.y, prev_value);
            SetUpperLimit(key, upper_y);
            SetLowerLimit(key, lower_y);
        } else {
            assert(upper_y >= last.position);
            const double upper_y_diff = upper_y - last.position;
            if (ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff,
                                   x_diff, upper_y_diff) == Orientation::CW) {
                SetUpperLimit(key, upper_y);
            }

            const double lower_y_diff = lower_y - last.position;
            if (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff,
                                   x_diff, lower_y_diff) == Orientation::CCW) {
                SetLowerLimit(key, lower_y);
            }
        }

        RememberPreviousCDFPoint(key, position);
    }

    size_t max_error_;
    std::vector<SplinePoint> spline_points_;
    size_t curr_num_distinct_keys_ = 0;
    const std::vector<std::pair<KeyType, ValueType>>* data_ref_ = nullptr;

    Coord<KeyType> upper_limit_;
    Coord<KeyType> lower_limit_;
    Coord<KeyType> prev_point_;
};

// Main function to convert non-line segments to spline segments using
// GreedySplineCorridor
template <typename KeyType, typename ValueType>
std::vector<SplineSegment<KeyType, ValueType>> CreateSplineSegments(
    const std::vector<std::pair<KeyType, ValueType>>& data,
    size_t max_error = 32) {
    if (data.empty()) {
        return {};
    }

    GreedySplineCorridorBuilder<KeyType, ValueType> builder(max_error);
    auto spline_points = builder.BuildSplinePoints(data);

    if (spline_points.size() <= 1) {
        SplineSegment<KeyType, ValueType> segment;
        segment.start = data.front().first;
        segment.end = data.back().first;
        segment.is_line = false;
        segment.data = data;
        segment.num_keys_covered = data.size();
        return {segment};
    }

    std::vector<SplineSegment<KeyType, ValueType>> segments;
    size_t cumulative_offset = 0;  // Concecutive spline segments, the estimated
                                   // position must remove this offset

    for (size_t i = 0; i < spline_points.size() - 1; ++i) {
        KeyType segment_start = spline_points[i].key;
        KeyType segment_end = spline_points[i + 1].key;

        auto start_it = std::lower_bound(
            data.begin(), data.end(),
            std::make_pair(segment_start, ValueType{}),
            [](const auto& a, const auto& b) { return a.first < b.first; });

        auto end_it = std::upper_bound(
            data.begin(), data.end(), std::make_pair(segment_end, ValueType{}),
            [](const auto& a, const auto& b) { return a.first < b.first; });

        if (start_it != end_it) {
            SplineSegment<KeyType, ValueType> segment;
            segment.start = start_it->first;
            segment.end = (end_it - 1)->first;
            segment.is_line = false;
            segment.data.assign(start_it, end_it);
            segment.num_keys_covered = segment.data.size();

            segment.offset = cumulative_offset;
            cumulative_offset += segment.data.size();

            if (segment.data.size() >= 2) {
                double start_pos = spline_points[i].position;
                double end_pos = spline_points[i + 1].position;
                double key_diff = segment_end - segment_start;
                if (key_diff > 0) {
                    segment.slope = (end_pos - start_pos) / key_diff;
                    segment.intercept =
                        start_pos - segment.slope * segment_start;
                }
            }

            segments.push_back(std::move(segment));
        }
    }

    return segments;
}

}  // namespace spline_utils
