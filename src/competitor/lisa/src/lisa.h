

#ifndef __LISA_H__
#define __LISA_H__

#include <omp.h>
// #include <pthread.h>
#include <tbb/parallel_sort.h>

#include <atomic>

// #include "../logger.h"
#include "lisa_base.h"
#include "lisa_util.h"

namespace lisa {
template <class KeyType, class ValueType>
class LISA {
   public:
    using value_type = std::pair<KeyType, ValueType>;

    LISA(value_type min_key, value_type max_key, size_t num_radix_bits = 18,
         size_t max_error = 32, size_t alpha = 10)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(
              GetNumShiftBits(max_key.first - min_key.first, num_radix_bits)),
          max_error_(max_error),
          num_keys_(0),
          alpha_(alpha),
          hit_num_(0),
          avg_line_length_(0),
          hit_count_(0),
          miss_count_(0),
          lookup_count_(0) {}

    ~LISA() {
        // TODO: free memory
        std::cout << "deconstructor" << std::endl;
    }

    void SetParameters(value_type min_key, value_type max_key,
                       size_t num_radix_bits = 18, size_t max_error = 32,
                       size_t alpha = 10) {
        min_key_ = min_key;
        max_key_ = max_key;
        num_radix_bits_ = num_radix_bits;
        max_error_ = max_error;
        alpha_ = alpha;
        num_shift_bits_ =
            GetNumShiftBits(max_key.first - min_key.first, num_radix_bits);
    }

    // input: pointer of keys array and num of keys
    void bulk_load(const value_type* data, size_t num_keys) {
        if (num_keys <= 0) {
            return;
        }
        data_.clear();
        // Copy data
        for (size_t i = 0; i < num_keys; ++i) {
            data_.push_back(data[i]);
        }

        min_key_ = data[0];
        max_key_ = data[num_keys - 1];
        num_keys_ = num_keys;
        num_shift_bits_ =
            GetNumShiftBits(max_key_.first - min_key_.first, num_radix_bits_);

        radix_table_.clear();
        // Initialize radix table, needs to contain all prefixes up to the
        // largest
        const uint32_t max_prefix =
            (max_key_.first - min_key_.first) >> num_shift_bits_;
        radix_table_.resize(max_prefix + 2, 0);
        spline_points_.clear();

        build_model(data_, spline_points_, radix_table_);
        data_.clear();
    }

    bool Lookup(const KeyType& key) {
        // lookup_count_++;
        size_t index = 1;
        if (key > min_key_.first && key <= max_key_.first) {
            index = GetSplineSegment(key);
        } else if (key > max_key_.first) {
            index = spline_points_.size();  // don't need 'up'
        }

        Coord<KeyType>& down = spline_points_[index - 1];
        extra<KeyType, ValueType>* ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(down));
        if (needCheckBuffer(down)) {
            auto it = ext->buffer->find(key);
            if (it != ext->buffer->end()) {
                // check if deleted
                if (it->second != (KeyType)(-1)) {
                    return true;
                }
                return false;
            }
        }
        size_t estimate = 0;
        if (key >= min_key_.first && key <= max_key_.first) {
            Coord<KeyType>& up = spline_points_[index];
            // Compute slope.
            const double x_diff = up.x - down.x;
            const double y_diff = up.y - down.y;
            const double slope = y_diff / x_diff;

            // Interpolate.
            const double key_diff = key - down.x;
            estimate = std::fma(key_diff, slope, down.y);

            if ((LTYPE)getType(down) == kRsSeg) {
                estimate = estimate - (size_t)(down.rsseg_start_pos);
            } else {
                estimate = estimate - (size_t)(down.y);
            }
        } else {
            return false;
        }
        if ((*(ext->data))[estimate].first == key) {
            // hit_count_++;
            return true;
        }

        if ((LTYPE)getType(down) == kLine) {
            // if the searching key belongs to a line but not found, no need to
            // do following BS
            return false;
        }

        const size_t begin =
            (estimate < max_error_) ? 0 : (estimate - max_error_);
        const size_t end = (estimate + max_error_ + 2 > ext->data_size)
                               ? ext->data_size
                               : (estimate + max_error_ + 2);

        size_t pos = BinarySearch(ext->data, begin, end, key);

        if ((*(ext->data))[pos].first == key) {
            return true;
        }

        return false;
    }

    bool Lookup(const KeyType& key, ValueType& value) {
        //   lookup_count_++;
        size_t index = 1;
        if (key > min_key_.first && key <= max_key_.first) {
            index = GetSplineSegment(key);
        } else if (key > max_key_.first) {
            index = spline_points_.size();
        }

        Coord<KeyType>& down = spline_points_[index - 1];

        extra<KeyType, ValueType>* ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(down));
        if (needCheckBuffer(down)) {
            auto it = ext->buffer->find(key);
            if (it != ext->buffer->end()) {
                // check if deleted
                if (it->second != (KeyType)(-1)) {
                    value = it->second;
                    return true;
                }
                return false;
            }
        }
        size_t estimate = 0;
        if (key >= min_key_.first && key <= max_key_.first) {
            Coord<KeyType>& up = spline_points_[index];

            // Compute slope.
            const double x_diff = up.x - down.x;
            const double y_diff = up.y - down.y;
            const double slope = y_diff / x_diff;

            // Interpolate.
            const double key_diff = key - down.x;
            estimate = std::fma(key_diff, slope, down.y);

            if ((LTYPE)getType(down) == kRsSeg) {
                estimate = estimate - (size_t)(down.rsseg_start_pos);
            } else {
                estimate = estimate - (size_t)(down.y);
            }
        } else {
            return false;
        }

        if ((*(ext->data))[estimate].first == key) {
            // hit_count_++;
            value = (*(ext->data))[estimate].second;
            return true;
        }

        if ((LTYPE)getType(down) == kLine) {
            // if the searching key belongs to a line but not found, no need to
            // do following BS
            return false;
        }

        const size_t begin =
            (estimate < max_error_) ? 0 : (estimate - max_error_);
        const size_t end = (estimate + max_error_ + 2 > ext->data_size)
                               ? ext->data_size
                               : (estimate + max_error_ + 2);

        size_t pos = BinarySearch(ext->data, begin, end, key);

        if ((*(ext->data))[pos].first == key) {
            value = (*(ext->data))[pos].second;
            return true;
        }
        return false;
    }

    // insert the key; point.second is meaningless now
    bool Insert(const std::pair<KeyType, ValueType>& point) {
        KeyType key = point.first;
        uint64_t index = 1;
        if (key > min_key_.first && key <= max_key_.first) {
            index = GetSplineSegment(key);
        } else if (key > max_key_.first) {
            index = spline_points_.size();
        }

        Coord<KeyType>& down = spline_points_[index - 1];

        extra<KeyType, ValueType>* ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(down));

        (*(ext->buffer))[key] = point.second;
        setBufferHasDataBit(down);

        return true;
    }

    // compact buffer with the line/seg below it
    void DoCompactAll() {
        extra<KeyType, ValueType>* prev_seg_data = nullptr;
        std::vector<value_type> merged_vec;

        std::vector<Coord<KeyType>> merged_spline_points;

        for (size_t i = 0; i < spline_points_.size(); ++i) {
            auto item = spline_points_[i];
            // LTYPE dtype = (LTYPE)getType (item);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(item));

            if (ext == prev_seg_data) {
                continue;
            }

            prev_seg_data = ext;
            if (needCheckBuffer(item)) {
                std::vector<value_type> buffer_vec = GetBufferVec(ext->buffer);
                bool is_include_first = (i == 0) ? true : false;
                std::vector<value_type> data_vec = GetDataVec(
                    ext->data, ext->data_size, buffer_vec, is_include_first);

                // remove deleted keys (tombstone: -1)
                for (size_t k = 0; k < buffer_vec.size(); ++k) {
                    if (buffer_vec[k].second != (ValueType)(-1)) {
                        merged_vec.push_back(buffer_vec[k]);
                    }
                }

                for (auto it : data_vec) {
                    merged_vec.push_back(it);
                }
                clearBufferHasDataBit(item);
                // clear current buffer
                ext->buffer->clear();  // for tlx::btree_map
                ext->buffer = nullptr;
            } else {
                // for last spline point, there is no data; may only have buffer
                if (i != spline_points_.size() - 1) {
                    if (i == 0) {  // first spline points: add all
                        for (size_t k = 0; k < ext->data_size; ++k) {
                            merged_vec.push_back(
                                std::make_pair((*(ext->data))[k].first,
                                               (*(ext->data))[k].second));
                        }
                    } else {  // exclude the start point
                        for (size_t k = 1; k < ext->data_size; ++k) {
                            merged_vec.push_back(
                                std::make_pair((*(ext->data))[k].first,
                                               (*(ext->data))[k].second));
                        }
                    }
                }
            }
            // clear data
            if (ext->data != nullptr) {
                delete ext->data;
                ext->data = nullptr;
            }
            ext->data_size = 0;
            if (ext != nullptr) {
                delete ext;
                ext = nullptr;
            }
        }

        FormatKeys<value_type>(merged_vec);

        // update basic info
        std::cout << "[DoCompactAll] Before compaction, min_key=("
                  << min_key_.first << "," << min_key_.second << "), max_key=("
                  << max_key_.first << "," << max_key_.second
                  << "), num_keys=" << num_keys_ << std::endl;
        if (merged_vec.empty()) {
            min_key_ = std::make_pair(0, 0);
            max_key_ = std::make_pair(0, 0);
            num_keys_ = 0;
            num_shift_bits_ = 0;
            radix_table_.clear();
            radix_table_.resize(2, 0);
            spline_points_.clear();
            spline_points_.shrink_to_fit();
            std::cerr << "[DoCompactAll] Error !" << std::endl;
            // TODO: need to deal with this case separately
            // since RS needs min-key, max-key to build model
            return;
        } else {
            min_key_ = merged_vec.front();
            max_key_ = merged_vec.back();
            num_keys_ = merged_vec.size();

            num_shift_bits_ = GetNumShiftBits(max_key_.first - min_key_.first,
                                              num_radix_bits_);

            std::vector<uint32_t> merged_radix_table;
            const uint32_t max_prefix =
                (max_key_.first - min_key_.first) >> num_shift_bits_;
            merged_radix_table.resize(max_prefix + 2, 0);

            build_model(merged_vec, merged_spline_points, merged_radix_table);

            // replace spline_points_ and radix_table_ with new ones
            radix_table_.clear();
            // Initialize radix table, needs to contain all prefixes up to the
            // largest
            radix_table_.resize(max_prefix + 2, 0);
            spline_points_.clear();
            //{ std::vector<Coord<KeyType>> ().swap (spline_points_); }
            // spline_points_.assign (merged_spline_points.begin (),
            // merged_spline_points.end ());
            for (const auto item : merged_spline_points) {
                spline_points_.push_back(item);
            }
            radix_table_.assign(merged_radix_table.begin(),
                                merged_radix_table.end());
        }
        {
            std::vector<value_type>().swap(merged_vec);
        }
        std::cout
            << "[CompactAll] Compact All Finished! After compaction, min_key=("
            << min_key_.first << "," << min_key_.second << "), max_key=("
            << max_key_.first << "," << max_key_.second
            << "), num_keys=" << num_keys_ << std::endl;
    }

    // must be called after insert
    void Statistics() {
        std::cout << "---------- [Lisa-Statistics] ----------" << std::endl;
        std::cout << "min_key: " << min_key_.first << std::endl;
        std::cout << "max_key: " << max_key_.first << std::endl;
        std::cout << "num_keys(us): " << num_keys_
                  << std::endl;  // the total num of keys that are written to
                                 // the model
        std::cout << "num_radix_bits: " << num_radix_bits_ << std::endl;
        std::cout << "num_shift_bits: " << num_shift_bits_ << std::endl;
        std::cout << "max_error: " << max_error_ << std::endl;
        std::cout << "radix_table size: " << radix_table_.size() << std::endl;
        std::cout << "spline_points size: " << spline_points_.size()
                  << std::endl;
        std::cout << "# of segments: " << spline_points_.size() - 1
                  << std::endl;
        std::cout << "lookup_count_: " << lookup_count_ << std::endl;
        std::cout << "hit_count_(Lookup): " << hit_count_ << std::endl;
        std::cout << "miss_count_: " << miss_count_ << std::endl;
        std::cout << "hit_ratio: " << GetLookupHitRatio() << std::endl;

        // statistics of diff type
        // statistics of diff type
        size_t type_line_num = 0;
        size_t line_buffer_num = 0;
        size_t seg_buffer_num = 0;
        size_t rs_seg_num = 0;
        size_t bs_num = 0;
        size_t btree_num = 0;
        size_t invalid_num = 0;
        std::unordered_set<size_t>
            points_on_line_set;  // points on line(line has buffer)
        std::unordered_set<size_t>
            points_on_line_buffer_set;  // points on line's buffer(line has
                                        // buffer)
        std::unordered_set<size_t>
            points_on_line_nobuffer_set;  // points on line(line has no buffer)
        std::unordered_set<size_t> points_on_seg_nobuffer_set;
        std::unordered_set<size_t> points_on_seg_set;
        std::unordered_set<size_t> points_on_seg_buffer_set;
        std::unordered_set<size_t>
            points_on_last_sp_buffer_set;  // last spline point's buffer
        std::unordered_set<size_t> points_on_btree_set;

        // deal with the first line/seg
        if (!spline_points_.empty()) {
            LTYPE dtype = (LTYPE)getType(spline_points_[0]);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(
                    getAddr(spline_points_[0]));
            if (dtype == kLine) {
                type_line_num++;

                if (hasBuffer(spline_points_[0])) {
                    line_buffer_num++;

                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        if (it->second != (ValueType)(-1)) {
                            points_on_line_buffer_set.insert(it->first);
                        }
                    }
                    // line w/ buffer
                    for (size_t i = 0; i < ext->data_size; ++i) {
                        points_on_line_set.insert((*(ext->data))[i].first);
                    }
                } else {
                    // line w/o buffer
                    for (size_t i = 0; i < ext->data_size; ++i) {
                        points_on_line_nobuffer_set.insert(
                            (*(ext->data))[i].first);
                    }
                }
            } else if (dtype == kRsSeg) {
                rs_seg_num++;

                if (hasBuffer(spline_points_[0])) {
                    seg_buffer_num++;
                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        if (it->second != (ValueType)(-1)) {
                            points_on_seg_buffer_set.insert(it->first);
                        }
                    }
                }
                for (size_t i = 0; i < ext->data_size; ++i) {
                    points_on_seg_set.insert((*(ext->data))[i].first);
                }
            }
        }

        for (size_t i = 1; i < spline_points_.size(); ++i) {
            auto& start_point = spline_points_[i];
            LTYPE dtype = (LTYPE)getType(start_point);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(
                    getAddr(start_point));

            if (dtype == kLine) {
                type_line_num++;

                if (hasBuffer(start_point)) {
                    line_buffer_num++;

                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        if (it->second != (ValueType)(-1)) {
                            points_on_line_buffer_set.insert(it->first);
                        }
                    }
                    // exclude the first key(it's been dealt with on previous
                    // line/seg)
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_line_set.insert((*(ext->data))[j].first);
                    }
                } else {
                    // exclude the first key(it's been dealt with on previous
                    // line/seg)
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_line_nobuffer_set.insert(
                            (*(ext->data))[j].first);
                    }
                }
            } else if (dtype == kRsSeg) {
                rs_seg_num++;

                if (hasBuffer(start_point)) {
                    seg_buffer_num++;
                    if (i == spline_points_.size() - 1) {
                        for (auto it = ext->buffer->begin();
                             it != ext->buffer->end(); ++it) {
                            if (it->second != (ValueType)(-1)) {
                                points_on_last_sp_buffer_set.insert(it->first);
                            }
                        }
                    }

                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        if (it->second != (ValueType)(-1)) {
                            points_on_seg_buffer_set.insert(it->first);
                        }
                    }
                }
                // exclude the first point
                for (size_t j = 1; j < ext->data_size; ++j) {
                    points_on_seg_set.insert((*(ext->data))[j].first);
                }
            } else {
                invalid_num++;
            }
        }

        // debug
        // SaveSet2File (points_on_seg_buffer_set, "./test_seg_buffer.csv");

        std::cout << "[line type num] "
                  << "," << type_line_num << "," << line_buffer_num << ","
                  << rs_seg_num << "," << seg_buffer_num << "," << bs_num << ","
                  << btree_num << "," << invalid_num << std::endl;
        std::cout << "num of lines(buffer + nobuffer): " << type_line_num
                  << std::endl;
        std::cout << "num of lines(buffer): " << line_buffer_num << std::endl;
        std::cout << "num of lines(nobuffer): "
                  << type_line_num - line_buffer_num << std::endl;
        std::cout << "num of seg(buffer + nobuffer): " << rs_seg_num
                  << std::endl;
        std::cout << "num of seg(buffer): " << seg_buffer_num << std::endl;
        std::cout << "num of seg(nobuffer): " << rs_seg_num - seg_buffer_num
                  << std::endl;
        std::cout << "num of line/seg/btree: " << spline_points_.size()
                  << " == " << type_line_num + rs_seg_num + btree_num
                  << std::endl;
        std::cout << "---------------" << std::endl;
        std::cout << "points on pure line(nobuffer): "
                  << points_on_line_nobuffer_set.size() << std::endl;
        std::cout << "points on line(has buffer): " << points_on_line_set.size()
                  << std::endl;
        std::cout << "points on line's buffer: "
                  << points_on_line_buffer_set.size() << std::endl;
        std::cout << "points on pure seg(nobuffer): "
                  << points_on_seg_nobuffer_set.size() << std::endl;
        std::cout << "points on seg(has buffer): " << points_on_seg_set.size()
                  << std::endl;
        std::cout << "points on seg's buffer: "
                  << points_on_seg_buffer_set.size() << std::endl;
        std::cout << "points on last sp's buffer: "
                  << points_on_last_sp_buffer_set.size() << std::endl;
        std::cout << "Total unique keys num="
                  << points_on_line_set.size() +
                         points_on_line_buffer_set.size() +
                         points_on_line_nobuffer_set.size() +
                         points_on_seg_set.size() +
                         points_on_seg_buffer_set.size() +
                         points_on_seg_nobuffer_set.size() +
                         points_on_btree_set.size()
                  << std::endl;
    }

    // Returns the size in bytes.
    size_t GetIndexSize() const {
        size_t metadata = sizeof(*this);
        size_t radix_table_size = radix_table_.size() * sizeof(uint32_t);
        size_t data_size = data_.size() * sizeof(value_type);
        size_t spline_points_size = 0;
        // spline_points_.size () * sizeof (Coord<KeyType>);
        size_t delta_index_size = 0;
        for (size_t i = 0; i < spline_points_.size(); ++i) {
            auto point = spline_points_[i];
            spline_points_size +=
                sizeof(Coord<KeyType>);  // metadata of spline point
            spline_points_size +=
                sizeof(extra<KeyType, ValueType>);  // metadata of extra

            LTYPE type = (LTYPE)getType(point);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(point));
            if (hasBuffer(point)) {
                // calculate the index size of the buffer btree
                spline_points_size += ext->buffer->index_size();
                delta_index_size += ext->buffer->index_size();
            }
        }

        size_t total_index_size =
            metadata + radix_table_size + spline_points_size + data_size;
        std::cout << "[GetIndexSize] debug: metadata=" << metadata << std::endl;
        std::cout << "radix_table_size=" << radix_table_size
                  << "; size=" << radix_table_.size() << " * 4" << std::endl;
        std::cout << "spline_points_.size()=" << spline_points_.size()
                  << std::endl;
        std::cout << "spline_points_size=" << spline_points_size
                  << ", size of coord = " << sizeof(Coord<KeyType>)
                  << " * size=" << spline_points_.size() << std::endl;
        std::cout << "data_size=" << data_size
                  << ", size of value_type = " << sizeof(value_type)
                  << "; data_.size=" << data_.size() << std::endl;
        std::cout << "total index size = " << total_index_size << std::endl;
        std::cout << "size of extra=" << sizeof(extra<KeyType, ValueType>)
                  << std::endl;
        std::cout << "Delta index size=" << delta_index_size << std::endl;
        return total_index_size;
    }

    size_t GetTotalSize() const {
        size_t metadata = sizeof(*this);
        size_t radix_table_size = radix_table_.size() * sizeof(uint32_t);
        size_t data_size = data_.size() * sizeof(value_type);
        size_t spline_points_size = 0;
        // spline_points_.size () * sizeof (Coord<KeyType>);
        size_t delta_total_size = 0;
        size_t ext_total_datasize = 0;
        size_t ext_total_data = 0;
        for (size_t i = 0; i < spline_points_.size(); ++i) {
            auto point = spline_points_[i];
            spline_points_size +=
                sizeof(Coord<KeyType>);  // metadata of spline point
            spline_points_size +=
                sizeof(extra<KeyType, ValueType>);  // metadata of extra

            LTYPE type = (LTYPE)getType(point);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(point));
            if (hasBuffer(point)) {
                // calculate the index size of the buffer btree
                spline_points_size += ext->buffer->total_size();
                spline_points_size += (ext->data_size * sizeof(value_type));

                delta_total_size += ext->buffer->total_size();
                ext_total_datasize += ext->data_size;
                ext_total_data += (ext->data_size * sizeof(value_type));
            }
        }

        size_t total_index_size =
            metadata + radix_table_size + spline_points_size + data_size;
        std::cout << "[GetTotalSize] debug: metadata=" << metadata << std::endl;
        std::cout << "radix_table_size=" << radix_table_size
                  << "; size=" << radix_table_.size() << " * 4" << std::endl;
        std::cout << "spline_points_.size()=" << spline_points_.size()
                  << std::endl;
        std::cout << "spline_points_size=" << spline_points_size
                  << ", size of coord = " << sizeof(Coord<KeyType>)
                  << " * size=" << spline_points_.size() << std::endl;
        std::cout << "data_size=" << data_size
                  << ", size of value_type = " << sizeof(value_type)
                  << "; data_.size=" << data_.size() << std::endl;
        std::cout << "total index size = " << total_index_size << std::endl;
        std::cout << "size of extra=" << sizeof(extra<KeyType, ValueType>)
                  << std::endl;
        std::cout << "Delta total size=" << delta_total_size << std::endl;
        std::cout << "Extra total datasize=" << ext_total_datasize << std::endl;
        std::cout << "Extra data size=" << ext_total_data << std::endl;
        return total_index_size;
    }
    size_t GetHitCount() { return hit_num_; }
    double GetHitRatio() { return ((double)hit_num_ / num_keys_) * 100; }
    size_t GetAvgLineLength() { return avg_line_length_; }
    size_t GetLookupHitCount() { return hit_count_; }
    double GetLookupHitRatio() {
        return !lookup_count_ ? 0 : ((double)hit_count_ / lookup_count_) * 100;
    }
    size_t GetSegmentNum() { return spline_points_.size() - 1; }

    // Save the lines in the model into a file (*.csv)
    // Dubug function (Only work for small demo trace)
    void SaveLinesInfoToFile(const std::string& output,
                             bool is_print_points = false) {
        std::ofstream log(output);
        if (!log.is_open()) {
            std::cout << "[LISA SaveLinesInfoToFile] Fail to open the file: "
                      << output << std::endl;
            return;
        }

        log << "sp_index,sp_x,sp_y,sp_rsseg_start_pos,sp_type_pointer,sp_line_"
               "type,sp_data_size\n";
        for (size_t i = 0; i < spline_points_.size(); ++i) {
            auto sp = spline_points_[i];
            LTYPE type = (LTYPE)getType(sp);
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(sp));
            std::string line_type = (type == kLine) ? "Line" : "RSSeg";
            log << i << "," << sp.x << "," << sp.y << "," << sp.rsseg_start_pos
                << "," << sp.type_pointer << "," << line_type << ","
                << ext->data_size << "\n";
            if (is_print_points && ext->data_size > 0) {
                log << "Points on this line:\n";
                for (size_t j = 0; j < ext->data_size; ++j) {
                    log << "(" << (*(ext->data))[j].first << ","
                        << (*(ext->data))[j].second << "),";
                    if ((j + 1) % 10 == 0) {
                        log << "\n";
                    }
                }
                log << "\n";
            }
        }
        std::cout << "[SaveLinesInfoToFile] Successfully save lines info into "
                  << output << std::endl;
        log.close();
    }

    /* void Compact () {
        auto merged_spline_points = DoCompactAll ();
        update_num_keys (merged_spline_points);
        std::cout << "[Compact] After compaction, min_key=(" << min_key_.first
    << ","
                  << min_key_.second << "), max_key=(" << max_key_.first << ","
    << max_key_.second
                  << "), num_keys=" << num_keys_ << std::endl;
        std::vector<uint32_t> merged_radix_table;
        const uint32_t max_prefix = (max_key_.first - min_key_.first) >>
    num_shift_bits_; merged_radix_table.resize (max_prefix + 2, 0);
        update_radix_table (merged_spline_points, merged_radix_table);
        spline_points_.clear ();
        spline_points_ = merged_spline_points;
    } */

    void update_num_keys(std::vector<Coord<KeyType>>& merged_spline_points) {
        std::set<value_type> all_keys_set;

        for (size_t i = 0; i < merged_spline_points.size(); ++i) {
            auto start_point = merged_spline_points[i];
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(
                    getAddr(start_point));
            if (hasBuffer(start_point)) {
                /* std::vector<value_type> buffer_data;
                start_point.buffer->scan_to_end (0, buffer_data);
                for (auto item : buffer_data) {
                    all_keys_set.insert (item);
                } */
                for (auto it = ext->buffer->begin(); it != ext->buffer->end();
                     ++it) {
                    all_keys_set.insert(std::make_pair(it->first, it->second));
                }
            } else {
                for (size_t j = 0; j < ext->data_size; ++j) {
                    all_keys_set.insert(std::make_pair(
                        (*(ext->data))[j].first, (*(ext->data))[j].second));
                }
            }
        }
        num_keys_ = all_keys_set.size();
        min_key_ = *(all_keys_set.begin());
        max_key_ = *(std::prev(all_keys_set.end()));
    }

   private:
    void SaveVec2File(const std::vector<KeyType>& keys,
                      const std::string output) {
        std::ofstream log(output);
        if (!log.is_open()) {
            std::cout << "[SaveVec2File] Fail to open the file: " << output
                      << std::endl;
            return;
        }

        for (size_t i = 0; i < keys.size(); ++i) {
            log << i << "," << keys[i] << "\n";
        }
        std::cout << "[SaveVec2File] Successfully save lines info into "
                  << output << std::endl;
        log.close();
    }

    void SaveSet2File(const std::unordered_set<size_t>& keys,
                      const std::string output) {
        std::ofstream log(output);
        if (!log.is_open()) {
            std::cout << "[SaveSet2File] Fail to open the file: " << output
                      << std::endl;
            return;
        }
        size_t i = 0;
        for (auto it = keys.begin(); it != keys.end(); ++it) {
            log << *it << "\n";
        }
        std::cout << "[SaveSet2File] Successfully save lines info into "
                  << output << std::endl;
        log.close();
    }

    // Returns the number of shift bits based on the `diff` between the largest
    // and the smallest key. KeyType == uint32_t.
    static size_t GetNumShiftBits(uint32_t diff, size_t num_radix_bits) {
        const uint32_t clz = __builtin_clz(diff);
        if ((32 - clz) < num_radix_bits) return 0;
        return 32 - num_radix_bits - clz;
    }
    // KeyType == uint64_t.
    static size_t GetNumShiftBits(uint64_t diff, size_t num_radix_bits) {
        const uint32_t clzl = __builtin_clzl(diff);
        if ((64 - clzl) < num_radix_bits) return 0;
        return 64 - num_radix_bits - clzl;
    }
    // KeyType == int64_t.
    static size_t GetNumShiftBits(int64_t diff, size_t num_radix_bits) {
        const uint32_t clzl =
            __builtin_clzl(static_cast<uint64_t>(std::floor(diff)));
        if ((64 - clzl) < num_radix_bits) return 0;
        return 64 - num_radix_bits - clzl;
    }

    // Returns the index of the spline point that marks the end of the spline
    // segment that contains the `key`: `key` ∈ (spline[index - 1],
    // spline[index]]
    size_t GetSplineSegment(const KeyType key) const {
        // Narrow search range using radix table.
        const KeyType prefix = (key - min_key_.first) >> num_shift_bits_;
        assert(prefix + 1 < radix_table_.size());
        const uint32_t begin = radix_table_[prefix];
        const uint32_t end = radix_table_[prefix + 1];

        if (end - begin < 32) {
            // Do linear search over narrowed range.
            uint32_t current = begin;
            while (spline_points_[current].x < key) {
                ++current;
            }
            return current;
        }
        // Do binary search over narrowed range.
        const auto lb = std::lower_bound(
            spline_points_.begin() + begin, spline_points_.begin() + end, key,
            [](const Coord<KeyType>& coord, const KeyType key) {
                return coord.x < key;
            });

        return std::distance(spline_points_.begin(), lb);
    }

    static size_t BinarySearch(const std::vector<value_type>* array,
                               const size_t& begin, const size_t& end,
                               const KeyType& target) {
        int lower = begin;
        int upper = end;
        while (lower <= upper) {
            size_t mid = lower + (upper - lower) / 2;
            if ((*array)[mid].first == target) {
                return mid;
            } else if ((*array)[mid].first < target) {
                lower = mid + 1;
            } else {
                upper = mid - 1;
            }
        }
        return lower;
    }

    static size_t BinarySearch(const std::vector<value_type> array,
                               const size_t& begin, const size_t& end,
                               const KeyType& target) {
        int lower = begin;
        int upper = end;
        while (lower <= upper) {
            size_t mid = lower + (upper - lower) / 2;
            if (array[mid].first == target) {
                return mid;
            } else if (array[mid].first < target) {
                lower = mid + 1;
            } else {
                upper = mid - 1;
            }
        }
        return lower;
    }

    std::vector<value_type> GetBufferVec(
        tlx::btree_map<KeyType, ValueType>* const buffer) {
        std::vector<value_type> buffer_vec;
        for (auto it = buffer->begin(); it != buffer->end(); ++it) {
            buffer_vec.push_back(std::make_pair(it->first, it->second));
        }
        // return std::move (buffer_vec);
        return buffer_vec;
    }

    // get data from *data and remove element that in the buffer
    std::vector<value_type> GetDataVec(std::vector<value_type>* const data,
                                       const size_t size,
                                       const std::vector<value_type> buffer,
                                       bool is_first_point) {
        std::vector<value_type> data_vec;

        // create an unordered_set to store keys from the buffer vector
        std::unordered_set<KeyType> buffer_set;
        for (const auto& item : buffer) {
            buffer_set.insert(item.first);
        }
        // if it's not first spine point, exclude the start point(since it's
        // already dealt with as the end point in previous spline point)
        size_t start_index = (is_first_point) ? 0 : 1;
        for (size_t i = start_index; i < size; ++i) {
            auto key = (*data)[i].first;

            // check if the key exists in the buffer_set
            if (buffer_set.find(key) == buffer_set.end()) {
                data_vec.push_back((*data)[i]);
            }
        }
        return data_vec;
    }

    // build model with unique and sorted keys
    // input: unique and sorted keys
    // output: model(line+rsseg/btree) based on input keys (spline_points and
    // radix_table)
    bool build_model(std::vector<value_type>& input,
                     std::vector<Coord<KeyType>>& spline_points,
                     std::vector<uint32_t>& radix_table) {
        if (input.empty()) {
            return true;
        }
        // step1: find lines and sand(only contain two endpoints) within keys
        std::vector<Line<KeyType, ValueType>> line_vec =
            find_model_line(input, 0);

        // step2: merge consecutive sand together
        std::vector<Line<KeyType, ValueType>> line_sand_vec =
            merge_sand(line_vec, alpha_);
        /* std::vector<Line<KeyType, ValueType>>* line_sand_ptr =
            new std::vector<Line<KeyType, ValueType>>(); */

        // step3: categorize clumping sand into kRsSeg or kBtree
        // and get spline points
        spline_points = get_spline_points(line_sand_vec, 0, true);

        // debug
        // show_spline_points (spline_points);

        // step4: update radix_table
        // traverse spline_points and update radix_table
        update_radix_table(spline_points, radix_table);

        return true;
    }

    // make keys unique and sorted
    template <class FType>
    void FormatKeys(std::vector<FType>& keys) {
        tbb::parallel_sort(keys.begin(), keys.end());
        size_t before_size = keys.size();
        auto it = std::unique(keys.begin(), keys.end());
        keys.erase(it, keys.end());
        std::cout << "[FormatKeys] Before: " << before_size
                  << "; after: " << keys.size() << std::endl;
    }

    // make keys unique and sorted
    /* void FormatKeys (std::vector<value_type>& keys) {
        tbb::parallel_sort (keys.begin (), keys.end ());
        size_t before_size = keys.size ();
        auto it = std::unique (keys.begin (), keys.end ());
        keys.erase (it, keys.end ());
        std::cout << "[FormatKeys] Before: " << before_size << "; after: " <<
    keys.size ()
                  << std::endl;
    } */

    // input: input keys (unique and sorted); begin position of the first key
    // output: vector with line and nonline(only two points)
    std::vector<Line<KeyType, ValueType>> find_model_line(
        std::vector<value_type>& input, const double begin_pos) {
        KeyType prev_key = 0;
        double prev_pos = begin_pos;
        KeyType curr_key = 0;
        double curr_pos = begin_pos;
        size_t prev_difference = 0;

        Line<KeyType, ValueType> current_line;
        std::vector<Line<KeyType, ValueType>> line_nonline_vec;

        for (size_t i = 0; i < input.size(); ++i) {
            curr_key = input[i].first;

            if (0 == i) {
                prev_key = curr_key;
                prev_pos = curr_pos;
                // current_line.start = {curr_key, curr_pos, 0, 0};
                current_line.start.x = curr_key;
                current_line.start.y = curr_pos;
                current_line.start.rsseg_start_pos = 0;
                current_line.start.type_pointer = 0;
                // current_line.end = {curr_key, curr_pos, 0, 0};
                current_line.end.x = curr_key;
                current_line.end.y = curr_pos;
                current_line.end.rsseg_start_pos = 0;
                current_line.end.type_pointer = 0;
                current_line.is_line = false;
                current_line.length = 1;
                current_line.points.push_back(
                    std::make_pair(curr_key, input[i].second));
                curr_pos += 1;
                continue;
            } else if (1 == i) {
                // current_line.end = {curr_key, curr_pos, 0, 0};
                current_line.end.x = curr_key;
                current_line.end.y = curr_pos;
                current_line.end.rsseg_start_pos = 0;
                current_line.end.type_pointer = 0;
                current_line.length++;
                current_line.points.push_back(
                    std::make_pair(curr_key, input[i].second));
                prev_difference = curr_key - prev_key;
                prev_key = curr_key;
                prev_pos = curr_pos;
                curr_pos += 1;
                continue;
            }

            size_t curr_difference = curr_key - prev_key;
            if (curr_difference == prev_difference) {
                current_line.end.x = curr_key;
                current_line.end.y = curr_pos;
                current_line.length++;
                current_line.points.push_back(
                    std::make_pair(curr_key, input[i].second));
            } else {
                if (current_line.length >= alpha_) {
                    current_line.is_line = true;
                }
                line_nonline_vec.push_back(current_line);

                // make the line adjacent -> easier to look up the line
                current_line.start.x = prev_key;
                current_line.start.y = prev_pos;
                current_line.end.x = curr_key;
                current_line.end.y = curr_pos;
                current_line.is_line = false;
                current_line.length = 2;
                current_line.points.clear();
                current_line.points.push_back(
                    std::make_pair(prev_key, input[i - 1].second));
                current_line.points.push_back(
                    std::make_pair(curr_key, input[i].second));
                prev_difference = curr_difference;
            }

            prev_pos = curr_pos;
            prev_key = curr_key;
            curr_pos += 1;
        }

        // deal with the last line
        if (current_line.length >= alpha_) {
            current_line.is_line = true;
        }
        line_nonline_vec.push_back(current_line);
        // return std::move (line_nonline_vec);
        return line_nonline_vec;
    }

    // input: Lines with line and sand
    // merge consecutive sand into one clumping sand
    std::vector<Line<KeyType, ValueType>> merge_sand(
        std::vector<Line<KeyType, ValueType>>& line_vec, const size_t alpha) {
        std::vector<Line<KeyType, ValueType>>
            line_sand_vec;  // contain line or clumping sand
        Line<KeyType, ValueType>
            clumping_sand;  // clumping sand: not actual line
        bool is_start_set = false;
        for (size_t i = 0; i < line_vec.size(); ++i) {
            if (line_vec[i].length < alpha)  //
            {
                if (!is_start_set) {
                    // deal the first clumping sand
                    is_start_set = true;
                    // clumping_sand.start = line_vec[i].start;
                    clumping_sand.start.x = line_vec[i].start.x;
                    clumping_sand.start.y = line_vec[i].start.y;
                    clumping_sand.start.rsseg_start_pos =
                        line_vec[i].start.rsseg_start_pos;
                    clumping_sand.start.type_pointer =
                        line_vec[i].start.type_pointer;
                    // clumping_sand.end = line_vec[i].end;
                    clumping_sand.end.x = line_vec[i].end.x;
                    clumping_sand.end.y = line_vec[i].end.y;
                    clumping_sand.end.rsseg_start_pos =
                        line_vec[i].end.rsseg_start_pos;
                    clumping_sand.end.type_pointer =
                        line_vec[i].end.type_pointer;
                    clumping_sand.points.insert(clumping_sand.points.end(),
                                                line_vec[i].points.begin(),
                                                line_vec[i].points.end());
                    clumping_sand.length += line_vec[i].length;
                } else if (clumping_sand.end.x != line_vec[i].start.x) {
                    // clumping_sand.end = line_vec[i].end;
                    clumping_sand.end.x = line_vec[i].end.x;
                    clumping_sand.end.y = line_vec[i].end.y;
                    clumping_sand.end.rsseg_start_pos =
                        line_vec[i].end.rsseg_start_pos;
                    clumping_sand.end.type_pointer =
                        line_vec[i].end.type_pointer;
                    clumping_sand.points.insert(clumping_sand.points.end(),
                                                line_vec[i].points.begin(),
                                                line_vec[i].points.end());
                    clumping_sand.length += line_vec[i].length;
                } else {
                    // clumping_sand.end = line_vec[i].end;
                    clumping_sand.end.x = line_vec[i].end.x;
                    clumping_sand.end.y = line_vec[i].end.y;
                    clumping_sand.end.rsseg_start_pos =
                        line_vec[i].end.rsseg_start_pos;
                    clumping_sand.end.type_pointer =
                        line_vec[i].end.type_pointer;
                    clumping_sand.points.insert(clumping_sand.points.end(),
                                                line_vec[i].points.begin() + 1,
                                                line_vec[i].points.end());
                    clumping_sand.length =
                        clumping_sand.length + line_vec[i].length - 1;
                }

                // deal with the last clumping sand
                if (i == line_vec.size() - 1) {
                    clumping_sand.is_line = false;
                    line_sand_vec.push_back(clumping_sand);
                }
            } else  // deal with lines (with 100% hit ratio)
            {
                // deal with the previous vibrated segment
                if (clumping_sand.length != 0) {
                    clumping_sand.is_line = false;
                    line_sand_vec.push_back(clumping_sand);
                }
                line_vec[i].is_line = true;
                line_sand_vec.push_back(line_vec[i]);
                clumping_sand.Clear();
                is_start_set = false;
            }
        }
        // return std::move (line_sand_vec);
        return line_sand_vec;
    }

    // get spline points
    // offset: indicates the offset of first spline point in spline_points;
    // it's the value of previous spline
    // is_last_lineseg: add new btree only when this value is true
    std::vector<Coord<KeyType>> get_spline_points(
        std::vector<Line<KeyType, ValueType>>& line_sand_vec,
        const double offset, const bool is_last_lineseg = true) {
        std::vector<Coord<KeyType>> spline_points;

        for (size_t i = 0; i < line_sand_vec.size(); ++i) {
            if (line_sand_vec[i].is_line)  // lines
            {
                extra<KeyType, ValueType>* ext =
                    new extra<KeyType, ValueType>();
                std::vector<value_type>* data = new std::vector<value_type>();
                tlx::btree_map<KeyType, ValueType>* buffer =
                    new tlx::btree_map<KeyType, ValueType>();
                // copy data
                for (size_t j = 0; j < line_sand_vec[i].points.size(); ++j) {
                    data->push_back(line_sand_vec[i].points[j]);
                }
                ext->data = data;
                ext->data_size = line_sand_vec[i].points.size();
                ext->buffer = buffer;

                setType(kLine, line_sand_vec[i].start);
                setAddr((void*)ext, line_sand_vec[i].start);
                setBufferBit(line_sand_vec[i].start);
                spline_points.push_back(line_sand_vec[i].start);

            } else  // nonline
            {
                extra<KeyType, ValueType>* ext =
                    new extra<KeyType, ValueType>();
                std::vector<value_type>* data = new std::vector<value_type>();
                tlx::btree_map<KeyType, ValueType>* buffer =
                    new tlx::btree_map<KeyType, ValueType>();
                // copy data
                for (size_t j = 0; j < line_sand_vec[i].points.size(); ++j) {
                    data->push_back(line_sand_vec[i].points[j]);
                }
                ext->data = data;
                ext->data_size = line_sand_vec[i].points.size();
                ext->buffer = buffer;
                if (line_sand_vec[0].start.y != offset) {
                    std::cout << "[get_spline_points] warning! start.y="
                              << line_sand_vec[0].start.y
                              << "; offset=" << offset << std::endl;
                }
                double start_pos = offset;
                if (!spline_points.empty()) {
                    auto last_point_on_prev_line = spline_points.back();
                    LTYPE dtype = (LTYPE)getType(last_point_on_prev_line);
                    if (dtype != kLine) {
                        std::cerr
                            << "[get_spline_points] Error! dtype=" << dtype
                            << std::endl;
                    }
                    extra<KeyType, ValueType>* ext =
                        reinterpret_cast<extra<KeyType, ValueType>*>(
                            getAddr(last_point_on_prev_line));
                    start_pos = last_point_on_prev_line.y + ext->data_size - 1;
                }
                bool is_handle_last = false;
                std::vector<Coord<KeyType>> rs_spline_points =
                    rs_gsc(line_sand_vec[i].points, start_pos, is_handle_last);

                // for (auto item : rs_spline_points)
                double rsseg_start_pos = 0;
                for (size_t j = 0; j < rs_spline_points.size(); ++j) {
                    auto& item = rs_spline_points[j];
                    if (j == 0) {
                        rsseg_start_pos = item.y;
                    }
                    item.rsseg_start_pos = rsseg_start_pos;
                    setType(kRsSeg, item);
                    setAddr((void*)ext, item);
                    setBufferBit(item);
                    spline_points.push_back(item);
                }
            }  // end of nonline
        }  // end of for

        // deal with the last endpoint
        if (spline_points.back().x != line_sand_vec.back().end.x) {
            if (is_last_lineseg) {
                extra<KeyType, ValueType>* ext =
                    new extra<KeyType, ValueType>();
                tlx::btree_map<KeyType, ValueType>* buffer =
                    new tlx::btree_map<KeyType, ValueType>();
                ext->buffer = buffer;
                setType(kRsSeg, line_sand_vec.back().end);
                setAddr((void*)ext, line_sand_vec.back().end);
                setBufferBit(line_sand_vec.back().end);
                spline_points.push_back(line_sand_vec.back().end);
            } else {
                spline_points.push_back(line_sand_vec.back().end);
            }
        }

        // return std::move (spline_points);
        return spline_points;
    }

    // new rs's tools
    // use greedysplinecorridor algorithm to generate spline points
    // input: Keys; output: vector of spline points
    std::vector<Coord<KeyType>> rs_gsc(std::vector<value_type>& keys,
                                       const double start_pos,
                                       const bool is_handle_last = true) {
        // make sure keys are monotonically increasing
        RSBuilder<KeyType> rsbuilder(keys.front().first, keys.back().first,
                                     num_radix_bits_, max_error_);
        double pos = start_pos;
        for (const auto& item : keys) {
            KeyType key = item.first;
            rsbuilder.AddKey(key, pos);
            pos = pos + 1;
        }
        return rsbuilder.Finalize(is_handle_last);
    }

    // generate radix_table from spline_points
    void update_radix_table(std::vector<Coord<KeyType>>& spline_points,
                            std::vector<uint32_t>& radix_table) {
        KeyType prev_prefix = 0;
        for (size_t i = 0; i < spline_points.size(); ++i) {
            KeyType key = spline_points[i].x;
            const KeyType curr_prefix =
                (key - min_key_.first) >> num_shift_bits_;
            if (curr_prefix != prev_prefix) {
                const uint32_t curr_index = i;
                for (KeyType prefix = prev_prefix + 1; prefix <= curr_prefix;
                     ++prefix) {
                    radix_table[prefix] = curr_index;
                }
                prev_prefix = curr_prefix;
            }
        }
        ++prev_prefix;
        const uint32_t num_spline_points = spline_points.size();
        for (; prev_prefix < radix_table.size(); ++prev_prefix) {
            radix_table[prev_prefix] = num_spline_points;
        }
    }

    double calculate_slope(const Coord<KeyType>& p1, const Coord<KeyType>& p2) {
        return (p2.y - p1.y) / (p2.x - p1.x);
    }

    // check if line p1-p2 and line p3-p4 are on same line
    bool are_on_same_line(const Coord<KeyType>& p1, const Coord<KeyType>& p2,
                          const Coord<KeyType>& p3, const Coord<KeyType>& p4) {
        if (p1.x == p2.x || p3.x == p4.x) {
            return false;
        }
        double slope_p12 = calculate_slope(p1, p2);
        double slope_p34 = calculate_slope(p3, p4);
        if (std::abs(slope_p12 - slope_p34) < 1e-9) {
            if (p2.x == p3.x) {
                return true;
            } else {
                double intercept_p12 = p1.y - slope_p12 * p1.x;
                double intercept_p34 = p3.y - slope_p34 * p3.x;
                if (std::abs(intercept_p12 - intercept_p34) < 1e-9) {
                    return true;
                }
            }
        }

        return false;
    }

    // merge spline_points vector
    std::vector<Coord<KeyType>> merge_sp_vec(
        std::vector<std::vector<Coord<KeyType>>>& sp_vec) {
        std::vector<Coord<KeyType>> merged_spline_points;
        for (size_t i = 0; i < sp_vec.size() - 1; ++i) {
            auto& front_sp = sp_vec[i];
            // get the last line/seg of front_sp
            if (front_sp.size() < 2) {
                std::cerr << "[merge_sp_vec] Error! front_sp.size="
                          << front_sp.size() << std::endl;
                continue;
            }
            auto front_sp_lastline_start = front_sp[front_sp.size() - 2];
            auto front_sp_lastline_end = front_sp[front_sp.size() - 1];
            auto front_sp_lastline_type = getType(front_sp_lastline_start);
            if (front_sp_lastline_type != kLine) {
                // TODO: remove the last invalid point;
                front_sp.pop_back();
                continue;
            }

            // TODO: may miss some short line. should solve it

            for (size_t j = i + 1; j < sp_vec.size(); ++j) {
                auto& back_sp = sp_vec[j];

                // get the first line/seg of back_sp
                if (back_sp.size() < 2) {
                    std::cerr << "[merge_sp_vec] Error! back_sp.size="
                              << back_sp.size() << std::endl;
                    continue;
                }
                auto back_sp_firstline_start = back_sp[0];
                auto back_sp_firstline_end = back_sp[1];
                auto back_sp_firstline_type = getType(back_sp_firstline_start);
                if (back_sp_firstline_type != kLine) {
                    // remove the last invalid point in front_sp
                    front_sp.pop_back();
                    break;
                }

                // test
                if (front_sp_lastline_end.x != back_sp_firstline_start.x) {
                    std::cout
                        << "[merge_sp_vec] warning! they should be equal. p2.x="
                        << front_sp_lastline_end.x
                        << ", p3.x=" << back_sp_firstline_start.x
                        << "; p2.y=" << front_sp_lastline_end.y
                        << ", p3.y=" << back_sp_firstline_start.y << std::endl;
                    return merged_spline_points;
                }
                bool are_same_line = are_on_same_line(
                    front_sp_lastline_start, front_sp_lastline_end,
                    back_sp_firstline_start, back_sp_firstline_end);
                if (are_same_line) {
                    // show_spline_points (front_sp);
                    front_sp[front_sp.size() - 1] = back_sp_firstline_end;
                    // show_spline_points (front_sp);
                    //  update front_sp_lastline_start_ext's data
                    update_splicing_data(front_sp_lastline_start,
                                         back_sp_firstline_start);
                    // show_spline_points (front_sp);
                    // show_spline_points (back_sp);
                    //  remove the first spline point in the back_sp
                    back_sp.erase(back_sp.begin());
                    if (back_sp.size() < 2) {
                        back_sp.clear();
                    }
                    break;
                }
            }
            // show_spline_points (front_sp);
            //  TODO: add the last element of front_sp
            merged_spline_points.insert(merged_spline_points.end(),
                                        front_sp.begin(), front_sp.end());
        }
        std::cout << "Test: merged_spline_points.size="
                  << merged_spline_points.size() << std::endl;
        // remove duplicated element
        merged_spline_points.erase(std::unique(
            merged_spline_points.begin(), merged_spline_points.end(),
            [](const auto& a, const auto& b) { return a.x == b.x; }));
        // return std::move (merged_spline_points);
        return merged_spline_points;
    }

    // update splicing data
    void update_splicing_data(Coord<KeyType>& front_start,
                              Coord<KeyType>& back_start) {
        std::set<value_type> data_set;
        extra<KeyType, ValueType>* front_start_ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(front_start));
        for (size_t i = 0; i < front_start_ext->data_size; ++i) {
            data_set.insert((*(front_start_ext->data))[i]);
        }

        extra<KeyType, ValueType>* back_start_ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(back_start));
        for (size_t i = 0; i < back_start_ext->data_size; ++i) {
            data_set.insert((*(back_start_ext->data))[i]);
        }
        // clear old data and allocate new space to store new data
        front_start_ext->data->clear();
        for (const auto& element : data_set) {
            front_start_ext->data->push_back(element);
        }
    }

    // debug function
    void show_spline_points(std::vector<Coord<KeyType>> sp) {
        std::set<size_t> points_on_line_set;  // points on line(line has buffer)
        std::set<size_t> points_on_line_buffer_set;  // points on line's
                                                     // buffer(line has buffer)
        std::set<size_t>
            points_on_line_nobuffer_set;  // points on line(line has no buffer)
        std::set<size_t> points_on_seg_nobuffer_set;
        std::set<size_t> points_on_seg_set;
        std::set<size_t> points_on_seg_buffer_set;
        std::set<size_t> points_on_last_sp_buffer_set;

        for (size_t i = 0; i < sp.size(); ++i) {
            auto start_point = sp[i];
            extra<KeyType, ValueType>* ext =
                reinterpret_cast<extra<KeyType, ValueType>*>(
                    getAddr(start_point));
            LTYPE dtype = (LTYPE)getType(start_point);
            if (dtype == kLine) {
                if (hasBuffer(start_point)) {
                    /* std::vector<value_type> buffer_data;
                    start_point.buffer->scan_to_end (0, buffer_data);
                    for (auto item : buffer_data) {
                        points_on_line_buffer_set.insert (item.first);
                    } */
                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        points_on_line_buffer_set.insert(it->first);
                    }
                    // exclude the first key(it's been dealt with on previous
                    // line/seg)
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_line_set.insert((*(ext->data))[j].first);
                    }
                } else {
                    // exclude the first key(it's been dealt with on previous
                    // line/seg)
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_line_nobuffer_set.insert(
                            (*(ext->data))[j].first);
                    }
                }
            } else if (dtype == kRsSeg) {
                if (hasBuffer(start_point)) {
                    if (i == sp.size() - 1) {
                        /* std::vector<value_type> buffer_data;
                        start_point.buffer->scan_to_end (0, buffer_data);
                        for (auto item : buffer_data) {
                            points_on_last_sp_buffer_set.insert (item.first);
                        } */
                        for (auto it = ext->buffer->begin();
                             it != ext->buffer->end(); ++it) {
                            points_on_last_sp_buffer_set.insert(it->first);
                        }
                    }
                    /* std::vector<value_type> buffer_data;
                    start_point.buffer->scan_to_end (0, buffer_data);
                    for (auto item : buffer_data) {
                        points_on_seg_buffer_set.insert (item.first);
                    } */
                    for (auto it = ext->buffer->begin();
                         it != ext->buffer->end(); ++it) {
                        points_on_seg_buffer_set.insert(it->first);
                    }
                    // exclude the first point
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_seg_set.insert((*(ext->data))[j].first);
                    }
                } else {
                    // exclude the first point
                    for (size_t j = 1; j < ext->data_size; ++j) {
                        points_on_seg_set.insert((*(ext->data))[j].first);
                    }
                }
            }
        }

        std::cout << "points on pure line(nobuffer): "
                  << points_on_line_nobuffer_set.size() << std::endl;
        std::cout << "points on line(has buffer): " << points_on_line_set.size()
                  << std::endl;
        std::cout << "points on line's buffer: "
                  << points_on_line_buffer_set.size() << std::endl;
        std::cout << "points on pure seg(nobuffer): "
                  << points_on_seg_nobuffer_set.size() << std::endl;
        std::cout << "points on seg(has buffer): " << points_on_seg_set.size()
                  << std::endl;
        std::cout << "points on seg's buffer: "
                  << points_on_seg_buffer_set.size() << std::endl;
        std::cout << "points on last sp's buffer: "
                  << points_on_last_sp_buffer_set.size() << std::endl;
    }

    void show_single_spline_point(Coord<KeyType> point) {
        std::cout << "Show single spline point: " << point.x << "," << point.y
                  << "," << point.rsseg_start_pos << std::endl;
        extra<KeyType, ValueType>* ext =
            reinterpret_cast<extra<KeyType, ValueType>*>(getAddr(point));
        LTYPE dtype = (LTYPE)getType(point);
        std::string type_str = "";
        if (dtype == kLine) {
            type_str = "Line";
        } else {
            type_str = "RS-Seg";
        }
        if (hasBuffer(point)) {
            std::cout << "Has buffer: buffer size=" << ext->buffer->size()
                      << std::endl;
            for (auto it = ext->buffer->begin(); it != ext->buffer->end();
                 ++it) {
                std::cout << it->first << ",";
            }
            std::cout << std::endl;
        }
        std::string hasbuffer = "";
        if (hasBuffer(point)) {
            hasbuffer = "It has buffer.";
        } else {
            hasbuffer = "No buffer.";
        }
        std::string hasbufferData = "";
        if (hasBufferData(point)) {
            hasbufferData = "The buffer has data.";
        } else {
            hasbufferData = "The buffer is empty (no data in the buffer).";
        }
        std::string need_check_buffer = "";
        if (needCheckBuffer(point)) {
            need_check_buffer =
                "You need to check the buffer: has buffer + buffer has data.";
        } else {
            if (!hasBuffer(point)) {
                need_check_buffer = "No need to check the buffer: no buffer.";
            } else if (hasBufferData(point)) {
                need_check_buffer =
                    "Error! Don't check buffer but it has buffer and buffer "
                    "has data!";
            } else {
                need_check_buffer =
                    "No need to check the buffer: buffer is empty.";
            }
        }

        std::cout << "line type: " << type_str << std::endl;
        std::cout << "need check buffer: " << need_check_buffer << std::endl;
        std::cout << "buffer: " << hasbuffer << std::endl;
        std::cout << "buffer data: " << hasbufferData << std::endl;
        /* for (size_t i = 0; i < ext->data_size; ++i) {
            std::cout << "Data: " << (*(ext->data))[i].first << "," <<
        (*(ext->data))[i].second
                      << std::endl;
        } */
        std::cout << "End of show_single_spline_point" << std::endl;
    }

   private:
    value_type min_key_;
    value_type max_key_;
    size_t num_keys_;
    size_t num_radix_bits_;
    size_t num_shift_bits_;
    size_t max_error_;  // epsilon

    std::vector<uint32_t> radix_table_;
    std::vector<Coord<KeyType>> spline_points_;

    size_t alpha_;  // minimum sequence length

    size_t hit_num_;          // hit number on real line
    size_t avg_line_length_;  // average real line length
    size_t hit_count_;        // hit count when lookup
    size_t miss_count_;
    size_t lookup_count_;  // number of lookup operations

    std::vector<value_type> data_;  // key-value
};
}  // namespace lisa

#endif  // __LISA_H__