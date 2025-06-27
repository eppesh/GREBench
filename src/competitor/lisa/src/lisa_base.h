#ifndef __LISA_BASE_H__
#define __LISA_BASE_H__

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <unordered_set>
#include <vector>

// #include "concurrency.h"
#include "btree_map.hpp"

namespace lisa {

enum LTYPE { kInvalid = 0, kLine, kRsSeg, kDirBinSearch, kBtree };

/**
 *      type       hasBuffer(60)     BufferHasData(59)           Addr
 * -----------------------------------------------------------------------
 * |   3 bits    |     1 bit    |     1 bit            |                 |
 * -----------------------------------------------------------------------
 * 60-59th bit
 * 0 0 no buffer, buffer has no data (no need)
 * 1 0 has buffer, buffer has no data (no need) // after compaction
 * 1 1 has buffer, buffer has data (check buffer)
 */

template <class KeyType>
struct Coord {
    // using value_type = std::pair<KeyType, uint64_t>;
    KeyType x;
    double y;
    // first rs-seg's start position; to calculate offset during lookup
    double rsseg_start_pos;
    uint64_t type_pointer;  // highest 3 bits represent type
};

// static_assert (sizeof (Coord<uint64_t>) == 64);
//  auto i = sizeof (Coord<uint64_t, uint64_t>);

template <class KeyType, class ValueType>
struct extra {
    using value_type = std::pair<KeyType, ValueType>;
    // buffer: point to delta (a buffer for receiving new data)
    // data: point to data array which contains all points on this line/seg
    std::vector<value_type>* data;
    size_t data_size;
    tlx::btree_map<KeyType, ValueType>* buffer;
    extra() : data(nullptr), data_size(0), buffer(nullptr) {}
};

struct SearchBound {
    size_t begin;
    size_t end;  // Exclusive.
};

template <class KeyType, class ValueType>
struct Line {
    using value_type = std::pair<KeyType, ValueType>;
    Coord<KeyType> start;
    Coord<KeyType> end;
    bool is_line;  // Line is a line or not
    size_t length;
    std::vector<value_type> points;  // TODO: change KeyType to value_type

    Line() : is_line(false), length(0) {
        start.x = 0;
        start.y = 0;
        start.rsseg_start_pos = 0;
        start.type_pointer = 0;
        end.x = 0;
        end.y = 0;
        end.rsseg_start_pos = 0;
        end.type_pointer = 0;
        points.clear();
    }

    Line& operator=(const Line& other) {
        if (this != &other) {
            start = other.start;
            end = other.end;
            length = other.length;
            is_line = other.is_line;
            points = other.points;
        }
        return *this;
    }

    void Clear() {
        // start = {0, 0, 0, 0};
        start.x = 0;
        start.y = 0;
        start.rsseg_start_pos = 0;
        start.type_pointer = 0;
        // end = {0, 0, 0, 0};
        end.x = 0;
        end.y = 0;
        end.rsseg_start_pos = 0;
        end.type_pointer = 0;
        is_line = false;
        length = 0;
        points.clear();
    }
};

static uint64_t typeMask = 0x1FFFFFFFFFFFFFFF;
static uint64_t addrMask = 0xFFFFFFFFFFFF;
static uint64_t clearBufferHasDataBitMask = 0xF000000000000000;

template <class KeyType>
static void setType(LTYPE ltype, struct Coord<KeyType>& coord) {
    uint64_t type = ltype;
    coord.type_pointer = coord.type_pointer | (type << 61);
}
template <class KeyType>
static uint8_t getType(const struct Coord<KeyType>& coord) {
    return ((coord.type_pointer & (~typeMask)) >> 61);
}
template <class KeyType>
static void setAddr(void* addr, struct Coord<KeyType>& coord) {
    coord.type_pointer =
        (coord.type_pointer & (~typeMask)) | reinterpret_cast<uint64_t>(addr);
}
template <class KeyType>
static char* getAddr(const struct Coord<KeyType>& coord) {
    return reinterpret_cast<char*>((coord.type_pointer & addrMask));
}
template <class KeyType>
static void setBufferBit(struct Coord<KeyType>& coord) {
    coord.type_pointer = (static_cast<uint64_t>(1) << 60) | coord.type_pointer;
}

// true means it has buffer
template <class KeyType>
static bool hasBuffer(const struct Coord<KeyType>& coord) {
    return ((static_cast<uint64_t>(1) << 60) & coord.type_pointer) > 0;
}

template <class KeyType>
static void setBufferHasDataBit(struct Coord<KeyType>& coord) {
    coord.type_pointer = (static_cast<uint64_t>(1) << 59) | coord.type_pointer;
}

template <class KeyType>
static void clearBufferHasDataBit(struct Coord<KeyType>& coord) {
    coord.type_pointer = (clearBufferHasDataBitMask & coord.type_pointer);
}

// true means buffer has data
template <class KeyType>
static bool hasBufferData(const struct Coord<KeyType>& coord) {
    return ((static_cast<uint64_t>(1) << 59) & coord.type_pointer) > 0;
}

// true means need to check buffer
template <class KeyType>
static bool needCheckBuffer(const struct Coord<KeyType>& coord) {
    constexpr uint64_t mask = static_cast<uint64_t>(3) << 59;
    return (coord.type_pointer & mask) == mask;
}

}  // namespace lisa

#endif  // __LISA_BASE_H__