#include <iomanip>

#pragma once

struct Param { // for xindex
  size_t worker_num;
  uint32_t thread_id;
  std::string config_file; // for libox
  
  // hope parameters
  size_t node_capacity=1000; // for hope
  double top_k=0.05;
  // lisa parameters
  size_t num_radix_bits=18;
  size_t max_error=32;
  size_t alpha=10; // min line length

  Param(size_t worker_num, uint32_t thread_id, int node_capacity, std::string config="") : worker_num(worker_num), thread_id(thread_id), node_capacity(node_capacity), config_file(config) {}
};

struct BaseCompare {
  template<class T1, class T2>
  bool operator()(const T1 &x, const T2 &y) const {
    static_assert(
      std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
      "Comparison types must be numeric.");
    return x < y;
  }
};

template<class KEY_TYPE, class PAYLOAD_TYPE, class KeyComparator=BaseCompare>
class indexInterface {
public:
  virtual void bulk_load(std::pair <KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param = nullptr) = 0;

  virtual bool get(KEY_TYPE key, PAYLOAD_TYPE &val, Param *param = nullptr) = 0;

  virtual bool put(KEY_TYPE key, PAYLOAD_TYPE value, Param *param = nullptr) = 0;

  virtual bool update(KEY_TYPE key, PAYLOAD_TYPE value, Param *param = nullptr) = 0;

  virtual bool remove(KEY_TYPE key, Param *param = nullptr) = 0;

  virtual size_t scan(KEY_TYPE key_low_bound, size_t key_num, std::pair<KEY_TYPE, PAYLOAD_TYPE> *result,
                      Param *param = nullptr) = 0;

  virtual void init(Param *param = nullptr) = 0;

  virtual long long memory_consumption() = 0; // bytes
};