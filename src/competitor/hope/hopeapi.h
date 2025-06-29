#include "../indexInterface.h"
#include "./src/hope.h"

template <class KEY_TYPE, class PAYLOAD_TYPE>
class HopeInterface : public indexInterface<KEY_TYPE, PAYLOAD_TYPE> {
   public:
    HopeInterface() {}

    ~HopeInterface() {}

    void init(Param *param = nullptr) {
        hope_.SetParameters(param->node_capacity, param->top_k,
                            param->max_error, param->alpha, param->use_radix, param->temp_node_cap);
    }

    void bulk_load(std::pair<KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num,
                   Param *param = nullptr);

    bool get(KEY_TYPE key, PAYLOAD_TYPE &val, Param *param = nullptr);

    bool put(KEY_TYPE key, PAYLOAD_TYPE value, Param *param = nullptr);

    bool update(KEY_TYPE key, PAYLOAD_TYPE value, Param *param = nullptr);

    bool remove(KEY_TYPE key, Param *param = nullptr);

    size_t scan(KEY_TYPE key_low_bound, size_t key_num,
                std::pair<KEY_TYPE, PAYLOAD_TYPE> *result,
                Param *param = nullptr);

    long long memory_consumption() { return 0; }

   private:
    hopens::Hope<KEY_TYPE, PAYLOAD_TYPE> hope_;
};

template <class KEY_TYPE, class PAYLOAD_TYPE>
void HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::bulk_load(
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param) {
    hope_.BulkLoad(key_value, num);
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::get(KEY_TYPE key, PAYLOAD_TYPE &val,
                                                Param *param) {
    auto res = hope_.Lookup(key, val);
    return res;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::put(KEY_TYPE key,
                                                PAYLOAD_TYPE value,
                                                Param *param) {
    auto res = hope_.Insert(key, value);
    return res;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::update(KEY_TYPE key,
                                                   PAYLOAD_TYPE value,
                                                   Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::remove(KEY_TYPE key, Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
size_t HopeInterface<KEY_TYPE, PAYLOAD_TYPE>::scan(
    KEY_TYPE key_low_bound, size_t key_num,
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *result, Param *param) {
    return 0;
}