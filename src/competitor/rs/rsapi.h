#include "../indexInterface.h"
#include "./src/builder.h"
#include "./src/radix_spline.h"

template <class KEY_TYPE, class PAYLOAD_TYPE>
class RSInterface : public indexInterface<KEY_TYPE, PAYLOAD_TYPE> {
   public:
    RSInterface() {}

    ~RSInterface() {}

    void init(Param *param = nullptr) {}

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
    rs::RadixSpline<KEY_TYPE> rs_;
    std::vector<KEY_TYPE> data_;
};

template <class KEY_TYPE, class PAYLOAD_TYPE>
void RSInterface<KEY_TYPE, PAYLOAD_TYPE>::bulk_load(
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param) {
    KEY_TYPE min_key = key_value[0].first;
    KEY_TYPE max_key = key_value[num - 1].first;
    rs::Builder<KEY_TYPE> rs_builder_(min_key, max_key, param->num_radix_bits,
                                      param->max_error);
    for (size_t i = 0; i < num; ++i) {
        KEY_TYPE key = key_value[i].first;
        rs_builder_.AddKey(key);
        data_.push_back(key);
    }
    rs_ = rs_builder_.Finalize();
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool RSInterface<KEY_TYPE, PAYLOAD_TYPE>::get(KEY_TYPE key, PAYLOAD_TYPE &val,
                                              Param *param) {
    auto res = rs_.Lookup(key, &data_);
    if (res == -1) {
        return false;
    } else {
        val = key;
        return true;
    }
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool RSInterface<KEY_TYPE, PAYLOAD_TYPE>::put(KEY_TYPE key, PAYLOAD_TYPE value,
                                              Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool RSInterface<KEY_TYPE, PAYLOAD_TYPE>::update(KEY_TYPE key,
                                                 PAYLOAD_TYPE value,
                                                 Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool RSInterface<KEY_TYPE, PAYLOAD_TYPE>::remove(KEY_TYPE key, Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
size_t RSInterface<KEY_TYPE, PAYLOAD_TYPE>::scan(
    KEY_TYPE key_low_bound, size_t key_num,
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *result, Param *param) {
    return 0;
}