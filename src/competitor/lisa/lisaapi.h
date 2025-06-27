#include "../indexInterface.h"
#include "./src/lisa.h"

template <class KEY_TYPE, class PAYLOAD_TYPE>
class LisaInterface : public indexInterface<KEY_TYPE, PAYLOAD_TYPE> {
   public:
    LisaInterface() {
        std::pair<KEY_TYPE, PAYLOAD_TYPE> min_key = std::make_pair(0, 0);
        std::pair<KEY_TYPE, PAYLOAD_TYPE> max_key = std::make_pair(0, 0);
        lisa_ = new lisa::LISA<KEY_TYPE, PAYLOAD_TYPE>(min_key, max_key, 18, 32,
                                                       10);
    }

    ~LisaInterface() {
        if (lisa_) {
            delete lisa_;
            lisa_ = nullptr;
        }
    }

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
    lisa::LISA<KEY_TYPE, PAYLOAD_TYPE> *lisa_;
};

template <class KEY_TYPE, class PAYLOAD_TYPE>
void LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::bulk_load(
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param) {
    lisa_->SetParameters(key_value[0], key_value[num - 1],
                         param->num_radix_bits, param->max_error, param->alpha);
    lisa_->bulk_load(key_value, num);
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::get(KEY_TYPE key, PAYLOAD_TYPE &val,
                                                Param *param) {
    auto res = lisa_->Lookup(key, val);
    return res;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::put(KEY_TYPE key,
                                                PAYLOAD_TYPE value,
                                                Param *param) {
    auto res = lisa_->Insert(std::make_pair(key, value));
    return res;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::update(KEY_TYPE key,
                                                   PAYLOAD_TYPE value,
                                                   Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::remove(KEY_TYPE key, Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
size_t LisaInterface<KEY_TYPE, PAYLOAD_TYPE>::scan(
    KEY_TYPE key_low_bound, size_t key_num,
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *result, Param *param) {
    return 0;
}