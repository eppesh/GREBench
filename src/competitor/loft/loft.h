#include <vector>

#include "./src/LOFT_impl.h"
#include "./src/util.h"
#include "./src/LOFT.h"
#include "../indexInterface.h"

namespace loft {
    template<class KEY_TYPE>
    class Key {
        typedef std::array<double, 1> model_key_t;

    public:
        static constexpr size_t

        model_key_size() { return 1; }

        static Key max() {
            static Key max_key(std::numeric_limits<KEY_TYPE>::max());
            return max_key;
        }

        static Key min() {
            static Key min_key(std::numeric_limits<KEY_TYPE>::min());
            return min_key;
        }

        Key() : key(0) {}

        Key(KEY_TYPE key) : key(key) {}

        Key(const Key<KEY_TYPE> &other) { key = other.key; }

        Key &operator=(const Key<KEY_TYPE> &other) {
            key = other.key;
            return *this;
        }

        model_key_t to_model_key() const {
            model_key_t model_key;
            model_key[0] = key;
            return model_key;
        }

        friend bool operator<(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key < r.key; }

        friend bool operator>(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key > r.key; }

        friend bool operator>=(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key >= r.key; }

        friend bool operator<=(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key <= r.key; }

        friend bool operator==(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key == r.key; }

        friend bool operator!=(const Key<KEY_TYPE> &l, const Key<KEY_TYPE> &r) { return l.key != r.key; }

        KEY_TYPE key;
    } PACKED;
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
class LOFTInterface : public indexInterface<KEY_TYPE, PAYLOAD_TYPE> {
public:
    void bulk_load(std::pair <KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param);

    bool get(KEY_TYPE key, PAYLOAD_TYPE &val, Param *param);

    bool put(KEY_TYPE key, PAYLOAD_TYPE value, Param *param);

    bool update(KEY_TYPE key, PAYLOAD_TYPE value, Param *param = nullptr);

    bool remove(KEY_TYPE key, Param *param);

    size_t
    scan(KEY_TYPE key_low_bound, size_t key_num, std::pair <KEY_TYPE, PAYLOAD_TYPE> *result, Param *param);

    void init(Param *param);

    long long memory_consumption() {
        return 0;
    }

    ~LOFTInterface() {
        delete this->index;
    }

    LOFTInterface() {}

    size_t worker_num;
    size_t bg_n;
private :
    // loft::LOFT <loft::Key<KEY_TYPE>, PAYLOAD_TYPE> *index;
    loft::LOFT <KEY_TYPE, PAYLOAD_TYPE> *index;
    size_t core_num;
};

template<class KEY_TYPE, class PAYLOAD_TYPE>
void LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::init(Param *param) {
    worker_num = param->worker_num;
    bg_n = param->worker_num / 12 + 1;
//    core_num = param->worker_num;
//    if(core_num == 1) {
//        worker_num = 1;
//        bg_n = 1;
//        return;
//    }
//    worker_num = core_num - (core_num / 12 + 1 - (core_num % 12 == 0));
//    bg_n = core_num - worker_num;
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
void LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::bulk_load(std::pair <KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num,
                                                        Param *param) {
    // std::vector <loft::Key<KEY_TYPE>> key;
    std::vector <KEY_TYPE> key;
    std::vector <PAYLOAD_TYPE> value;
    bool random_insert = false;

    KEY_TYPE last_key;
    for (int i = 0; i < num; i++) {
        // if (i != 0 && key_value[i].first == last_key) { continue; }
        // key.push_back(loft::Key<KEY_TYPE>(key_value[i].first));
        key.push_back(key_value[i].first);
        value.push_back(key_value[i].second);
        // last_key = key_value[i].first;
    }

    printf("worker_num: %llu, bg_n: %llu\n", worker_num, bg_n);
    // index = new loft::LOFT<loft::Key<KEY_TYPE>, PAYLOAD_TYPE>(key, value, worker_num, bg_n);
    index = new loft::LOFT<KEY_TYPE, PAYLOAD_TYPE>(key, value, worker_num, bg_n);

    // for (int i = num / 2; i < num; i++) {
    //     this->put(key_value[i].first, key_value[i].second, param);
    // }

}

template<class KEY_TYPE, class PAYLOAD_TYPE>
bool LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::get(KEY_TYPE key, PAYLOAD_TYPE &val, Param *param) {
    // auto ret = index->query(loft::Key<KEY_TYPE>(key), val, param->thread_id);
    uint8_t tid = static_cast<uint8_t>(param->thread_id);
    auto ret = index->query(key, val, tid);
    return ret;
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
bool LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::put(KEY_TYPE key, PAYLOAD_TYPE value, Param *param) {
    // return index->insert(loft::Key<KEY_TYPE>(key), value, param->thread_id);
    uint8_t tid = static_cast<uint8_t>(param->thread_id);
    return index->insert(key, value, tid);
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
bool LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::update(KEY_TYPE key, PAYLOAD_TYPE value, Param *param) {
    // return index->update(loft::Key<KEY_TYPE>(key), value, param->thread_id);
    uint8_t tid = static_cast<uint8_t>(param->thread_id);
    return index->update(key, value, tid);
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
bool LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::remove(KEY_TYPE key, Param *param) {
    // return index->remove(loft::Key<KEY_TYPE>(key), param->thread_id);
    uint8_t tid = static_cast<uint8_t>(param->thread_id);
    return index->remove(key, tid);
}

template<class KEY_TYPE, class PAYLOAD_TYPE>
size_t LOFTInterface<KEY_TYPE, PAYLOAD_TYPE>::scan(KEY_TYPE key_low_bound, size_t key_num,
                                                     std::pair <KEY_TYPE, PAYLOAD_TYPE> *result,
                                                     Param *param) {
    // std::vector <std::pair<loft::Key<KEY_TYPE>, PAYLOAD_TYPE>> result_temp;
    // auto scan_num = index->scan(loft::Key<KEY_TYPE>(key_low_bound), key_num, result_temp, param->thread_id);
    std::vector <std::pair<KEY_TYPE, PAYLOAD_TYPE>> result_temp;
    uint8_t tid = static_cast<uint8_t>(param->thread_id);
    auto scan_num = index->scan(key_low_bound, key_num, result_temp, tid);
    PAYLOAD_TYPE accumulator = 0;
    for (auto i : result_temp) {
        accumulator += i.second;
    }
    return scan_num;
}