#include "../indexInterface.h"
#include "./src/libox.h"

template <class KEY_TYPE, class PAYLOAD_TYPE>
class LiBoxInterface : public indexInterface<KEY_TYPE, PAYLOAD_TYPE> {
   public:
    LiBoxInterface() {
        vector<int> taskCoreIDs = {
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23};  // Adjust this as needed
        for (size_t i = 0; i < taskCoreIDs.size(); i++) {
            task_pool_.emplace_back([this]() {
                std::packaged_task<void()> task;
                while (split_task_queue_.popTask(task)) {
                    task();
                }
            });
            setThreadAffinity(task_pool_.back(), taskCoreIDs[i]);
        }
    }

    ~LiBoxInterface() {
        split_task_queue_.shutdown();
        for (auto &thr : task_pool_) {
            if (thr.joinable()) thr.join();
        }
    }

    void init(Param *param = nullptr) {
        libox_.init(0.5, 0.1, &split_task_queue_, param->worker_num);
        libox_.loadConfigByFile(param->config_file);
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

    long long memory_consumption() { return libox_.get_total_size(); }

   private:
    liboxns::LiBox<KEY_TYPE, PAYLOAD_TYPE> libox_;
    TaskQueue split_task_queue_;
    vector<std::thread> task_pool_;
};

template <class KEY_TYPE, class PAYLOAD_TYPE>
void LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::bulk_load(
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *key_value, size_t num, Param *param) {
    libox_.bulk_load(key_value, num);
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::get(KEY_TYPE key,
                                                 PAYLOAD_TYPE &val,
                                                 Param *param) {
    auto res = libox_.searchKey(key);
    val = res.value;
    return (res.status == liboxns::SearchStatus::SUCCESS);
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::put(KEY_TYPE key,
                                                 PAYLOAD_TYPE value,
                                                 Param *param) {
    auto res = libox_.insertKeyValue(key, value);
    return (res.status == liboxns::InsertStatus::SUCCESS);
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::update(KEY_TYPE key,
                                                    PAYLOAD_TYPE value,
                                                    Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
bool LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::remove(KEY_TYPE key,
                                                    Param *param) {
    return false;
}

template <class KEY_TYPE, class PAYLOAD_TYPE>
size_t LiBoxInterface<KEY_TYPE, PAYLOAD_TYPE>::scan(
    KEY_TYPE key_low_bound, size_t key_num,
    std::pair<KEY_TYPE, PAYLOAD_TYPE> *result, Param *param) {
    return 0;
}