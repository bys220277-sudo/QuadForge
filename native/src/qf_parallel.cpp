// Пул потоков (fork-join) без зависимостей от OpenMP.
#include "qf_common.h"

#include <chrono>
#include <condition_variable>
#include <cstdarg>
#include <mutex>
#include <thread>

namespace qf {

namespace {

class Pool {
public:
    Pool() { resize(default_threads()); }
    ~Pool() { shutdown(); }

    static int default_threads() {
        const char* env = std::getenv("QUADFORGE_THREADS");
        if (env) {
            int v = std::atoi(env);
            if (v > 0) return v;
        }
        unsigned hc = std::thread::hardware_concurrency();
        return hc ? (int)hc : 4;
    }

    int size() const { return nthreads_; }

    void resize(int n) {
        if (n < 1) n = 1;
        if (n > 256) n = 256;
        shutdown();
        nthreads_ = n;
        stop_ = false;
        gen_ = 0;
        for (int i = 1; i < n; ++i) workers_.emplace_back([this, i] { worker_loop(i); });
    }

    // Запускает fn(thread_index) на всех потоках, включая вызывающий.
    void run(const std::function<void(int)>& fn) {
        if (nthreads_ == 1 || in_parallel_) {
            fn(0);
            return;
        }
        std::lock_guard<std::mutex> run_lock(run_mutex_);
        {
            std::lock_guard<std::mutex> lk(m_);
            job_ = &fn;
            pending_ = nthreads_ - 1;
            ++gen_;
        }
        cv_.notify_all();
        in_parallel_ = true;
        fn(0);
        {
            std::unique_lock<std::mutex> lk(m_);
            done_cv_.wait(lk, [this] { return pending_ == 0; });
            job_ = nullptr;
        }
        in_parallel_ = false;
    }

private:
    void worker_loop(int idx) {
        u64 seen = 0;
        for (;;) {
            const std::function<void(int)>* job = nullptr;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [&] { return stop_ || gen_ != seen; });
                if (stop_) return;
                seen = gen_;
                job = job_;
            }
            if (job) (*job)(idx);
            {
                std::lock_guard<std::mutex> lk(m_);
                if (--pending_ == 0) done_cv_.notify_one();
            }
        }
    }

    void shutdown() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto& t : workers_) t.join();
        workers_.clear();
    }

    int nthreads_ = 1;
    std::vector<std::thread> workers_;
    std::mutex m_, run_mutex_;
    std::condition_variable cv_, done_cv_;
    const std::function<void(int)>* job_ = nullptr;
    int pending_ = 0;
    u64 gen_ = 0;
    bool stop_ = false;
    std::atomic<bool> in_parallel_{false};
};

Pool& pool() {
    static Pool* p = new Pool();  // намеренно не уничтожается: избегаем гонок при выгрузке DLL
    return *p;
}

}  // namespace

int num_threads() { return pool().size(); }
void set_num_threads(int n) { pool().resize(n); }

void parallel_for(i64 n, const std::function<void(i64, i64)>& body, i64 grain) {
    if (n <= 0) return;
    if (grain < 1) grain = 1;
    const i64 nchunks = (n + grain - 1) / grain;
    Pool& p = pool();
    if (nchunks <= 1 || p.size() == 1) {
        for (i64 c = 0; c < nchunks; ++c) body(c * grain, std::min(n, (c + 1) * grain));
        return;
    }
    std::atomic<i64> next{0};
    p.run([&](int) {
        for (;;) {
            i64 c = next.fetch_add(1, std::memory_order_relaxed);
            if (c >= nchunks) break;
            body(c * grain, std::min(n, (c + 1) * grain));
        }
    });
}

double parallel_sum(i64 n, const std::function<double(i64, i64)>& partial, i64 grain) {
    if (n <= 0) return 0.0;
    const i64 nchunks = (n + grain - 1) / grain;
    std::vector<double> parts((size_t)nchunks, 0.0);
    parallel_for(nchunks, [&](i64 cb, i64 ce) {
        for (i64 c = cb; c < ce; ++c) parts[(size_t)c] = partial(c * grain, std::min(n, (c + 1) * grain));
    }, 1);
    double s = 0.0;
    for (double v : parts) s += v;
    return s;
}

double now_seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

void Ctx::logf(const char* fmt, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    log += buf;
    log += '\n';
    if (verbose) std::fprintf(stderr, "[qf] %s\n", buf);
}

}  // namespace qf
