#include "headers/boltzi_config.h"
#include "headers/util/pool.h"
#include <thread>


namespace boltzi {
namespace util {
ParallelRunner pool(config::CPUS);
void ParallelRunner::Worker::loop(Worker* self) {
    while(self->active) {
        self->owner->dibs.lock();
        if(self->owner->jobs.empty()) {
            self->owner->dibs.unlock();
            continue;
        }
        self->job = self->owner->jobs.front();
        self->owner->jobs.pop();
        self->running = true;
        self->owner->dibs.unlock();
        self->job();
        self->running = false;
    }
    self->owner->dibs.unlock();
}
void ParallelRunner::add_job(std::function<void()> job) {
    if(terminated) throw BoltziException("Pool's closed.");
    jobs.push(job);
}

void ParallelRunner::terminate() {
    if(terminated) return;
    for(auto& worker : workers) worker.deactivate();
    dibs.unlock();
    for(auto& thread : threads) thread.join();
    terminated=true;
    workers.clear();
    threads.clear();
}

void ParallelRunner::exec() {
    if(jobs.empty()) return;
    if(terminated) throw BoltziException("Pool's closed.");
    try {
        dibs.unlock();  // lets threads dequeue jobs until the queue is empty 
        while(true){
            dibs.lock();
            if(jobs.empty()) break;
            dibs.unlock();
            std::this_thread::yield();
        }  // wait until all jobs are finished without blocking
    }
    catch(...) {
        dibs.lock();  // lock must be held by the pool outside of exec 
        throw;
    }
    for(auto& worker : workers) 
        while(worker.running) std::this_thread::yield();
}
}}
