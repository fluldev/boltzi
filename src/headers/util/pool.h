#pragma once

#include <cstddef>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>

#include <iostream>

#include "../boltzi_exception.h"

// Ein Mutex pro pool.
//
// Falls Job-Queue leer:
// Mutex wird versucht von allen Workern und dem Pool selbst versucht zu locken.
// Worker releasen das Mutex wenn sie merken, dass die Job-Queue leer ist.
// Pool released erst, wenn die Queue nicht mehr leer ist.
// Dadurch warten alle Threads ohne CPU-Ressourcen zu verbrauchen und der 
// Pool besitzt das Mutex.
//
// Falls Job-Queue nicht leer:
// Pool versucht nicht mehr das Mutex zu locken.
// Worker, der das Mutex besitzt dequeued einen Job von der Job-Queue und unlocked das
// Mutex wieder (so dass andere Threads auch Jobs bekommen)
// Dieser Prozess läuft ab, bis keine Jobs mehr vorhanden sind, dann setzt der 
// Job-Queue leer State wieder ein.


namespace boltzi{
namespace util{
class ParallelRunner{
    private:
    std::mutex dibs;
    std::queue<std::function<void()>> jobs; 
    volatile bool terminated = false;
    unsigned cpus;

    class Worker {
        private:
        std::function<void()> job;
        ParallelRunner* owner;
        volatile bool active = true;

        public:
        volatile bool running = false;
        void deactivate() {active=false;}
        Worker() = delete;
        Worker(ParallelRunner* owner) : owner(owner) {}
        static void loop(Worker* self);
    };
    std::vector<Worker> workers;
    std::vector<std::thread> threads;

    public:
    unsigned n_jobs() const {
        return jobs.size();
    }
    ParallelRunner() = delete;
    ParallelRunner(unsigned cpus) : cpus(cpus) {
        workers.reserve(cpus);
        for(unsigned i = 0; i<cpus; ++i) 
            workers.push_back(Worker(this));

        dibs.lock();  // important! otherwise the threads are locking, unlocking 
                      // dibs continuously until the pool receives a job
                      
        threads.reserve(cpus);
        for(auto& worker : workers)
            threads.push_back(std::thread(Worker::loop, &worker));
    }
    void add_job(std::function<void()> job);

    void terminate();

    void exec();
    
    template<typename T>
    void add_container_job(std::function<void(typename T::iterator)> job, T& work_domain) {
        size_t elements_per_worker = work_domain.size() / cpus;
        if((work_domain.size() % cpus) > 0) elements_per_worker+=1;
        for(unsigned i=0; i<cpus; ++i) {
            const auto beg = work_domain.begin() + std::min(elements_per_worker*i, work_domain.size());
            const auto end = work_domain.begin() + std::min(elements_per_worker*(i+1), work_domain.size());
            if(beg == work_domain.end()) break;
            add_job(
                [beg,end,job]()->void{
                    for(auto it=beg;it!=end;++it) job(it); 
                }
            );
        }
    }
    template<typename T>
    void add_container_job(std::function<void(typename T::iterator)> job, T& work_domain, unsigned elements_per_worker) {
        auto job_cursor = work_domain.begin();
        unsigned n_added = 0;
        while(job_cursor!=work_domain.end()) {
            auto end = work_domain.end();
            if((job_cursor-work_domain.begin() + elements_per_worker) < work_domain.size())
                end = job_cursor+elements_per_worker;
            add_job(
                [job_cursor,end,job]()->void{
                    for(auto it=job_cursor;it!=end;++it) job(it); 
                }
            );
            n_added++;
            job_cursor = end;
        }
    }

    ~ParallelRunner() {terminate();}
}; 

extern ParallelRunner pool;
}
}
