#include <iostream>
#include <chrono>
#include <ratio>
#include <thread>
#include <vector>
#include "util/pool.h"

using namespace boltzi::util;

int main() {
    ParallelRunner runner(8);
    auto somejob = []{
        for(unsigned i=0;i<10;++i) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout<<i<<std::endl;
        }
    };
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.add_job(somejob);
    runner.exec();

    std::cout<<"Test container parallelization."<<std::endl;
    std::vector<int> test_v;
    for(unsigned i=0;i<200;++i) {
        test_v.push_back(i);
    } 
    runner.add_container_job(
        [](std::vector<int>::iterator it) {
            std::cout<<(*it)<<std::endl; 
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }, 
        test_v
    );
    runner.exec();
    return 0;
}
