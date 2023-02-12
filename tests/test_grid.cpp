#include<iostream>

#include "util/containers.h"


int main() {
    using namespace std;
    using namespace boltzi::util;

    Grid<3, int> test_grid{5, 5, 5};
    cout<<"Initialized."<<endl;
    test_grid.at({1, 1, 1}) = 5;
    cout<<"Written."<<endl;
    cout<<"Read: "<<test_grid.at({1, 1, 1})<<"."<<endl;
    for(auto it = test_grid.begin(); it!= test_grid.end(); ++it)
        (*it).second= it-test_grid.begin();
    cout<<"Iterator."<<std::endl;
    for(auto it = test_grid.begin(); it!= test_grid.end(); ++it) {
        auto [idx, val] = *it;
        for(unsigned i=0;i<idx.size();++i)
            std::cout<<idx[i]<<", ";
        std::cout<<" : "<<val<<std::endl;
    }
    return 0;
}
