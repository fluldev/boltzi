#include <iostream>
#include <vector>

#include "boltzi_util.h"
#include <armadillo>


int main() {
    using namespace std;
    using namespace boltzi;
    using namespace arma;

    vector<int> for_draw{1, 2, 3, 4, 5};

    cout<<util::random.draw_real()<<endl;
    cout<<util::random.draw_uint(5)<<endl;
    cout<<util::random.draw_element(for_draw)<<endl;
    cout<<util::random.draw_index_weighted(for_draw)<<endl;
    util::random.random_dir().print();
    util::random.random_dir_area(vec3{0,0,1}, vec3{1,0,0}, vec3{0,1,0}).print();
    
    util::BinDist bd{[] (double x) -> double {return x;}, 0, 1, 10};
    cout<<bd()<<endl;
    return 0;
}

