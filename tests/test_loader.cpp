#include <iostream>

#include "sample_geometry_factory.h"

using namespace boltzi;
int main() {
    auto res = load_geometry("../test_loader.obj");
    std::cout<<"Succesfully loaded sample geometry."<<std::endl;
    std::cout<<"Box:\n"<<res.box.first<<std::endl<<res.box.second; 
    std::cout<<"Volume: "<<res.volume()<<std::endl; 
    std::cout<<"Area: "<<res.area()<<std::endl; 
    return 0;
}
