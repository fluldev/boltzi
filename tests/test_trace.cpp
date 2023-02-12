#include "data_tracer.h"

using namespace boltzi;

int main() {
    DataTracer<int, int ,int> tracer1("tracer1.output");
    tracer1.trace(std::make_tuple(1,2,3));
    tracer1.trace(std::make_tuple(4,5,6));
    tracer1.save();
    
    DataTracer<double, double ,double> tracer2("tracer2.output");
    tracer2.trace(std::make_tuple(1.0,2.0,3.0));
    tracer2.trace(std::make_tuple(4.0,5.0,6.0));
    tracer2.save();
    return 0;
}
