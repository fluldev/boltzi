#include <iostream>
#include "util/lindisp.h"

using namespace std;
using namespace boltzi;

int main() {
    double omega_d = 80e10;
    double temp = 1;
    double v_g = 1;
    cout<<"should be 0: "<<
        util::lindisp::canonical_dist_integral<2>(v_g, 0, omega_d)
        <<endl;
    
    cout<<"should be x: "<<
        util::lindisp::canonical_dist_integral<2>(v_g, 1, omega_d)
        <<endl;

    cout<<"should be 8.29043e+43: "<<
        util::lindisp::canonical_dist_integral<3>(v_g, temp, omega_d)
        <<endl;

    temp = 10;
    cout<<"should be "<<temp<<": "<<util::lindisp::canonical_dist_pseudo_temp<3>(
        util::lindisp::canonical_dist_integral<3>(v_g, temp, omega_d),
        v_g, temp, omega_d
    )<<endl;

    temp = 1;
    double deviational_temp = .5;
    cout<<util::lindisp::canonical_dist_max<3>(temp, deviational_temp, v_g)<<endl;
    deviational_temp = .1;
    cout<<util::lindisp::canonical_dist_max<3>(temp, deviational_temp, v_g)<<endl;
    return 0;
}
