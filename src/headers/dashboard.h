#pragma once
#include <iostream>
#include <graphs.hpp>
#include <algorithm>
#include <vector>

namespace boltzi{
void dashboard(
        double t, 
        const std::vector<std::vector<double>>& n_phonons, 
        unsigned n_steps,
        unsigned n_total_steps,
        const std::vector<std::vector<double>>& mean_temperature
) {
    double max_val = 0;
    for(auto& ele : n_phonons)
        if(ele[1]>max_val) max_val = ele[1];
    double max_temp = 0;
    for(auto& ele : mean_temperature)
        if(ele[1]>max_temp) max_temp = ele[1];
    graphs::array(80,150, -1, t, -1, max_val*1.2, n_phonons);
    graphs::array(80,150, -1, t, -1, max_temp*1.2, mean_temperature);
}
}
