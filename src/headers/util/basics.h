#pragma once

#include<cmath>

#include "../boltzi_types.h"
#include "../boltzi_config.h"
#include "constants.h"

namespace boltzi {
namespace util{

    inline vec3 sphere_to_cart(vec2 sphere) {
        return vec3{
            std::sin(sphere[1])*std::cos(sphere[0]),    
            std::sin(sphere[1])*std::sin(sphere[0]),    
            std::cos(sphere[1])
        }; 
    }

    template<unsigned int power> 
    constexpr inline double ipow(double x) {
        if constexpr(power == 0) return 1;
        double res = 1;
        // hopefully gets unwrapped for small powers
        for(unsigned int i=0;i<power;++i) res*=x;
        return res;
    }
    
    constexpr inline double ipow(double x, unsigned power) {
        if (power == 0) return 1;
        double res = 1;
        for(unsigned int i=0;i<power;++i) res*=x;
        return res;
    }

    inline double bose_einst_dist(double omega, double temp) {
        return 1/(std::exp(HBAR*omega / (BMC*temp))-1);
    }

    inline double euler_method(std::function<double(double)> fn, std::function<double(double)> deriv, double guess) {
        double x = guess;
        double lx;
        do {
            lx = x;
            x -= fn(x) / deriv(x) * config::EULER_STEP;
        } while(std::abs((x-lx)/guess)>config::EULER_DELTA);
        return x;
    }

    constexpr unsigned factorial(unsigned n) {
        unsigned res = 1;
        for(;n>0;--n) res*=n;
        return res;
    }

    inline double modf(double val) {
        return val-static_cast<int>(val);
    }
}}

