#pragma once

#include "../boltzi_types.h"
#include "../boltzi_config.h"
#include <boost/math/special_functions/lambert_w.hpp>
#include "basics.h"
#include "constants.h"

#include <boost/math/special_functions/math_fwd.hpp>
#include <cmath>
#include <limits>

// TODO: check size of n everywhere otherweise segfault


namespace boltzi{
namespace util{
namespace lindisp { 
        
// Calculation with std::riemann_zeta takes too long
const double riemann_zeta_tbl[] = {-0.5, std::numeric_limits<double>::infinity(), 1.64493, 1.20206, 1.08232, 1.03693, 1.01734, \
    1.00835, 1.00408, 1.00201, 1.00099, 1.00049, 1.00025, 1.00012, \
    1.00006, 1.00003, 1.00002, 1.00001, 1., 1., 1.
};

inline double dos(double omega, double v_g) {
    return ipow<2>(omega) / (2*ipow<2>(PI)*ipow<3>(v_g));
}

inline double total_phonon_density(double omega_debye, double v_g) {
    return ipow<3>(omega_debye) / (6*ipow<2>(PI) * ipow<3>(v_g));
}

// table starts with n=2 -> indexing with n-2
// more values can be calculated with scipts/canonical_integral_fit.py
const unsigned n_max = 15;  // max allowed n
const double tanh_fitparams[][3] = {
{1.0037481896798244, -0.2152410157282248, 0.32655869522522984},
{1.0035642715391897, -0.43421186447099613, 0.2837488763649314},
{1.004581595163449, -0.6011055661786162, 0.2540362801966956},
{1.0077550537710915, -0.7282826821685925, 0.22963192178783357},
{1.0156350947803914, -0.8163435845334377, 0.20687424286715453},
{1.013886800307991, -0.9260205413534944, 0.19438265538712168},
{1.0126252063544088, -1.0250047502218136, 0.18362968480788708},
{1.0116802924923554, -1.1155978256344652, 0.17426008482098573},
{1.0109527070933384, -1.1993635905304123, 0.16600903197317435},
{1.0103799762162806, -1.2774111794268252, 0.1586751918486879},
{1.0099207798839385, -1.3505609954626792, 0.1521031011165501},
{1.0095467835760137, -1.4194430012182737, 0.14617101853992365},
{1.009238012298972, -1.4845563224173859, 0.14078234691541183},
{1.008980060792296, -1.5463065680410892, 0.13585948659023667},
};

const double deviational_max_fitparams[][3] = {
    {-4.41991345, 5.05155082, 2.10369747},
    {-30.22823601, 31.64157007, 3.1280662},
    {-294.28044079, 299.04982306, 4.15358734},
    {-3649.16948733, 3670.34672258, 5.16741721}
};

template<unsigned int n>
inline double canonical_dist(double omega, double temp, double deviational_temp, double v_g) {
    return  1/(2*ipow<2>(PI)*ipow<3>(v_g)) * ipow<n>(omega) 
            * std::abs(
                    1/(std::exp(HBAR*omega / (BMC*temp))-1) 
                    - 1/(std::exp(HBAR*omega / (BMC*deviational_temp))-1)
            ); 
}

inline double canonical_dist(unsigned n, double omega, double temp, double deviational_temp, double v_g) {
    return  1/(2*ipow<2>(PI)*ipow<3>(v_g)) * ipow(omega, n) 
            * std::abs(
                    1/(std::exp(HBAR*omega / (BMC*temp))-1) 
                    - 1/(std::exp(HBAR*omega / (BMC*deviational_temp))-1)
            ); 
}

template<unsigned int n>
inline std::function<double(double)> canonical_dist_factory(double temp, double deviational_temp, double v_g) {
    return [temp, deviational_temp, v_g] (double omega) -> double {
        return canonical_dist<n>(omega, temp, deviational_temp, v_g);
    };
}

inline std::function<double(double)> canonical_dist_factory(unsigned n, double temp, double deviational_temp, double v_g) {
    return [temp, deviational_temp, v_g, n] (double omega) -> double {
        return canonical_dist(n, omega, temp, deviational_temp, v_g);
    };
}

template<unsigned int n>
inline double canonical_dist_max(double temp, double deviational_temp, double v_g) {
    static_assert(n>=2 && n<6, "Unsupported wavevector power.");
    using boost::math::lambert_w0;

    double a = deviational_max_fitparams[n-2][0];
    double b = deviational_max_fitparams[n-2][1];
    double c = deviational_max_fitparams[n-2][2];
    double temp_max = std::max(temp, deviational_temp);
    double temp_min = std::min(temp, deviational_temp);
    if(temp_min==0 || temp_max/temp_min>3)
        return canonical_dist<n>(
                temp_max*BMC/HBAR*(lambert_w0(-(int)n*std::exp(-(int)n))+n), 
                temp, deviational_temp, v_g);
    return  ipow<n>(temp_max*BMC/HBAR)*(a+b*(1-std::exp(-c*temp_max/temp_min)))
            /(2*ipow<2>(PI)*ipow<3>(v_g));
}

inline double canonical_dist_max(unsigned n, double temp, double deviational_temp, double v_g) {
    using boost::math::lambert_w0;

    double a = deviational_max_fitparams[n-2][0];
    double b = deviational_max_fitparams[n-2][1];
    double c = deviational_max_fitparams[n-2][2];
    double temp_max = std::max(temp, deviational_temp);
    double temp_min = std::min(temp, deviational_temp);
    if(temp_min==0 || temp_max/temp_min>3)
        return canonical_dist(
                n,
                temp_max*BMC/HBAR*(lambert_w0(-(int)n*std::exp(-(int)n))+n), 
                temp, deviational_temp, v_g);
    return  ipow(temp_max*BMC/HBAR, n)*(a+b*(1-std::exp(-c*temp_max/temp_min)))
            /(2*ipow<2>(PI)*ipow<3>(v_g));
}

template<unsigned int n>
inline double canonical_dist_pseudo_temp(double s, double v_g, double temp, double omega_d) {
    // without constant factors
    static_assert(n>=2 && n<6, "Unsupported wavevector power.");
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 
    
    // since delta must be smaller than 1 only corrections to higher temperatures are possible (energy density gets larger with higher temperatures)
    // therefore if the temperature guess yields an delta >= 1 we can just use the initial guess
    if(a*std::tanh(b+c*HBAR*omega_d/(BMC*temp))>=1) return temp;

    double res = euler_method(
            [&a,&b,&c,&omega_d,&s,&v_g](double x)->double {
                return  a*std::tanh(b+x*c)/ipow<n+1>(x) 
                        - s*(2*ipow<2>(PI)*ipow<3>(v_g)) / (factorial(n)*riemann_zeta_tbl[n+1]*ipow<n+1>(omega_d));
            }, 
            [&a,&b,&c,&omega_d](double x)->double {
                return  a/ipow<2>(std::cosh(b+x*c))*c/ipow<n+1>(x) 
                        - a*std::tanh(b+x*c)/ipow<n+2>(x)*(n+1);
            }, 
            HBAR*omega_d/temp/BMC
    );
    return HBAR*omega_d/(res*BMC);
}

inline double canonical_dist_pseudo_temp(unsigned n,double s, double v_g, double temp, double omega_d) {
    // without constant factors
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 
    
    // since delta must be smaller than 1 only corrections to higher temperatures are possible (energy density gets larger with higher temperatures)
    // therefore if the temperature guess yields an delta >= 1 the initial guess suffices (if calculated with the omega_debye=infinity formula)
    if(a*std::tanh(b+c*HBAR*omega_d/(BMC*temp))>=1) return temp;

    double res = euler_method(
            [a,b,c,omega_d,s,v_g,n](double x)->double {
                return  a*std::tanh(b+x*c)/ipow(x,n+1) 
                        - s*(2*ipow<2>(PI)*ipow<3>(v_g)) / (factorial(n)*riemann_zeta_tbl[n+1]*ipow(omega_d,n+1));
            }, 
            [a,b,c,omega_d,n](double x)->double {
                return  a/ipow<2>(std::cosh(b+x*c))*c/ipow(x,n+1) 
                        - a*std::tanh(b+x*c)/ipow(x,n+2)*(n+1);
            }, 
            HBAR*omega_d/temp/BMC
    );
    return HBAR*omega_d/(res*BMC);
}

inline double canonical_dist_temp(double s, double v_g, double omega_d) {
    double temp_guess = std::pow(s/(factorial(3)*riemann_zeta_tbl[4])*2*ipow<2>(PI)*ipow<3>(v_g)*ipow<4>(HBAR/BMC),1./4);
    return canonical_dist_pseudo_temp<3>(s, v_g, temp_guess, omega_d);
}

inline double branch_weight(double omega_debye, double v_g) {
    return ipow<3>(omega_debye/v_g);
}

template<unsigned int n>
inline double canonical_dist_integral(double v_g, double temp, double omega_d) {
    // without constant factor
    static_assert(n>=2 && n<6, "Unsupported wavevector power.");
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 

    double delta = std::min(a*std::tanh(b+c*HBAR/BMC*omega_d/temp), 1.0); 
    return  factorial(n)*riemann_zeta_tbl[n+1] * delta 
            / (2*ipow<2>(PI)*ipow<3>(v_g)) * ipow<n+1>(BMC*temp/HBAR);
}

inline double canonical_dist_integral(unsigned n, double v_g, double temp, double omega_d) {
    // without constant factor
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 

    double delta = std::min(a*std::tanh(b+c*HBAR/BMC*omega_d/temp), 1.0); 
    return  factorial(n)*riemann_zeta_tbl[n+1] * delta 
            / (2*ipow<2>(PI)*ipow<3>(v_g)) * ipow(BMC*temp/HBAR,n+1);
}

template<unsigned int n>
inline double total_spawnrate_spatial(double gradient, double v_g, double temp, double omega_d) {
    // without constant factor
    static_assert(n>=2 && n<6, "Unsupported wavevector power.");
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 

    double delta = std::min(a*std::tanh(b+c*HBAR/BMC*omega_d/temp), 1.0); 
    return  std::abs(gradient*ipow<n>(temp) / (4*ipow<2>(PI*v_g)) * ipow<n+1>(BMC/HBAR)
            * (
                factorial(n+1)*riemann_zeta_tbl[n+1] * delta
                - ipow<n+1>(HBAR*omega_d/(BMC*temp))/(std::exp(HBAR*omega_d/(BMC*temp))-1)
            ));
}

inline double total_spawnrate_spatial(unsigned n, double gradient, double v_g, double temp, double omega_d) {
    // without constant factor
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 

    double delta = std::min(a*std::tanh(b+c*HBAR/BMC*omega_d/temp), 1.0); 
    return  std::abs(gradient*ipow(temp,n) / (4*ipow<2>(PI*v_g)) * ipow(BMC/HBAR,n+1)
            * (
                factorial(n+1)*riemann_zeta_tbl[n+1] * delta
                - ipow(HBAR*omega_d/(BMC*temp),n+1)/(std::exp(HBAR*omega_d/(BMC*temp))-1)
            ));
}

template<unsigned int n>
inline double spawnrate_dist_spatial(double omega, double temp) {
    static_assert(n==2, "n!=2 not supported.");
    return  ipow<3>(HBAR*omega / (BMC*temp)) 
            * std::exp(HBAR*omega/(BMC*temp)) / ipow<2>(std::exp(HBAR*omega/(BMC*temp))-1);
}

template<unsigned int n>
inline double spawnrate_dist_spatial_max(double temp) {
    static_assert(n==2, "n!= not supported.");
    return spawnrate_dist_spatial<n>(2.57567891 * BMC*temp/HBAR, temp);
}

inline double total_integral_temp_derivative(unsigned n, double temp, double omega_d) {
    double a = tanh_fitparams[n-2][0]; 
    double b = tanh_fitparams[n-2][1]; 
    double c = tanh_fitparams[n-2][2]; 
    
    double delta = std::min(a*std::tanh(b+c*HBAR/BMC*omega_d/temp), 1.0); 
    return  factorial(n) * riemann_zeta_tbl[n] * delta
            - ipow(HBAR/BMC*omega_d/temp, n) / (std::exp(HBAR/BMC*omega_d/temp)-1);
}

inline double linearized_post_scatter_dist(unsigned n, double omega, double equ_temp) {
    // function omits constant factors
    double x = HBAR*omega/(BMC * equ_temp);
    return ipow(omega,n+1) * std::exp(x) / ipow<2>(std::exp(x)-1);
}

// starting from n=3
const double linearized_post_scatter_max_tbl[] = {
    2.57568, 3.83002, 4.92812, 5.96941, 6.98708, 7.99461, 8.99777, 9.99909
};

inline double linearized_post_scatter_dist_max(unsigned n, double equ_temp) {
    //TODO: check if n is in range [3,10]
    double omega_max = linearized_post_scatter_max_tbl[n-2] * BMC * equ_temp / HBAR; 
    return linearized_post_scatter_dist(n, omega_max, equ_temp);
}
}}}
