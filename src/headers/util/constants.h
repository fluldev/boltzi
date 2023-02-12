#pragma once

#include <cmath>

namespace boltzi{
namespace util {
    // constants
    const double PI = std::atan(1)*4;
    const double HBAR = 1.05457182e-34;  // reduced Planck constant m^2 kg s^-1
    const double BMC = 1.380649e-23;  // Boltzmann constant m^2 kg s^-2 K^-1
    const double lambertw[] = {  // scipy.special.lambertw(x) for x in range(0, 10)
        0.0,
        0.5671432904097838,
        0.8526055020137254,
        1.04990889496404,
        1.2021678731970429,
        1.3267246652422002,
        1.4324047758983003,
        1.5243452049841444,
        1.6058119963201776,
        1.6790164197855981
    };
}}
