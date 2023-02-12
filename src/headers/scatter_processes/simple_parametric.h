#pragma once

#include <armadillo>
#include <boost/math/tools/config.hpp>
#include <cmath>

#include "../boltzi_util.h"
#include "../util/lindisp.h"
#include "../material.h"

namespace boltzi {
namespace scatter {
    class SimpleParametricTwoPhononProcess : public TwoPhononProcess {
        public:
        const unsigned int temp_pow;
        const unsigned int k_pow;
        const double fit_constant;
        const double v_g;
        const double omega_debye;

        SimpleParametricTwoPhononProcess() = delete;
        SimpleParametricTwoPhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye) 
            : temp_pow(temp_pow), k_pow(k_pow), fit_constant(fit_constant), v_g(v_g), omega_debye(omega_debye) {}
       
        double scatter_rate(const PhononData& phonon, double temp) const override {
            return fit_constant * util::ipow(arma::norm(phonon.k),k_pow) * util::ipow(temp,temp_pow);
        }
    };
    
    class SimpleParametricThreePhononProcess : public ThreePhononProcess {
        public:
        const unsigned int temp_pow;
        const unsigned int k_pow;
        const double fit_constant;
        const double v_g;
        const double omega_debye;
        const bool is_umklapp;

        SimpleParametricThreePhononProcess() = delete;
        SimpleParametricThreePhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye, bool is_umklapp) 
            : temp_pow(temp_pow), k_pow(k_pow), fit_constant(fit_constant), v_g(v_g), omega_debye(omega_debye), is_umklapp(is_umklapp) {}
       
        double scatter_rate(const PhononData& phonon, double temp) const override {
            return fit_constant * util::ipow(arma::norm(phonon.k),k_pow) * util::ipow(temp,temp_pow);
        }

        double pseudo_temp(double scatter_sum, double temp) const override {
            scatter_sum = std::max(scatter_sum, 0.);
            scatter_sum *= util::ipow(v_g, k_pow);
            scatter_sum /= fit_constant * util::ipow(temp, temp_pow);
            double temp_guess = std::pow(
                    scatter_sum/(util::factorial(3+k_pow)*util::lindisp::riemann_zeta_tbl[4+k_pow])\
                    *2*util::ipow<2>(util::PI)*util::ipow<3>(v_g)*util::ipow(util::HBAR/util::BMC,4+k_pow),1./(4+k_pow));
            return util::lindisp::canonical_dist_pseudo_temp(
                k_pow+3,
                scatter_sum, v_g, temp_guess, omega_debye
            );
        }
        
        double scattered_omegasum(double scatter_temp, double temp) const override {
            return util::lindisp::canonical_dist_integral(3+k_pow, v_g, temp, omega_debye) / util::ipow(v_g, k_pow) \
                * fit_constant * util::ipow(scatter_temp,temp_pow);
        } 
        
        vec3 draw_k(PhononData& phonon, double scatter_temp, double pseudo_temp, double deviational_temp) const override {
            // scatter_temp not needed here since it is eliminated by the normalisation of the distribution
            #ifndef BOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING
            return (is_umklapp?util::random.random_dir():arma::normalise(phonon.k)) * util::random.accept_reject(
                util::lindisp::canonical_dist_factory(3+k_pow, pseudo_temp, deviational_temp, v_g), 
                0, omega_debye, util::lindisp::canonical_dist_max(3+k_pow, pseudo_temp, deviational_temp, v_g)
            ) / v_g; 
            #else
            return util::random.random_dir() * util::random.accept_reject(
                util::lindisp::canonical_dist_factory(3+k_pow, pseudo_temp, deviational_temp, v_g), 
                0, omega_debye, util::lindisp::canonical_dist_max(3+k_pow, pseudo_temp, deviational_temp, v_g)
            ) / v_g; 
            #endif
        } 
        
        vec3 draw_k_linearized(PhononData& phonon, double equ_temp) const override {
            #ifndef BOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING
            return (is_umklapp?util::random.random_dir():arma::normalise(phonon.k)) * util::random.accept_reject(
                [equ_temp,this](double omega) {return util::lindisp::linearized_post_scatter_dist(2+k_pow, omega, equ_temp);}, 
                0, omega_debye,
                util::lindisp::linearized_post_scatter_dist_max(2+k_pow, equ_temp)
            ) / v_g; 
            #else
            return util::random.random_dir() * util::random.accept_reject(
                [equ_temp,this](double omega) {return util::lindisp::linearized_post_scatter_dist(2+k_pow, omega, equ_temp);}, 
                0, omega_debye,
                util::lindisp::linearized_post_scatter_dist_max(2+k_pow, equ_temp)
            ) / v_g; 
            #endif
        } 
    };
    
    class SimpleParametricExpTwoPhononProcess : public TwoPhononProcess {
        public:
        const unsigned int temp_pow;
        const unsigned int k_pow;
        const double fit_constant;
        const double v_g;
        const double omega_debye;
        const double expfac;

        SimpleParametricExpTwoPhononProcess() = delete;
        SimpleParametricExpTwoPhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye, double expfac) 
            : temp_pow(temp_pow), k_pow(k_pow), fit_constant(fit_constant), v_g(v_g), omega_debye(omega_debye), expfac(expfac) {}
        SimpleParametricExpTwoPhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye) 
            : SimpleParametricExpTwoPhononProcess(temp_pow, k_pow, fit_constant, v_g, omega_debye, util::HBAR * omega_debye / util::BMC){}
       
        double scatter_rate(const PhononData& phonon, double temp) const override {
            return fit_constant * util::ipow(arma::norm(phonon.k),k_pow) * util::ipow(temp,temp_pow) * std::exp(-expfac/temp);
        }
    };
    
    class SimpleParametricExpThreePhononProcess : public ThreePhononProcess {
        public:
        const unsigned int temp_pow;
        const unsigned int k_pow;
        const double fit_constant;
        const double v_g;
        const double omega_debye;
        const double expfac;
        const bool is_umklapp;

        SimpleParametricExpThreePhononProcess() = delete;
        SimpleParametricExpThreePhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye, bool is_umklapp, double expfac) 
            : temp_pow(temp_pow), k_pow(k_pow), fit_constant(fit_constant), v_g(v_g), omega_debye(omega_debye), expfac(expfac), is_umklapp(is_umklapp) {}
        SimpleParametricExpThreePhononProcess(unsigned temp_pow, unsigned k_pow, double fit_constant, double v_g, double omega_debye, bool is_umklapp) 
            : SimpleParametricExpThreePhononProcess(temp_pow, k_pow, fit_constant, v_g, omega_debye, is_umklapp, util::HBAR*omega_debye/util::BMC) {}
       
        double scatter_rate(const PhononData& phonon, double temp) const override {
            return fit_constant * util::ipow(arma::norm(phonon.k),k_pow) * util::ipow(temp,temp_pow) * std::exp(-expfac/temp);
        }

        double pseudo_temp(double scatter_sum, double temp) const override {
            scatter_sum = std::max(scatter_sum, 0.);
            scatter_sum *= util::ipow(v_g, k_pow);
            scatter_sum /= fit_constant * util::ipow(temp, temp_pow) * std::exp(-expfac/temp);
            double temp_guess = std::pow(
                    scatter_sum/(util::factorial(3+k_pow)*util::lindisp::riemann_zeta_tbl[4+k_pow])\
                    *2*util::ipow<2>(util::PI)*util::ipow<3>(v_g)*util::ipow(util::HBAR/util::BMC,4+k_pow),1./(4+k_pow));
            return util::lindisp::canonical_dist_pseudo_temp(
                k_pow+3,
                scatter_sum, v_g, temp_guess, omega_debye
            );
        }
        
        double scattered_omegasum(double scatter_temp, double temp) const override {
            return util::lindisp::canonical_dist_integral(3+k_pow, v_g, temp, omega_debye) / util::ipow(v_g, k_pow) \
                * fit_constant * util::ipow(scatter_temp,temp_pow)\
                * std::exp(-expfac/scatter_temp);
        } 
        
        vec3 draw_k(PhononData& phonon, double scatter_temp, double pseudo_temp, double deviational_temp) const override {
            // scatter_temp not needed here since it is eliminated by the normalisation of the distribution
            #ifndef BOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING
            return (is_umklapp?util::random.random_dir():arma::normalise(phonon.k)) * util::random.accept_reject(
                util::lindisp::canonical_dist_factory(3+k_pow, pseudo_temp, deviational_temp, v_g), 
                0, omega_debye, util::lindisp::canonical_dist_max(3+k_pow, pseudo_temp, deviational_temp, v_g)
            ) / v_g; 
            #else
            return util::random.random_dir() * util::random.accept_reject(
                util::lindisp::canonical_dist_factory(3+k_pow, pseudo_temp, deviational_temp, v_g), 
                0, omega_debye, util::lindisp::canonical_dist_max(3+k_pow, pseudo_temp, deviational_temp, v_g)
            ) / v_g; 
            #endif
        } 
        
        vec3 draw_k_linearized(PhononData& phonon, double equ_temp) const override {
            #ifndef BOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING
            return (is_umklapp?util::random.random_dir():arma::normalise(phonon.k)) * util::random.accept_reject(
                [equ_temp,this](double omega) {return util::lindisp::linearized_post_scatter_dist(2+k_pow, omega, equ_temp);}, 
                0, omega_debye,
                util::lindisp::linearized_post_scatter_dist_max(2+k_pow, equ_temp)
            ) / v_g; 
            #else
            return util::random.random_dir() * util::random.accept_reject(
                [equ_temp,this](double omega) {return util::lindisp::linearized_post_scatter_dist(2+k_pow, omega, equ_temp);}, 
                0, omega_debye,
                util::lindisp::linearized_post_scatter_dist_max(2+k_pow, equ_temp)
            ) / v_g; 
            #endif
        } 
    };
}}
