#pragma once

#include <vector>

#include "../boltzi_exception.h"
#include "../util/basics.h"
#include "vacuum_boundary.h"

namespace boltzi {
namespace boundaries {
    class BKVacuumBoundary : public VacuumBoundary {
        const double eta;
        public:
        BKVacuumBoundary(const std::vector<Primitive>& primitives, double std_dev) : VacuumBoundary(primitives), eta(std_dev) {}

        double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override  {
            return std::exp(
                -util::ipow<2>(2*eta*arma::norm(phonon.k)
                *arma::dot(primitives[coll->primitive].normal,arma::normalise(phonon.k)))
            );
        } 

        vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override {
            double theta = std::asin(util::random.draw_real()); 
            vec3 dir = util::sphere_to_cart(vec2{util::random.draw_real()*2*util::PI,theta});
            return  (primitives[coll->primitive].normal * dir[2] + primitives[coll->primitive].e_1 * dir[0] + primitives[coll->primitive].e_2 * dir[1])
                    * arma::norm(phonon.k);
        }
    };
}}
