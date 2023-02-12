#pragma once

#include <cmath>
#include <vector>

#include "../boltzi_exception.h"
#include "../util/basics.h"
#include "vacuum_boundary.h"

namespace boltzi {
namespace boundaries {
    class EvenVacuumBoundary : public VacuumBoundary {
        public:
        EvenVacuumBoundary(const std::vector<Primitive>& primitives) : VacuumBoundary(primitives) {}

        double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override  {
            return 0;        
        } 

        vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override {
            double theta = std::acos(1-util::random.draw_real()); 
            vec3 dir = util::sphere_to_cart(vec2{util::random.draw_real()*2*util::PI,theta});
            return  (primitives[coll->primitive].normal * dir[2] + primitives[coll->primitive].e_1 * dir[0] + primitives[coll->primitive].e_2 * dir[1]) 
                    * arma::norm(phonon.k);
        }
    };
}}
