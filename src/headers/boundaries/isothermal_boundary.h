#pragma once

#include "../boundary.h"
#include <memory>

namespace boltzi {
    class IsothermalBoundary : public Boundary {
        public:
        IsothermalBoundary() = delete;
        IsothermalBoundary(const std::vector<Primitive>& primitives) : Boundary(primitives) {}
        
        double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const final  {
            throw BoltziException("isothermal boundaries do not reflect.");
        } 

        vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const final {
            throw BoltziException("isothermal boundaries do not reflect.");
        }
        bool is_active() const final {return true;}
        bool is_generator() const final {return true;}
        double p_absorb(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const final {return 1;};
    };
}
