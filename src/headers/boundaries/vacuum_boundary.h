#pragma once

#include "../boundary.h"
#include "../boltzi_exception.h"

namespace boltzi {
    class VacuumBoundary : public Boundary {
        public:
        VacuumBoundary() = delete;
        VacuumBoundary(const std::vector<Primitive>& primitives) : Boundary(primitives) {}
        bool is_active() const final {return true;}
        bool is_generator() const final {return false;}
        double p_absorb(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const final {return 0;};
        virtual double surface_temperature(vec3 pos) const{
           throw BoltziException("VacuumBoundary does not have a temperature."); 
        };
    };
}
