#pragma once

#include "../boundary.h"
#include "../boltzi_exception.h"

namespace boltzi {
namespace boundaries {
    class MeasureBoundary : public Boundary {
        public:
        MeasureBoundary() = delete;
        MeasureBoundary(const std::vector<Primitive>& primitives) : Boundary(primitives) {}
        bool is_active() const final {return false;}
        bool is_generator() const final {return false;}
        double p_absorb(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const final {return 0;};
        double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override  {return 0;}
        virtual double surface_temperature(vec3 pos) const override{
           throw BoltziException("MeasureBoundary does not have a temperature."); 
        }
        vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override {
           throw BoltziException("MeasureBoundary does not allow collisions."); 
        }

    };
}
}
