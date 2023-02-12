#pragma once

#include <vector>

#include "../boltzi_exception.h"
#include "vacuum_boundary.h"

namespace boltzi {
namespace boundaries {
    class SpecularVacuumBoundary : public VacuumBoundary {
        public:
        SpecularVacuumBoundary(const std::vector<Primitive>& primitives) : VacuumBoundary(primitives) {}

        double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override  {
            return 1;
        } 

        vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const override {
            throw BoltziException("diffusive_draw not supported for SpecularVacuumBoundary.");
        }
    };
}}
