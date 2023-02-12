#pragma once

#include <vector>

#include "../boltzi_exception.h"
#include "isothermal_boundary.h"

namespace boltzi {
namespace boundaries {
    class SimpleIsothermalBoundary : public IsothermalBoundary {
        public:
        const double temperature;
        SimpleIsothermalBoundary(const std::vector<Primitive>& primitives, double temperature) : IsothermalBoundary(primitives), temperature(temperature) {}
        virtual double surface_temperature(vec3 pos) const override {
            return temperature;
        }
    };
}}
