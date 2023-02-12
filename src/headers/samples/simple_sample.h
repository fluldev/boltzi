#pragma once

#include "../sample.h"

namespace boltzi {
namespace samples {
    class SimpleSample : public Sample {
        public:
        SimpleSample() = delete;
        SimpleSample(
            double n_eff,
            std::tuple<unsigned, unsigned, unsigned> n_bins,
            const SampleGeometry& geometry,
            const Material& material,
            std::function<double(vec3 r)> deviational_temp = [] (vec3 r) -> double {return 0;},
            std::function<vec3(vec3 r)> deviational_temp_gradient = [] (vec3 r) -> vec3 {return vec3{0,0,0};}
        ) : Sample(n_eff, n_bins, geometry, std::vector<Material>{material}, deviational_temp, deviational_temp_gradient) {init_cells();}

        const Material& associated_material(const PhononData& phonon) const override {return materials.front();}

        const Material& associated_material(vec3 r) const override {return materials.front();}

        unsigned material_index(vec3 r) const override {
            return 0;
        }
    };
}}
