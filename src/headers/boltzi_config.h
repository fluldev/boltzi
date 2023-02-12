#pragma once

namespace boltzi {
namespace config {
    extern const unsigned int BINDIST_BINS; 
    extern const unsigned int CELL_OVERLAP_SAMPLES; 
    extern const double WELD_RADIUS;
    extern const double EULER_DELTA;
    extern const double EULER_STEP;
    extern const double INSIDE_CHECK_VOLUME_MARGIN;
    extern const double INSIDE_CHECK_AREA_MARGIN;
    extern const unsigned CPUS;
    extern const double SAVE_MOVE_STEP;
    extern const unsigned MINIMUM_TEMPERATURE_PHONONS;
    extern const double MIN_ARTIFACT_VOLUME;
    extern const unsigned INFORMER_SLEEP_DURATION;
    extern const unsigned MAX_ACCEPT_REJECT_ATTEMPTS;
    extern const double DOUBLE_COLL_IGNORE_RADIUS;
}}
