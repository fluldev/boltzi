#include "boltzi_config.h"

namespace boltzi {
namespace config {
    const unsigned int BINDIST_BINS = 1000; 
    const unsigned int CELL_OVERLAP_SAMPLES = 100; 
    const double WELD_RADIUS = .01;
    const double EULER_DELTA = 1e-6;
    const double EULER_STEP = .2;
    const double INSIDE_CHECK_VOLUME_MARGIN = 1.005;
    const double INSIDE_CHECK_AREA_MARGIN = 1.005;
    const unsigned CPUS = 8;  // amount of cpus - 1 seems to be the quickest
    const double SAVE_MOVE_STEP = 1e-3;
    const unsigned MINIMUM_TEMPERATURE_PHONONS = 10;
    const double MIN_ARTIFACT_VOLUME=.03;
    const unsigned INFORMER_SLEEP_DURATION=100;
    const unsigned MAX_ACCEPT_REJECT_ATTEMPTS=50000;
    const double DOUBLE_COLL_IGNORE_RADIUS = 1e-9;
}}
