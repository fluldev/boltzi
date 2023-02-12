#pragma once

#include <string>

#include "boundary.h"

namespace boltzi{
SampleGeometry load_geometry(std::string fname, double scale=1.0);
}
