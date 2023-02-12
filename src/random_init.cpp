#include "util/random.h"

// this file must be included in every build of boltzi but NOT when building the library

namespace boltzi {
namespace util {
    thread_local RandomUtil random;
}}
