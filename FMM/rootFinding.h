//#define BOOST_MATH_INSTRUMENT

#include "grid.h"

#include <boost/math/tools/roots.hpp>
//using boost::math::policies::policy;
//using boost::math::tools::newton_raphson_iterate;
//using boost::math::tools::halley_iterate; //
using boost::math::tools::eps_tolerance; // Binary functor for specified number of bits.
using boost::math::tools::bracket_and_solve_root;
//using boost::math::tools::toms748_solve;

#include <boost/math/special_functions/next.hpp> // For float_distance.

#include <iostream>
//using std::cout; using std::endl;
#include <iomanip>
//using std::setw; using std::setprecision;
#include <limits>
//using std::numeric_limits;
#include <algorithm>
#include "math.h"

double solver(double guess, grid& g, int i, int j);