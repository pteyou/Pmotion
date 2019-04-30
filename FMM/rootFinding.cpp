#include "rootFinding.h"

struct my_functor_noderiv
{ 
  //  cube root of x using only function - no derivatives.
  my_functor_noderiv(grid &g, int i, int j) : g_(g), i_(i), j_(j) { }
  double operator()(double T)
  {
    double res {0};
    res += pow(std::max( (T - g_.get(i_-1, j_).T_) / g_.dx_, 0.0), 2.0);
    res += pow(std::min( (g_.get(i_+1, j_).T_ - T) / g_.dx_, 0.0), 2.0);
    res += pow(std::max( (T - g_.get(i_, j_-1).T_) / g_.dy_, 0.0), 2.0);
    res += pow(std::min( (g_.get(i_, j_+1).T_ - T) / g_.dy_, 0.0), 2.0);
    res -= 1.0;
    return res;
  }
private:
  grid &g_;
  int i_, j_;
};

struct my_functor_noderiv2
{ 
  //  cube root of x using only function - no derivatives.
  my_functor_noderiv2(grid &g, int i, int j) : g_(g), i_(i), j_(j) { }
  double operator()(double T)
  {
    double res {0};
    res += pow(std::max( std::max((T - g_.get(i_-1, j_).T_) / g_.dx_, 0.0), 
                        -std::min((g_.get(i_+1, j_).T_ - T) / g_.dx_, 0.0)), 2.0);

    res += pow(std::max( std::max((T - g_.get(i_, j_-1).T_) / g_.dy_, 0.0), 
                        -std::min((g_.get(i_, j_+1).T_ - T) / g_.dy_, 0.0)), 2.0);
    res -= 1.0;
    return res;
  }
private:
  grid &g_;
  int i_, j_;
};


double solver(double guess, grid& g, int i, int j)
{ 
  double factor = 2;                                 // How big steps to take when searching.

  const boost::uintmax_t maxit = 100;            // Limit to maximum iterations.
  boost::uintmax_t it = maxit;                  // Initally our chosen max iterations, but updated with actual.
  bool is_rising = true;                        // So if result if guess^3 is too low, then try increasing guess.
  int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  // Some fraction of digits is used to control how accurate to try to make the result.
  int get_digits = digits - 5;                  // We have to have a non-zero interval at each step, so
                                                // maximum accuracy is digits - 1.  But we also have to
                                                // allow for inaccuracy in f(x), otherwise the last few
                                                // iterations just thrash around.
  eps_tolerance<double> tol(get_digits);             // Set the tolerance.
  std::pair<double, double> r = bracket_and_solve_root(my_functor_noderiv(g, i, j), guess, factor, is_rising, tol, it);
  if (it == maxit) std::cout << " solver diverged at indices " << i << "\t" << j << std::endl;
  return r.first + (r.second - r.first)/2;      // Midway between brackets is our result, if necessary we could
                                                // return the result as an interval here.
}


