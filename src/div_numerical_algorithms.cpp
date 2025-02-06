#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Numerical Gradient Approximation Using Richardson Extrapolation
//'
//' @title Approximate Numerical Gradient for a Function \code{f}
//' @name numerical_gradient
//'
//' @param x A \code{double}. The point where \code{f} is to be evaluated.
//' @param f A \code{function}. If the function have more than 1 argument, ensure that these are fixed.
//' @param eps A \code{double}. The delta value in which the gradient are approximated.
//' @param r An \code{integer}. The number of Iterations utilized for the Richardson extrapolation algorithm.
//'
//' @details
//' Given a particular vale of \code{x}, approximates the gradient of \code{f(x)} at this point.
//'
//' @return A \code{double} that signify the approximated gradient of \code{f} evaluated at \code{x}.
//' @export
//'
//' @examples
//' example_fun <- function(x, a = 3){return((x - a)^2)}
//' gradient_at_minimum <- numerical_gradient(3, function(x) example_fun(3, 3))
//'
//' # Expected to be approximately '0'
//' print(gradient_at_minimum)

// [[Rcpp::export]]
double numerical_gradient(double x, const Function f, double eps = 1e-4, int r = 4) {
  std::vector<double> D(r, 0.0);
  double h = eps;

  for(int i = 0; i < r; ++i) {
    double step = h / pow(2, i);
    double f_plus = as<double>(f(x + step));
    double f_minus = as<double>(f(x - step));
    D[i] = (f_plus - f_minus) / (2 * step);

    // Richardson acceleration
    for(int j = i; j > 0; --j)
      D[j-1] = (4.0 * D[j] - D[j-1]) / 3.0;
  }
  return D[0];
}

//' Numerical Minimization Using Brent's Algorithm
//'
//' @title Obtain the value of \code{x} so that \code{f(x)} is a local minimum
//' @name brent_min
//'
//' @param f A \code{function}. Must be defined and continuous in the domain \code{[a, b]}. If the function have more than 1 argument, ensure that these are fixed.
//' @param a A \code{double}. The lower bound for the search interval.
//' @param b A \code{double}. The upper bound for the search interval.
//' @param tol A \code{double}. The desired accuracy.
//' @param max_iter An \code{integer}. The maximum number of iteration that should be performed before the algorithm forces to quit.
//'
//' @details
//' Given a particular function \code{f}, obtains a local minimum of \code{f} in the domain \code{[a, b]}. Brent's algorithm uses a combination of golden section search and quadratic interpolation to obtain the minimum.
//'
//' @return A \code{double} that signify the value of \code{x} that results in local minimum of \code{f}.
//' @export
//'
//' @examples
//' example_fun <- function(x, a = 3){return((x - a)^2)}
//' the_minimum <- brent_min(f = function(x) example_fun(x, 3), a = -3, b = 6)
//'
//' # Expected to be approximately '3'
//' print(the_minimum)
//'
//' # The value of f at the minimum
//' print(example_fun(the_minimum))

// [[Rcpp::export]]
double brent_min(Function f, double a, double b, double tol = 1.22e-4, int max_iter = 100) {
  const double golden = (3 - sqrt(5)) / 2;
  double x = a + golden * (b - a);
  double w = x, v = x, fx = as<double>(f(x)), fw = fx, fv = fx, d = 0, e = 0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double m = 0.5 * (a + b);
    double tol_act = tol * abs(x) + tol;
    if (abs(x - m) <= 2 * tol_act - 0.5 * (b - a)) return x;

    bool use_para = (abs(e) > tol_act);
    double p = 0, q = 0, r_val = 0;
    if (use_para) {
      r_val = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r_val;
      q = 2 * (q - r_val);
      if (q > 0) p = -p;
      else q = -q;
      use_para = (abs(p) < abs(0.5 * q * e) && p > q * (a - x) && p < q * (b - x));
      if (use_para) {
        e = d;
        d = p / q;
        double u = x + d;
        if (u - a < 2 * tol_act || b - u < 2 * tol_act)
          d = (x < m) ? tol_act : -tol_act;
      }
    }
    if (!use_para) {
      e = (x < m) ? b - x : a - x;
      d = golden * e;
    }

    double u = x + ((abs(d) >= tol_act) ? d : (d > 0 ? tol_act : -tol_act));
    double fu = as<double>(f(u));

    if (fu <= fx) {
      (u >= x) ? (a = x) : (b = x);
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      (u < x) ? (a = u) : (b = u);
      if (fu <= fw || w == x) {
        v = w; w = u;
        fv = fw; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

//' Numerical Initial Minimization Using a Stepwise Approach
//'
//' @title Obtain the Minimum Value of \code{x} Where \code{f(x)} Reaches a Minimum
//' @name first_min
//'
//' @param f A \code{function}. Must be defined and continuous in the domain \code{[a, b]}. If the function have more than 1 argument, ensure that these are fixed.
//' @param a A \code{double}. The lower bound for the search interval.
//' @param b A \code{double}. The upper bound for the search interval.
//' @param tol A \code{double}. The desired accuracy.
//' @param step_init A \code{double}. The initial step size multiplier. This number of multiplied by the length of \code{[a, b]} to get the initial step size \eqn{\Delta_0}.
//' @param max_iter An \code{integer}. The maximum number of iteration that should be performed before the algorithm forces to quit.
//'
//' @details
//' Given a particular function \code{f}, obtains the first local minimum of \code{f} in the domain \code{[a, b]}.
//' For each iteration, we check whether \eqn{f'(\cdot) > 0}.
//' The i-th search interval is given by \eqn{[a + \Delta_0 \sum_{j=0}^{i-2} 2^{j}, a + \Delta_0 \sum_{j=0}^{i-1} 2^{j}\,]}.
//' If a minimum is not found before \eqn{a + \Delta_0 \sum_{j=0}^{i-1} 2^{j} > b}, \code{brent_min} is called instead.
//' Note that this algorithm works best if \code{f} is a differentiable and smooth function. Note also that If \eqn{f'(a) > 0}, this algorithm will
//' return \code{a}, so this algorithm is most useful if \eqn{f'(x) < 0} in the start of the interval \code{[a, b]}.
//'
//' @return A \code{double} that signify the value of \code{x} that results in the initial minimum of \code{f}.
//' @export
//'
//' @examples
//' example_fun <- function(x, a = 3){return((x - a)^2)}
//' the_first_minimum <- first_min(f = function(x) example_fun(x, 3), a = -3, b = 6)
//'
//' # Expected to be approximately '3'
//' print(the_first_minimum)
//'
//' # The value of f at the minimum
//' print(example_fun(the_first_minimum))


// [[Rcpp::export]]
double first_min(Function f, double a, double b, double tol = 1e-8, double step_init = 1e-4, int max_iter = 1000) {
  double x = a, step = step_init * (b - a);
  double g = numerical_gradient(x, f);

  // Edge case: initial point is minimum
  if(g > 0) return x;

  // Bracketing phase
  while(x < b && step > 1e-15) {
    double x_next = std::min(x + step, b);
    double g_next = numerical_gradient(x_next, f);

    if(g_next > 0) {
      return brent_min(f, x, x_next, tol, max_iter);
    }

    step *= 2;

    // Prevent infinite loops
    if(x + step > b && x_next == b) {
      break;
    }
    x = x_next;
    g = g_next;
  }

  // Fallback to Brent's on full interval
  return brent_min(f, a, b, tol, max_iter);
}

