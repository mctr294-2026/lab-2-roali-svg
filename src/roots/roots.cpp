#include "roots.hpp"

#include <cmath>

// these are the constants required by the lab 
const double TOL = 1e-6;  // answer good enough if f(x)<= tol
const int MAX_ITERS = 1000000; // cap 
const double EPS = 1e-14;

// bisection is js if f(a) and f(b) have opp signs its a root 
// somewhere between a and b it halves the interval , the half that contains the sign changes. 
// othereise i retuen false 
bool bisection(std::function<double(double)> f,
               double a, double b,
               double* root)
{
    // a root is invalid we cant write the ans
    if (!root) return false;
    if (a > b) std::swap(a, b);

    // evalute endpoints once 
    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= TOL) { *root = a; return true; }
    if (std::abs(fb) <= TOL) { *root = b; return true; }
    // Bisection requires a sign change across [a,b]
    // If fa and fb have the same sign, no guarantee of a root ie it fails 
    if (fa * fb > 0.0) return false;

    for (int i = 0; i < MAX_ITERS; i++) {
        // Midpoint splits the interval in half
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (std::abs(fc) <= TOL || std::abs(b - a) <= TOL) {
            *root = c;
            return true;
        }

        if (fa * fc < 0.0) { b = c; fb = fc; }
        else               { a = c; fa = fc; }
    }
    return false;
}


// for regula falsi its kinda the samething as bisection but instead of using the midppoint
 //we use the x-inter of sec line between 
 
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= TOL) { *root = a; return true; }
    if (std::abs(fb) <= TOL) { *root = b; return true; }
    // getting intial bracket 
    if (fa * fb > 0.0) return false;

    for (int i = 0; i < MAX_ITERS; i++) {
        //if deno is close to zero sec line is alomst flat 
        double denom = fb - fa;
        if (std::abs(denom) <= EPS) return false;

        double c = (a * fb - b * fa) / denom;
        double fc = f(c);

        if (std::abs(fc) <= TOL || std::abs(b - a) <= TOL) {
            *root = c;
            return true;
        }

        if (fa * fc < 0.0) { b = c; fb = fc; }
        else               { a = c; fa = fc; }
    }
    return false;
}

    // uses tangent line at x to jump closer to the root 
    // also converges fast if the starting guess is good 
    // x_next = x - f(x) / f'(x)
    // If derivative is close to 0 or the guess leaves [a,b], I return false.
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);
    //c is our intial guess fails if its not inside interval 
    if (c < a || c > b) return false;

    double x = c;

    for (int i = 0; i < MAX_ITERS; i++) {
        double fx = f(x);
        // if its closer to zero we are done 
        if (std::abs(fx) <= TOL) { *root = x; return true; }

        double gx = g(x);
        if (std::abs(gx) <= EPS) return false;

        double x_next = x - fx / gx;

        
        if (!std::isfinite(x_next) || x_next < a || x_next > b) return false;

        if (std::abs(x_next - x) <= TOL) { *root = x_next; return true; }
        x = x_next;
    }
    return false;
}

    //secant 
// Like Newton, but we approximate derivative using two points.
// x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
// Also: fail if iteration leaves [a,b].
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);
    // c is initial guess must be inside [a,b]
    if (c < a || c > b) return false;

    double fa = f(a);
    double fb = f(b); // Evaluate function at the three points
    double fc = f(c);

    if (std::abs(fa) <= TOL) { *root = a; return true; }
    if (std::abs(fb) <= TOL) { *root = b; return true; }
    if (std::abs(fc) <= TOL) { *root = c; return true; }

    
        // If c forms a sign-change bracket with an endpoint, use that pair.
    // Otherwise fallback to (a,b).
    double x0, x1, f0, f1;
    if (fa * fc < 0.0)    // root likely between a and c
      { x0 = a; f0 = fa; x1 = c; f1 = fc; }
    else if (fb * fc < 0.0) 
    { x0 = b; f0 = fb; x1 = c; f1 = fc; } // root likely between b and c
    else                    
    { x0 = a; f0 = fa; x1 = b; f1 = fb; } // no bracket info, just use endpoints

    for (int i = 0; i < MAX_ITERS; i++) {
        // Denominator is the slope estimate (f1 - f0)
        double denom = f1 - f0;
        if (std::abs(denom) <= EPS) return false;

        double x2 = x1 - f1 * (x1 - x0) / denom;

        if (!std::isfinite(x2) || x2 < a || x2 > b) return false;

        double f2 = f(x2);

        if (std::abs(f2) <= TOL || std::abs(x2 - x1) <= TOL) {
            *root = x2;
            return true;
        }

        x0 = x1; f0 = f1;
        x1 = x2; f1 = f2;
    }
    return false;
}


