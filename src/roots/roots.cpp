#include "roots.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>

// Lab constants
static constexpr double TOL = 1e-6;
static constexpr std::size_t MAX_ITERS = 1'000'000;
static constexpr double EPS = 1e-14;  // protect divide-by-zero-ish cases

bool bisection(std::function<double(double)> f,
               double a, double b,
               double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= TOL) { *root = a; return true; }
    if (std::abs(fb) <= TOL) { *root = b; return true; }

    // Guaranteed only if signs differ
    if (fa * fb > 0.0) return false;

    for (std::size_t i = 0; i < MAX_ITERS; ++i) {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (std::abs(fc) <= TOL || std::abs(b - a) <= TOL) {
            *root = c;
            return true;
        }

        // Keep the sign-change bracket
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return false;
}

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

    if (fa * fb > 0.0) return false;

    for (std::size_t i = 0; i < MAX_ITERS; ++i) {
        double denom = (fb - fa);
        if (std::abs(denom) <= EPS) return false;

        // x-intercept of the line through (a,fa) and (b,fb)
        double c = (a * fb - b * fa) / denom;
        double fc = f(c);

        if (std::abs(fc) <= TOL || std::abs(b - a) <= TOL) {
            *root = c;
            return true;
        }

        // Keep the sign-change bracket
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);
    if (c < a || c > b) return false;

    double x = c;

    for (std::size_t i = 0; i < MAX_ITERS; ++i) {
        double fx = f(x);
        if (std::abs(fx) <= TOL) { *root = x; return true; }

        double gx = g(x);
        if (std::abs(gx) <= EPS) return false;

        double x_next = x - fx / gx;

        // Spec: fail if it leaves the interval
        if (!std::isfinite(x_next) || x_next < a || x_next > b) return false;

        if (std::abs(x_next - x) <= TOL) { *root = x_next; return true; }
        x = x_next;
    }
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double* root)
{
    if (!root) return false;
    if (a > b) std::swap(a, b);
    if (c < a || c > b) return false;

    double fa = f(a);
    double fb = f(b);
    double fc = f(c);

    if (std::abs(fa) <= TOL) { *root = a; return true; }
    if (std::abs(fb) <= TOL) { *root = b; return true; }
    if (std::abs(fc) <= TOL) { *root = c; return true; }

    // Choose the better initial pair using c:
    double x0, x1, f0, f1;

    if (fa * fc < 0.0) {          // bracket between a and c
        x0 = a; f0 = fa;
        x1 = c; f1 = fc;
    } else if (fb * fc < 0.0) {   // bracket between b and c
        x0 = b; f0 = fb;
        x1 = c; f1 = fc;
    } else {                      // fallback
        x0 = a; f0 = fa;
        x1 = b; f1 = fb;
    }

    for (std::size_t i = 0; i < MAX_ITERS; ++i) {
        double denom = (f1 - f0);
        if (std::abs(denom) <= EPS) return false;

        double x2 = x1 - f1 * (x1 - x0) / denom;

        // Spec says fail if it leaves the interval
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

