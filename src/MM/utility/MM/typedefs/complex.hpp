#pragma once
#include "header.hpp"

struct Complex {
    // Constructors
    Complex(double r = 0.0, double i = 0.0) : re(r), im(i) {}
    Complex(const Complex&) = default;
    Complex(Complex&&) = default;
    Complex& operator=(const Complex&) = default;
    Complex& operator=(Complex&&) = default;

    // Operators
    // Arithmetic with another complex
    inline Complex operator+(const Complex& o) const { return {re + o.re, im + o.im}; }
    inline Complex operator-(const Complex& o) const { return {re - o.re, im - o.im}; }
    inline Complex operator*(const Complex& o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
    inline Complex operator/(const Complex& o) const {
        double denom = o.re * o.re + o.im * o.im;
        if (denom == 0.0 || !std::isfinite(denom) || std::isnan(denom))
            throw std::runtime_error("Complex division by zero or invalid denominator");
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im))
            throw std::runtime_error("Complex division: numerator not finite");
        return {(re * o.re + im * o.im) / denom, (im * o.re - re * o.im) / denom};
    }
    // Arithmetic with real
    inline Complex operator+(double r) const { return {re + r, im}; }
    inline Complex operator-(double r) const { return {re - r, im}; }
    inline Complex operator*(double r) const { return {re * r, im * r}; }
    inline Complex operator/(double r) const {
        if (r == 0.0 || !std::isfinite(r) || std::isnan(r))
            throw std::runtime_error("Complex division by zero or invalid real");
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im))
            throw std::runtime_error("Complex division: numerator not finite");
        return {re / r, im / r};
    }
    // Arithmetic with integer types
    inline Complex operator+(int r) const { return {re + r, im}; }
    inline Complex operator-(int r) const { return {re - r, im}; }
    inline Complex operator*(int r) const { return {re * r, im * r}; }
    inline Complex operator/(int r) const {
        if (r == 0)
            throw std::runtime_error("Complex division by zero (int)");
        return {re / r, im / r};
    }
    friend Complex operator*(double r, const Complex& c) { return c * r; }
    friend Complex operator/(double r, const Complex& c) {
        double denom = c.re * c.re + c.im * c.im;
        if (denom == 0.0 || !std::isfinite(denom) || std::isnan(denom))
            throw std::runtime_error("Complex division by zero or invalid denominator");
        if (!std::isfinite(r) || std::isnan(r))
            throw std::runtime_error("Complex division: numerator not finite");
        return Complex(r) / c;
    }
    friend Complex operator*(int r, const Complex& c) { return c * r; }
    friend Complex operator/(int r, const Complex& c) {
        double denom = c.re * c.re + c.im * c.im;
        if (denom == 0.0 || !std::isfinite(denom) || std::isnan(denom))
            throw std::runtime_error("Complex division by zero or invalid denominator");
        return Complex(static_cast<double>(r)) / c;
    }
    // Compound assignment
    inline Complex& operator+=(const Complex& o) { re += o.re; im += o.im; return *this; }
    inline Complex& operator-=(const Complex& o) { re -= o.re; im -= o.im; return *this; }
    inline Complex& operator*=(const Complex& o) { *this = *this * o; return *this; }
    inline Complex& operator/=(const Complex& o) {
        double denom = o.re * o.re + o.im * o.im;
        if (denom == 0.0 || !std::isfinite(denom) || std::isnan(denom))
            throw std::runtime_error("Complex division by zero or invalid denominator");
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im))
            throw std::runtime_error("Complex division: numerator not finite");
        *this = *this / o;
        return *this;
    }
    // Compound assignment with real
    inline Complex& operator+=(double r) { re += r; return *this; }
    inline Complex& operator-=(double r) { re -= r; return *this; }
    inline Complex& operator*=(double r) { re *= r; im *= r; return *this; }
    inline Complex& operator/=(double r) {
        if (r == 0.0 || !std::isfinite(r) || std::isnan(r))
            throw std::runtime_error("Complex division by zero (real)");
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im))
            throw std::runtime_error("Complex division: numerator not finite");
        re /= r; im /= r; return *this;
    }
    // Comparison
    inline bool operator==(const Complex& o) const {
        constexpr double eps = 1e-12;
        if (!std::isfinite(re) || !std::isfinite(im) || !std::isfinite(o.re) || !std::isfinite(o.im)
            || std::isnan(re) || std::isnan(im) || std::isnan(o.re) || std::isnan(o.im))
            return false;
        return std::abs(re - o.re) < eps && std::abs(im - o.im) < eps;
    }
    inline bool operator!=(const Complex& o) const { return !(*this == o); }
    // Comparison with real numbers
    inline bool operator==(double r) const {
        constexpr double eps = 1e-12;
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im) || !std::isfinite(r) || std::isnan(r))
            return false;
        return std::abs(re - r) < eps && std::abs(im) < eps;
    }
    inline bool operator!=(double r) const { return !(*this == r); }
    // Unary operators
    inline Complex operator+() const { return *this; }
    inline Complex operator-() const { return {-re, -im}; }
    // operator[] for real/imag access
    double& operator[](int idx) {
        if (idx == 0) return re;
        else if (idx == 1) return im;
        else throw std::out_of_range("Complex index must be 0 (real) or 1 (imag)");
    }
    const double& operator[](int idx) const {
        if (idx == 0) return re;
        else if (idx == 1) return im;
        else throw std::out_of_range("Complex index must be 0 (real) or 1 (imag)");
    }
    // Explicit bool conversion
    explicit operator bool() const { return re != 0.0 || im != 0.0; }

    // Basic functions
    inline Complex exp() const {
        double e = std::exp(re);
        return {e * std::cos(im), e * std::sin(im)};
    }
    inline Complex log() const { return {std::log(this->abs()), this->arg()}; }
    inline Complex pow(const Complex& o) const {
        // z^w = exp(w * log(z))
        return (o * this->log()).exp();
    }
    inline Complex pow(double p) const {
        // z^p = exp(p * log(z))
        return (this->log() * p).exp();
    }
    inline Complex log(int branch) const { return {std::log(this->abs()), this->arg() + 2 * PI * branch}; }
    inline Complex pow(const Complex& o, int branch) const { return (this->log(branch) * o).exp(); }
    inline Complex pow(double p, int branch) const { return (this->log(branch) * p).exp(); }
    inline Complex sqrt() const {
        double r = this->abs();
        double theta = this->arg() / 2.0;
        return from_polar(std::sqrt(r), theta);
    }
    inline Complex inverse() const {
        double d = re * re + im * im;
        if (d == 0.0 || !std::isfinite(d) || std::isnan(d))
            throw std::runtime_error("Complex inverse division by zero or invalid denominator");
        if (!std::isfinite(re) || !std::isfinite(im) || std::isnan(re) || std::isnan(im))
            throw std::runtime_error("Complex inverse: numerator not finite");
        return {re / d, -im / d};
    }
    // n-th roots (returns principal root)
    inline Complex root(int n) const {
        if (n == 0) throw std::runtime_error("Complex root: n == 0");
        double r = std::pow(this->abs(), 1.0 / n);
        double theta = this->arg() / n;
        return from_polar(r, theta);
    }
    // Polar form output
    std::string to_polar() const {
        std::ostringstream oss;
        oss << abs() << " * exp(i*" << arg() << ")";
        return oss.str();
    }
    // Static abs/arg
    static double abs(const Complex& c) { return c.abs(); }
    static double arg(const Complex& c) { return c.arg(); }
    // Static polar form constructor
    static inline Complex from_polar(double r, double theta) { return {r * std::cos(theta), r * std::sin(theta)}; }
    // Template support for other floating-point types
    template<typename T>
    static Complex from(const T& c) {
        if constexpr (requires { c.real(); c.imag(); }) {
            return Complex(static_cast<double>(c.real()), static_cast<double>(c.imag()));
        } else if constexpr (requires { c.re; c.im; }) {
            return Complex(static_cast<double>(c.re), static_cast<double>(c.im));
        } else {
            static_assert(sizeof(T) == 0, "Type does not have real/imag or re/im members");
        }
    }

    // Trigonometric functions
    inline Complex sin() const { return {std::sin(re) * std::cosh(im), std::cos(re) * std::sinh(im)}; }
    inline Complex cos() const { return {std::cos(re) * std::cosh(im), -std::sin(re) * std::sinh(im)}; }
    inline Complex tan() const {
        Complex c = this->cos();
        double d = c.re * c.re + c.im * c.im;
        if (d == 0.0 || !std::isfinite(d) || std::isnan(d) || !std::isfinite(c.re) || std::isnan(c.re) || !std::isfinite(c.im) || std::isnan(c.im))
            throw std::runtime_error("Complex tan division by zero or invalid denominator");
        Complex s = this->sin();
        if (!std::isfinite(s.re) || std::isnan(s.re) || !std::isfinite(s.im) || std::isnan(s.im))
            throw std::runtime_error("Complex tan: numerator not finite");
        return s / c;
    }

    // Inverse trigonometric functions
    inline Complex asin() const {
        Complex iz(-im, re); // i*z
        Complex sqrt_term = (Complex(1.0, 0.0) - (*this) * (*this)).sqrt();
        Complex log_arg = iz + sqrt_term;
        Complex res = log_arg.log();
        return Complex(res.im, -res.re); // -i*log(...)
    }
    inline Complex acos() const {
        Complex sqrt_term = (Complex(1.0, 0.0) - (*this) * (*this)).sqrt();
        Complex log_arg = *this + sqrt_term * Complex(0.0, 1.0);
        Complex res = log_arg.log();
        return Complex(res.im, -res.re); // -i*log(...)
    }
    inline Complex atan() const {
        // atan(z) = (i/2) * [log(1 - iz) - log(1 + iz)]
        Complex iz(-im, re); // i*z
        Complex one(1.0, 0.0);
        Complex log1 = (one - iz).log();
        Complex log2 = (one + iz).log();
        return Complex(0.0, 0.5) * (log1 - log2);
    }

    // Reciprocal trigonometric functions
    inline Complex csc() const {
        Complex s = this->sin();
        double d = s.re * s.re + s.im * s.im;
        if (d == 0.0 || !std::isfinite(d))
            throw std::runtime_error("Complex csc division by zero or invalid denominator");
        return Complex(1.0) / s;
    }
    inline Complex sec() const {
        Complex c = this->cos();
        double d = c.re * c.re + c.im * c.im;
        if (d == 0.0 || !std::isfinite(d))
            throw std::runtime_error("Complex sec division by zero or invalid denominator");
        return Complex(1.0) / c;
    }
    inline Complex cot() const {
        Complex s = this->sin();
        double d = s.re * s.re + s.im * s.im;
        if (d == 0.0 || !std::isfinite(d) || std::isnan(d) || !std::isfinite(s.re) || std::isnan(s.re) || !std::isfinite(s.im) || std::isnan(s.im))
            throw std::runtime_error("Complex cot division by zero or invalid denominator");
        Complex c = this->cos();
        if (!std::isfinite(c.re) || std::isnan(c.re) || !std::isfinite(c.im) || std::isnan(c.im))
            throw std::runtime_error("Complex cot: numerator not finite");
        return c / s;
    }

    // Inverse reciprocals of trigonometric functions
    inline Complex acsc() const { return this->inverse().asin(); }
    inline Complex asec() const { return this->inverse().acos(); }
    inline Complex acot() const { return this->inverse().atan(); }

    // Hyperbolic functions
    inline Complex sinh() const { return {std::sinh(re) * std::cos(im), std::cosh(re) * std::sin(im)}; }
    inline Complex cosh() const { return {std::cosh(re) * std::cos(im), std::sinh(re) * std::sin(im)}; }
    inline Complex tanh() const { return this->sinh() / this->cosh(); }

    // Inverse hyperbolic functions
    inline Complex asinh() const {
        return ((*this) + ((*this) * (*this) + Complex(1.0, 0.0)).sqrt()).log();
    }
    inline Complex acosh() const {
        return ((*this) + ((*this) * (*this) - Complex(1.0, 0.0)).sqrt()).log();
    }
    inline Complex atanh() const {
        return (Complex(1.0, 0.0) + *this).log() - (Complex(1.0, 0.0) - *this).log();
    }

    // Reciprocal hyperbolic functions
    inline Complex csch() const {
        Complex s = this->sinh();
        double d = s.re * s.re + s.im * s.im;
        if (d == 0.0 || !std::isfinite(d))
            throw std::runtime_error("Complex csch division by zero or invalid denominator");
        return Complex(1.0) / s;
    }
    inline Complex sech() const {
        Complex c = this->cosh();
        double d = c.re * c.re + c.im * c.im;
        if (d == 0.0 || !std::isfinite(d))
            throw std::runtime_error("Complex sech division by zero or invalid denominator");
        return Complex(1.0) / c;
    }
    inline Complex coth() const {
        Complex t = this->tanh();
        double d = t.re * t.re + t.im * t.im;
        if (d == 0.0 || !std::isfinite(d))
            throw std::runtime_error("Complex coth division by zero or invalid denominator");
        return Complex(1.0) / t;
    }

    // Inverse reciprocals of hyperbolic functions
    inline Complex acsch() const { return this->inverse().asinh(); }
    inline Complex asech() const { return this->inverse().acosh(); }
    inline Complex acoth() const { return this->inverse().atanh(); }

    // Miscellaneous
    double re, im;
    inline Complex conj() const { return {re, -im}; }
    inline double norm() const {
        // Use hypot for better numerical stability
        return std::hypot(re, im);
    }
    inline double abs() const {
        return std::hypot(re, im);
    }
    inline double arg() const { return std::atan2(im, re); }
    inline double real() const { return re; }
    inline double imag() const { return im; }
    inline double phase_deg() const { return arg() * (180.0 / PI); }
    // NaN/Inf checks
    inline bool is_nan() const { return std::isnan(re) || std::isnan(im); }
    inline bool is_inf() const { return std::isinf(re) || std::isinf(im); }

    // Static constants
    static const Complex I;
    static const Complex Zero;
    static const Complex One;
};

// Stream output
inline std::ostream& operator<<(std::ostream& os, const Complex& c)
{
    os << '(' << c.re;
    if (c.im < 0)
        os << c.im << "i)";
    else
        os << "+" << c.im << "i)";
    if (!std::isfinite(c.re) || !std::isfinite(c.im) || std::isnan(c.re) || std::isnan(c.im))
        os << " [invalid]";
    return os;
}

// Stream input
inline std::istream& operator>>(std::istream& is, Complex& c)
{
    is >> c.re >> c.im;
    if (!is || !std::isfinite(c.re) || !std::isfinite(c.im) || std::isnan(c.re) || std::isnan(c.im))
    {
        c.re = 0.0;
        c.im = 0.0;
    }
    return is;
}

// Define the static constants
const Complex Complex::I    = Complex(0.0, 1.0);
const Complex Complex::Zero = Complex(0.0, 0.0);
const Complex Complex::One  = Complex(1.0, 0.0);