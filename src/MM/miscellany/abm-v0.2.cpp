#include <cstdlib>
#include <ginac/integral.h>
#include <iostream>
#include <ginac/ginac.h>

template<typename T>
using Vec        = std::vector<T>;
using dVec       = Vec<double>;
using wVec       = Vec<size_t>;
using symbol     = GiNaC::symbol;
using expression = GiNaC::ex;
using list       = GiNaC::lst;
using integrate  = GiNaC::integral;

inline expression lagrange_basis(const symbol& x, const dVec& xs, size_t i)
{
    size_t     n          = xs.size();
    expression polynomial = 1;

    for ( size_t j=0; j<n; ++j )
    {
        if ( j!=i )
        {
            polynomial *= (x - xs[j]) / (xs[i] - xs[j]);
        }
    }
    return polynomial;
};

inline list integrated_lagrange_basis(const symbol& x, const dVec& xs, const double& h)
{
    size_t     n          = xs.size();
    double     x0         = xs[0];
    double     xi         = xs[n-1];
    double     xf         = xi + h;
    expression polynomial = 0;
    expression integral   = 0;
    list       integrals;
    for ( size_t k=0; k<n; ++k )
    {
        polynomial = 1;
        for ( size_t j=0; j<n; j++)
        {
            if ( j!=k )
            {
                polynomial *= (x - xs[j]) / (xs[k] - xs[j]);
            }
        }
        integral = integrate(x, xi, xf, polynomial).eval_integ().normal();
        integrals.append(integral);
    }
    return integrals;
}

inline list integrated_lagrange_basis_implicit(const symbol& x, const dVec& xs)
{
    size_t     n          = xs.size();
    double     x0         = xs[0];
    double     xi         = xs[n-2];
    double     xf         = xs[n-1];
    expression polynomial = 0;
    expression integral   = 0;
    list       integrals;
    for ( size_t k=0; k<n; ++k )
    {
        polynomial = 1;
        for ( size_t j=0; j<n; j++)
        {
            if ( j!=k )
            {
                polynomial *= (x - xs[j]) / (xs[k] - xs[j]);
            }
        }
        integral = integrate(x, xi, xf, polynomial).eval_integ().normal();
        integrals.append(integral);
    }
    return integrals;
}

inline expression make_lagrange(const symbol& x, const dVec& xs, const dVec& ys)
{
    size_t     n          = xs.size();
    expression polynomial = 0;

    for ( size_t i=0; i<n; ++i )
    {
        expression basis = 1;
        for ( size_t j=0; j<n; ++j )
        {
            if ( i!=j )
            {
                basis *= (x - xs[j]) / (xs[i] - xs[j]);
            }
        }
        polynomial += ys[i] * basis;
    }
    return polynomial;
};

int main()
{
    std::cout << "Initializing the symbols\n";
    dVec xs_e = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    dVec xs_i = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double h = 1.0;
    symbol x("x");
    auto coefs_e = integrated_lagrange_basis(x, xs_e, h);
    auto coefs_i = integrated_lagrange_basis_implicit(x, xs_i);
    std::cout << "The explicit coefficients:\n";
    std::cout << coefs_e << '\n';

    std::cout << "The implicit coefficients:\n";
    std::cout << coefs_i << '\n';
}
