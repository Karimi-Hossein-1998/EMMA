#include <cstdlib>
#include <iostream>
#include <ginac/ginac.h>

template<typename T>
using Vec        = std::vector<T>;
using dVec       = Vec<double>;
using symbol     = GiNaC::symbol;
using expression = GiNaC::ex;

inline expression make_lagrange(const symbol& x, const dVec& xs, const dVec& ys)
{
    size_t     n          = xs.size();
    expression polynomial = 0;

    for ( size_t i=0; i<n; ++i )
    {
        expression basis = 1;
        for ( size_t j=0; j<n; ++j )
        {
            if ( i != j )
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
    dVec   xs = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    dVec   ys = {0.0, 1.0, -4.0, 9.0, -16.0, 25.0};
    symbol x("x");

    expression raw_form      = make_lagrange(x, xs, ys);
    expression expanded_form = raw_form.expand();
    std::cout << "Raw Form:\n" << raw_form << '\n';
    std::cout << "Expanded Form:\n" << expanded_form << '\n';

    std::cout << "P(0) = " << expanded_form.subs(x == 0.0) << '\n';
    std::cout << "P(0.5) = " << expanded_form.subs(x == 0.5) << '\n';
    std::cout << "P(1) = " << expanded_form.subs(x == 1.0) << '\n';
    std::cout << "P(1.5) = " << expanded_form.subs(x == 1.5) << '\n';
    std::cout << "P(2) = " << expanded_form.subs(x == 2.0) << '\n';
    std::cout << "P(2.5) = " << expanded_form.subs(x == 2.5) << '\n';
    std::cout << "P(3) = " << expanded_form.subs(x == 3.0) << '\n';
    std::cout << "P(3.5) = " << expanded_form.subs(x == 3.5) << '\n';
    std::cout << "P(4) = " << expanded_form.subs(x == 4.0) << '\n';
    std::cout << "P(4.5) = " << expanded_form.subs(x == 4.5) << '\n';
    std::cout << "P(5) = " << expanded_form.subs(x == 5.0) << '\n';
}
