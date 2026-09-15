"""
abm.py - Symbolic Adams-Bashforth-Moulton (ABM) Coefficient Generator

This script generates symbolic and numeric coefficients for Adams-Bashforth (predictor) and Adams-Moulton (corrector) methods of arbitrary order using sympy.

Usage:
    python abm.py <order>
    # or import and use abm_coeffs(order)

Example:
    python abm.py 4
    # Prints AB4/AM4 coefficients (symbolic and numeric)

Author: Hossein Karimi (with help from ChatGPT)
"""

import numpy as np
import sympy
import sys

t = sympy.symbols('t')

def P(t, s, m):
    res = 1
    for l in range(s):
        if l != m:
            res *= (t + s - 1 - l) / (m - l)
    return res

def Q(t, s, m):
    res = 1
    for l in range(s):
        if l != m:
            res *= (t + s - 2 - l) / (m - l)
    return res

def abm_coef_collector(s):
    # Adams-Bashforth (predictor)
    p = sympy.symarray('p', shape=s)
    b = sympy.symarray('b', shape=s)
    bnum = np.zeros(s, dtype=object)
    bdenom = np.zeros(s, dtype=object)
    ab = np.zeros(s, dtype=float)
    for m in range(s):
        p[m] = P(t, s, m)
        b[m] = sympy.Rational(sympy.integrate(p[m], [t, 0, 1]).simplify())
        bnum[m] = b[m].numerator
        bdenom[m] = b[m].denominator
        ab[m] = float(b[m])
    # Adams-Moulton (corrector)
    q = sympy.symarray('q', shape=s+1)
    c = sympy.symarray('c', shape=s+1)
    cnum = np.zeros(s+1, dtype=object)
    cdenom = np.zeros(s+1, dtype=object)
    abm = np.zeros(s+1, dtype=float)
    for n in range(s+1):
        q[n] = Q(t, s+1, n)
        c[n] = sympy.Rational(sympy.integrate(q[n], [t, 0, 1]).simplify())
        cnum[n] = c[n].numerator
        cdenom[n] = c[n].denominator
        abm[n] = float(c[n])
    return ab, abm, bnum, bdenom, cnum, cdenom

def print_abm(s):
    ab, abm, bnum, bdenom, cnum, cdenom = abm_coef_collector(s)
    print(f"Adams-Bashforth-Moulton coefficients for order {s}:")
    print("\nAdams-Bashforth (AB{}) coefficients:".format(s))
    for i in range(s):
        print(f"  b[{i}] = {bnum[i]}/{bdenom[i]}")
    print("Numeric (float):")
    print(list(ab))
    print("As fractions:")
    print([f"{bnum[i]}/{bdenom[i]}" for i in range(s)])
    print("\nAdams-Moulton (AM{}) coefficients:".format(s))
    for i in range(s+1):
        print(f"  b[{i}] = {cnum[i]}/{cdenom[i]}")
    print("Numeric (float):")
    print(list(abm))
    print("As fractions:")
    print([f"{cnum[i]}/{cdenom[i]}" for i in range(s+1)])
    print(f"\nNumerators (AB): {list(bnum)}")
    print(f"Denominators (AB): {list(bdenom)}")
    print(f"Numerators (AM): {list(cnum)}")
    print(f"Denominators (AM): {list(cdenom)}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python abm.py <order>")
        print("Example: python abm.py 4")
    else:
        s = int(sys.argv[1])
        print_abm(s)