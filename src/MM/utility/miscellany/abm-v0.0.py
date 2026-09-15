import numpy as np
import sympy

# s = np.int64(input("Please enter the degree of predictor-corrector method you want to be calculated: "))
s = 1
m = 1
t = sympy.symbols('t')
P = sympy.Function('P')(t, s, m)
def P(t,s,m):
    res = 1
    for l in range(s):
        if l!= m:
            res*= (t+s-1-l)/(m-l)
    return res
Q = sympy.Function('Q')(t, s, m)
def Q(t,s,m):
    res = 1
    # This loop could be from 1 to s or 0 to s-1, I chose the latter
    for l in range(s):
        if l!= m:
            res *= (t+s-2-l)/(m-l)
    return res

# This Function Collects the Adams-Bashforth and Adams-Moulton coefficients
def abm_coef_collector(s):
    # 'p's are the polynomial interpolations of the predictor
    p = sympy.symarray(('p'), shape=s)
    b = sympy.symarray(('b'),shape=s)
    bnum = np.zeros(s).reshape((1,s))
    bdenom = np.zeros(s).reshape((1,s))
    ab = np.zeros(s).reshape((1,s))
    for m in range(s):
        p[m]=P(t,s,m)
        b[m]=sympy.Rational(sympy.integrate(p[m],[t,0,1]).simplify())
        bnum[0,m] = b[m].numerator
        bdenom[0,m] = b[m].denominator
        ab[0,m] = 1.*np.int64(bnum[0,m])/np.int64(bdenom[0,m])
    # 'q's are the polynomial interpolations of the corrector
    q = sympy.symarray(('q'), shape=s+1)
    c = sympy.symarray(('c'),shape=s+1)
    cnum = np.zeros(s+1).reshape((1,s+1))
    cdenom = np.zeros(s+1).reshape((1,s+1))
    abm = np.zeros(s+1).reshape((1,s+1))
    for n in range(s+1):
        q[n]=Q(t,s+1,n)
        c[n]=sympy.Rational(sympy.integrate(q[n],[t,0,1]).simplify())
        cnum[0,n] = c[n].numerator
        cdenom[0,n] = c[n].denominator
        abm[0,n] = 1.*np.int64(cnum[0,n])/np.int64(cdenom[0,n])
    return ab,abm,bnum,bdenom,cnum,cdenom
ab,abm,bnum,bdenom,cnum,cdenom = abm_coef_collector(2)
print(ab,'\n')
print(abm,'\n')
print(bnum,'\n')
print(bdenom,'\n')
print(cnum,'\n')
print(cdenom,'\n')