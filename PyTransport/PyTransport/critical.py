import numpy as np
import sympy as sym

x, t1, t2, p1, p2, ps = sym.symbols("x a b p q z")
m1, m2, m3 = sym.symarray("m", 3)

trig_arg = m1*p1 + m2*p2 + m3/2
rtrig, itrig = 2*sym.cos(trig_arg), 2*sym.sin(trig_arg)

def z(m1_, m2_, m3_):
    
    str_ = "C({}|{}|{})".format(m1_, m2_, m3_)
    
    Cr = sym.symbols("r" + str_)
    Ci = sym.symbols("i" + str_)
    
    subs = dict(zip([m1, m2, m3], [m1_, m2_, m3_]))
    
    exp = Cr*rtrig.subs(subs) + Ci*itrig.subs(subs)
    
    return exp

half = sym.Rational(1, 2)

data = [
    [3*half, half, half, -1, -half,  half],
    [3*half, half, half, +1, +half, -half],
    [2, 1, 0, 0, 0, 0],
    [2, 0, 1, 0, 0, 0],
    [3, 1, 1, -2, -1, +1],
    [3, 1, 1, +2, +1, -1],
    [-2+sym.sqrt(28), 1, 1, 0, 0, 0]
]

exps = [
    sym.sqrt((1+sym.cos(t1)*(1-sym.cos(t2)))),
    sym.cos(t1) + sym.cos(t2),
    (1+sym.cos(t1)*(1-sym.cos(t2))),
    sym.cos(t1)*sym.cos(t2)
]

terms = [
    x**sym.Rational(3, 2) * 3 * sym.sqrt(3) / (4 * sym.pi **sym.Rational(3, 2) * (z(*data[0][3:]) + z(*data[1][3:]))),
    x**2 * 9 / (4 * sym.pi**sym.Ration(3, 2)) *
]