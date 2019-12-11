import os
import sympy as sym
import sys
tree = os.path.abspath(os.path.join(__file__, '../..'))
pytpath = os.path.join(tree, 'PyTransport', 'PyTransport')
sys.path.append(pytpath)
import PyTransSetup
from sympy import Rational as rat

# p[0] := Brane tension, T3
# p[1] := Flux source, a0
# p[2] := Tip warp fact, Q
# p[3] := Size of field space phi_uv
# p[4] := Distance sources, V_0
# p[5:32] := CLM coeffecients
# f[0] := radial coordinate, x
# f[1] := ang. coord, theta_1
# f[2] := ang. coord, theta_2
# f[3] := ang. coord, phi_1
# f[4] := ang. coord, phi_2
# f[5] := ang. coord, psi

p = sym.symarray('p',32)
f = sym.symarray('f',6)
s = [
	2*p[0]*p[1]**4*(-27*p[0]*p[1]**4/(32*sym.pi**2*p[3]**4*f[0]**4) + 1),
	p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])/3,
	9*p[5]*p[2]*p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[2])/(4*sym.pi**(3/2)),
	3*sym.sqrt(15)*p[6]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(10))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2])**2 - 1)/(8*sym.pi**(3/2)),
	9*p[7]*p[2]*p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1])/(4*sym.pi**(3/2)),
	9*sym.sqrt(3)*p[8]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(7))*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1])*sym.cos(f[2])/(4*sym.pi**(3/2)),
	9*sym.sqrt(5)*p[9]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(13))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2])**2 - 1)*sym.cos(f[1])/(8*sym.pi**(3/2)),
	3*sym.sqrt(15)*p[10]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(10))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1])**2 - 1)/(8*sym.pi**(3/2)),
	9*sym.sqrt(5)*p[11]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(13))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1])**2 - 1)*sym.cos(f[2])/(8*sym.pi**(3/2)),
	9*sym.sqrt(3)*p[2]*p[3]**2*f[0]**3*(-p[22]*sym.sin(f[3] - f[4] + f[5]) + p[12]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])*sym.sin(f[2])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
	9*sym.sqrt(5)*p[2]*p[3]**2*f[0]**5*(-p[23]*sym.sin(f[3] - f[4] + f[5]) + p[13]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[2]) + sym.sin(2*f[2]))*sym.sin(f[1])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
	9*sym.sqrt(5)*p[2]*p[3]**2*f[0]**5*(p[24]*sym.sin(f[3] - f[4] + f[5]) - p[14]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[1]) - sym.sin(2*f[1]))*sym.sin(f[2])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(3/2)*(-p[25]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[15]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(2*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	3*sym.sqrt(6)*p[2]*p[3]**2*f[0]**(7/2)*(-p[26]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[16]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	9*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(241)/2)*(-p[27]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[17]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(-5*sym.sin(f[2])**2 + 2*sym.cos(f[2]) + 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	3*sym.sqrt(6)*p[2]*p[3]**2*f[0]**(7/2)*(-p[28]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[18]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(193)/2)*(-p[29]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[19]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(9/2)*(-p[30]*sym.sin(sym.Rational(3, 2)*f[3] - sym.Rational(3, 2)*f[4] + 3*f[5]/2) + p[20]*sym.cos(sym.Rational(3, 2)*f[3] - sym.Rational(3, 2)*f[4] + 3*f[5]/2))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])**(3/2)*sym.sin(f[2])**(3/2)*sym.tan(f[2]/2)**(3/2)/(4*sym.pi**(3/2)*sym.tan(f[1]/2)**(3/2)),
	9*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(241)/2)*(p[31]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) - p[21]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(5*sym.sin(f[1])**2 + 2*sym.cos(f[1]) - 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
	0
]
V = sum(s)

G = sym.Matrix(
	[
		[f[0]**-2, 0, 0, 0, 0, 0],
		[0, rat(1, 6), 0, 0, 0, 0],
		[0, 0, rat(1, 6), 0, 0, 0],
		[0, 0, 0, rat(1, 9)*sym.cos(f[1])**2 + rat(1, 6)*sym.sin(f[1])**2, sym.cos(f[1])*sym.cos(f[2])*rat(1, 9), sym.cos(f[1])*rat(1, 9)],
		[0, 0, 0, sym.cos(f[1])*sym.cos(f[2])*rat(1, 9), rat(1, 9)*sym.cos(f[2])**2 + rat(1, 6)*sym.sin(f[2])**2, sym.cos(f[2])*rat(1, 9)],
		[0, 0, 0, sym.cos(f[1])*rat(1, 9), sym.cos(f[2])*rat(1, 9), rat(1, 9)]
	]
)
G *= p[3]**2*f[0]**2
PyTransSetup.potential(V,6,32,True,G,silent=False)
PyTransSetup.compileName('agarwal_6pt0_HALFSIMPLE',True)
