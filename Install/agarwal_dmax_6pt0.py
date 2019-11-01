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
	2*p[0]*p[1]**4*(-0.0854897486982225*p[0]*p[1]**4/(p[3]**4*f[0]**4) + 1),
	0.333333333333333*p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4]),
	0.404071024781625*p[5]*p[2]*p[3]**2*f[0]**2.0*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[2]),
	p[6]*p[2]*p[3]**2*f[0]**4.32455532033676*(2*p[0]*p[1]**4 + p[4])*(0.782480174832099*sym.cos(f[2])**2 - 0.260826724944033),
	0.404071024781625*p[7]*p[2]*p[3]**2*f[0]**2.0*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1]),
	0.699871544788197*p[8]*p[2]*p[3]**2*f[0]**3.29150262212918*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1])*sym.cos(f[2]),
	0.451765139574858*p[9]*p[2]*p[3]**2*f[0]**5.21110255092798*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2])**2 - 1)*sym.cos(f[1]),
	p[10]*p[2]*p[3]**2*f[0]**4.32455532033676*(2*p[0]*p[1]**4 + p[4])*(0.782480174832099*sym.cos(f[1])**2 - 0.260826724944033),
	1.0*p[11]*p[2]*p[3]**2*f[0]**5.21110255092798*(2*p[0]*p[1]**4 + p[4])*(1.35529541872457*sym.cos(f[1])**2 - 0.451765139574858)*sym.cos(f[2]),
	0.349935772394099*p[2]*p[3]**2*f[0]**3.0*(-p[22]*sym.sin(1.0*f[3] - 1.0*f[4] + 1.0*f[5]) + p[12]*sym.cos(1.0*f[3] - 1.0*f[4] + 1.0*f[5]))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])*sym.sin(f[2])*sym.tan(f[2]/2)/sym.tan(f[1]/2),
	0.451765139574858*p[2]*p[3]**2*f[0]**5.0*(-p[23]*sym.sin(1.0*f[3] - 1.0*f[4] + 1.0*f[5]) + p[13]*sym.cos(1.0*f[3] - 1.0*f[4] + 1.0*f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[2]) + sym.sin(2*f[2]))*sym.sin(f[1])*sym.tan(f[2]/2)/sym.tan(f[1]/2),
	0.451765139574858*p[2]*p[3]**2*f[0]**5.0*(p[24]*sym.sin(1.0*f[3] - 1.0*f[4] + 1.0*f[5]) - p[14]*sym.cos(1.0*f[3] - 1.0*f[4] + 1.0*f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[1]) - sym.sin(2*f[1]))*sym.sin(f[2])*sym.tan(f[2]/2)/sym.tan(f[1]/2),
	0.466581029858798*p[2]*p[3]**2*f[0]**1.5*(-p[25]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) + p[15]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
	0.329922610186159*p[2]*p[3]**2*f[0]**3.5*(-p[26]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) + p[16]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
	0.404071024781625*p[2]*p[3]**2*f[0]**5.76208734813001*(-p[27]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) + p[17]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*(-5*sym.sin(f[2])**2 + 2*sym.cos(f[2]) + 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
	0.329922610186159*p[2]*p[3]**2*f[0]**3.5*(-p[28]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) + p[18]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
	0.233290514929399*p[2]*p[3]**2*f[0]**4.9462219947249*(-p[29]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) + p[19]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
	0.233290514929399*p[2]*p[3]**2*f[0]**4.5*(-p[30]*sym.sin(1.5*f[3] - 1.5*f[4] + 1.5*f[5]) + p[20]*sym.cos(1.5*f[3] - 1.5*f[4] + 1.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])**(3/2)*sym.sin(f[2])**(3/2)*sym.tan(f[2]/2)**(3/2)/sym.tan(f[1]/2)**(3/2),
	0.404071024781625*p[2]*p[3]**2*f[0]**5.76208734813001*(p[31]*sym.sin(0.5*f[3] - 0.5*f[4] + 0.5*f[5]) - p[21]*sym.cos(0.5*f[3] - 0.5*f[4] + 0.5*f[5]))*(2*p[0]*p[1]**4 + p[4])*(5*sym.sin(f[1])**2 + 2*sym.cos(f[1]) - 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/sym.sqrt(sym.tan(f[1]/2)),
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
PyTransSetup.potential(V,6,32,False,G,silent=False)
PyTransSetup.compileName('agarwal_dmax_6pt0',True)