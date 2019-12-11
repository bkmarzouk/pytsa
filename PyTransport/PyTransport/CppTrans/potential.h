//#This file is part of PyTransport.

//#PyTransport is free software: you can redistribute it and/or modify
//#it under the terms of the GNU General Public License as published by
//#the Free Software Foundation, either version 3 of the License, or
//#(at your option) any later version.

//#PyTransport is distributed in the hope that it will be useful,
//#but WITHOUT ANY WARRANTY; without even the implied warranty of
//#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//#GNU General Public License for more details.

//#You should have received a copy of the GNU General Public License
//#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.

// This file contains a prototype of the potential.h file of PyTransport -- it is edited by the PyTransScripts module

#ifndef POTENTIAL_H  // Prevents the class being re-defined
#define POTENTIAL_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

// #Rewrite
// Potential file rewriten at Wed Dec 11 14:37:33 2019

class potential
{
private:
	int nF; // field number
	int nP; // params number which definFs potential
    
    
public:
	// flow constructor
	potential()
	{
// #FP
nF=6;
nP=32;

//        p.resize(nP);
        
// pdef

    }
	
    //void setP(vector<double> pin){
    //    p=pin;
    //}
	//calculates V()
	double V(vector<double> f, vector<double> p)
	{
		double sum ;
        
// Pot
  auto x0 = p[0]*std::pow(p[1], 4);
  auto x1 = 2*x0;
  auto x2 = std::pow(p[3], 2)*(p[4] + x1);
  auto x3 = std::cos(f[1]);
  auto x4 = std::pow(x3, 2);
  auto x5 = p[2]*x2;
  auto x6 = std::pow(f[0], 4.3245553203367599)*x5;
  auto x7 = std::cos(f[2]);
  auto x8 = std::pow(x7, 2);
  auto x9 = x5*x7;
  auto x10 = 0.40407102478162499*std::pow(f[0], 2.0);
  auto x11 = x3*x5;
  auto x12 = std::pow(f[0], 5.21110255092798);
  auto x13 = 1.0*f[3] - 1.0*f[4] + 1.0*f[5];
  auto x14 = std::cos(x13);
  auto x15 = std::sin(x13);
  auto x16 = std::sin(f[1]);
  auto x17 = std::tan((1.0/2.0)*f[1]);
  auto x18 = 1.0/x17;
  auto x19 = std::sin(f[2]);
  auto x20 = std::tan((1.0/2.0)*f[2]);
  auto x21 = x18*x19*x20*x5;
  auto x22 = x16*x21;
  auto x23 = 1.5*f[3] - 1.5*f[4] + 1.5*f[5];
  auto x24 = 0.5*f[3] - 0.5*f[4] + 0.5*f[5];
  auto x25 = std::cos(x24);
  auto x26 = std::sin(x24);
  auto x27 = std::sqrt(x16)*std::sqrt(x19)*std::sqrt(x20)*x5/std::sqrt(x17);
  auto x28 = 0.45176513957485798*std::pow(f[0], 5.0);
  auto x29 = x27*(3*x7 + 1);
  auto x30 = 0.32992261018615898*std::pow(f[0], 3.5);
  auto x31 = 3*x3 - 1;
  auto x32 = 0.40407102478162499*std::pow(f[0], 5.7620873481300103)*x27;
  sum=0.33333333333333298*std::pow(f[0], 2)*x2 + 0.46658102985879801*std::pow(f[0], 1.5)*x27*(p[15]*x25 - p[25]*x26) + 0.34993577239409901*std::pow(f[0], 3.0)*x22*(p[12]*x14 - p[22]*x15) + 0.69987154478819702*std::pow(f[0], 3.2915026221291801)*p[8]*x3*x9 + 0.23329051492939901*std::pow(f[0], 4.5)*x22*(p[20]*std::cos(x23) - p[30]*std::sin(x23)) + 0.23329051492939901*std::pow(f[0], 4.9462219947249002)*x29*x31*(p[19]*x25 - p[29]*x26) + p[10]*x6*(0.78248017483209897*x4 - 0.26082672494403297) + 1.0*p[11]*x12*x9*(1.35529541872457*x4 - 0.45176513957485798) + p[5]*x10*x9 + p[6]*x6*(0.78248017483209897*x8 - 0.26082672494403297) + p[7]*x10*x11 + 0.45176513957485798*p[9]*x11*x12*(3*x8 - 1) + x1*(1 - 0.085489748698222498*x0/(std::pow(f[0], 4)*std::pow(p[3], 4))) + x16*x18*x20*x28*x5*(x19 + std::sin(2*f[2]))*(p[13]*x14 - p[23]*x15) + x21*x28*(x16 - std::sin(2*f[1]))*(-p[14]*x14 + p[24]*x15) + x27*x30*x31*(p[18]*x25 - p[28]*x26) + x29*x30*(p[16]*x25 - p[26]*x26) + x32*(p[17]*x25 - p[27]*x26)*(-5*std::pow(x19, 2) + 2*x7 + 4) + x32*(-p[21]*x25 + p[31]*x26)*(5*std::pow(x16, 2) + 2*x3 - 4);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot
  auto x0 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x1 = std::cos(f[1]);
  auto x2 = std::pow(x1, 2);
  auto x3 = p[2]*x0;
  auto x4 = p[10]*x3;
  auto x5 = 4.3245553203367599*std::pow(f[0], 3.3245553203367599);
  auto x6 = std::cos(f[2]);
  auto x7 = std::pow(x6, 2);
  auto x8 = p[6]*x3;
  auto x9 = p[5]*x3;
  auto x10 = 0.80814204956324998*std::pow(f[0], 1.0);
  auto x11 = p[7]*x3;
  auto x12 = x1*x3;
  auto x13 = p[8]*x6;
  auto x14 = std::pow(f[0], 4.21110255092798);
  auto x15 = p[9]*(3*x7 - 1);
  auto x16 = x3*x6;
  auto x17 = p[11]*(1.35529541872457*x2 - 0.45176513957485798);
  auto x18 = std::tan((1.0/2.0)*f[1]);
  auto x19 = 1.0/x18;
  auto x20 = std::tan((1.0/2.0)*f[2]);
  auto x21 = x19*x20;
  auto x22 = 2.2588256978742898*std::pow(f[0], 4.0)*x21;
  auto x23 = 2*f[1];
  auto x24 = std::sin(f[1]);
  auto x25 = x24 - std::sin(x23);
  auto x26 = 1.0*f[3] - 1.0*f[4] + 1.0*f[5];
  auto x27 = std::sin(x26);
  auto x28 = std::cos(x26);
  auto x29 = -p[14]*x28 + p[24]*x27;
  auto x30 = std::sin(f[2]);
  auto x31 = x3*x30;
  auto x32 = x29*x31;
  auto x33 = x25*x32;
  auto x34 = x24*x3;
  auto x35 = p[13]*x28 - p[23]*x27;
  auto x36 = 2*f[2];
  auto x37 = x30 + std::sin(x36);
  auto x38 = x35*x37;
  auto x39 = std::pow(f[0], 3.5);
  auto x40 = 1.5*f[3] - 1.5*f[4] + 1.5*f[5];
  auto x41 = std::cos(x40);
  auto x42 = std::sin(x40);
  auto x43 = p[20]*x41 - p[30]*x42;
  auto x44 = x21*x43;
  auto x45 = x30*x34;
  auto x46 = std::pow(f[0], 2.0);
  auto x47 = p[12]*x28 - p[22]*x27;
  auto x48 = x21*x47;
  auto x49 = std::sqrt(x24);
  auto x50 = 0.5*f[3] - 0.5*f[4] + 0.5*f[5];
  auto x51 = std::cos(x50);
  auto x52 = std::sin(x50);
  auto x53 = p[15]*x51 - p[25]*x52;
  auto x54 = x49*x53;
  auto x55 = std::sqrt(x30);
  auto x56 = std::pow(x18, -1.0/2.0);
  auto x57 = std::sqrt(x20);
  auto x58 = x55*x56*x57;
  auto x59 = x3*x58;
  auto x60 = 0.69987154478819702*x59;
  auto x61 = 1.1547291356515563*std::pow(f[0], 2.5);
  auto x62 = 3*x6 + 1;
  auto x63 = p[16]*x51 - p[26]*x52;
  auto x64 = x49*x63;
  auto x65 = x62*x64;
  auto x66 = p[18]*x51 - p[28]*x52;
  auto x67 = x59*x66;
  auto x68 = 3*x1 - 1;
  auto x69 = x49*x68;
  auto x70 = -5*std::pow(x30, 2) + 2*x6 + 4;
  auto x71 = p[17]*x51 - p[27]*x52;
  auto x72 = x49*x59;
  auto x73 = x71*x72;
  auto x74 = 2.3282925396401293*std::pow(f[0], 4.7620873481300103);
  auto x75 = 2*x1 + 5*std::pow(x24, 2) - 4;
  auto x76 = -p[21]*x51 + p[31]*x52;
  auto x77 = x72*x76;
  auto x78 = p[19]*x51 - p[29]*x52;
  auto x79 = x62*x72;
  auto x80 = x68*x79;
  auto x81 = 0.40407102478162499*x46;
  auto x82 = std::pow(f[0], 5.21110255092798);
  auto x83 = 0.45176513957485798*x34;
  auto x84 = 0.69987154478819702*std::pow(f[0], 3.2915026221291801);
  auto x85 = x1*x24;
  auto x86 = 1.5649603496641979*std::pow(f[0], 4.3245553203367599);
  auto x87 = x12*x30;
  auto x88 = 0.34993577239409901*std::pow(f[0], 3.0);
  auto x89 = x48*x88;
  auto x90 = 0.23329051492939901*std::pow(f[0], 4.5);
  auto x91 = x44*x90;
  auto x92 = std::pow(f[0], 5.0);
  auto x93 = 0.45176513957485798*x92;
  auto x94 = x21*x93;
  auto x95 = std::pow(x24, 3.0/2.0);
  auto x96 = 0.98976783055847695*x39;
  auto x97 = x45*x88;
  auto x98 = x47*x97;
  auto x99 = std::pow(x18, 2);
  auto x100 = x20*(-1.0/2.0*x99 - 1.0/2.0)/x99;
  auto x101 = std::pow(f[0], 1.5);
  auto x102 = 0.23329051492939901*x101;
  auto x103 = x12*x58/x49;
  auto x104 = x45*x90;
  auto x105 = x104*x43;
  auto x106 = 0.46658102985879801*x101;
  auto x107 = x106*x54*x55;
  auto x108 = x3*x57;
  auto x109 = (-1.0/4.0*x99 - 1.0/4.0)/std::pow(x18, 3.0/2.0);
  auto x110 = x108*x109;
  auto x111 = x33*x93;
  auto x112 = x83*x92;
  auto x113 = x112*x38;
  auto x114 = std::pow(f[0], 5.7620873481300103);
  auto x115 = 0.40407102478162499*x114;
  auto x116 = std::pow(f[0], 4.9462219947249002);
  auto x117 = x116*x78;
  auto x118 = x117*x62;
  auto x119 = 0.16496130509307949*x39;
  auto x120 = x103*x119;
  auto x121 = x108*x64;
  auto x122 = 0.32992261018615898*x39;
  auto x123 = x122*x55;
  auto x124 = x66*x69;
  auto x125 = x123*x124;
  auto x126 = x70*x71;
  auto x127 = 0.20203551239081249*x114;
  auto x128 = x103*x127;
  auto x129 = x75*x76;
  auto x130 = x126*x49;
  auto x131 = x110*x55;
  auto x132 = x115*x131;
  auto x133 = x129*x49;
  auto x134 = 0.1166452574646995*x118;
  auto x135 = 0.23329051492939901*x118*x69;
  auto x136 = x30*x6;
  auto x137 = std::pow(x20, 2);
  auto x138 = x19*((1.0/2.0)*x137 + 1.0/2.0);
  auto x139 = x16*x24;
  auto x140 = x25*x94;
  auto x141 = x112*x21;
  auto x142 = std::pow(x30, 3.0/2.0)*x56;
  auto x143 = x16*x56*x57/x55;
  auto x144 = x3*x56*((1.0/4.0)*x137 + 1.0/4.0)/x57;
  auto x145 = x119*x143;
  auto x146 = x127*x143;
  auto x147 = x144*x55;
  auto x148 = x115*x147;
  auto x149 = 1.0*x27;
  auto x150 = p[12]*x149;
  auto x151 = 1.0*x28;
  auto x152 = p[22]*x151;
  auto x153 = x21*x97;
  auto x154 = 1.5*p[20]*x42;
  auto x155 = 1.5*p[30]*x41;
  auto x156 = x104*x21;
  auto x157 = 0.5*x52;
  auto x158 = p[15]*x157;
  auto x159 = 0.5*x51;
  auto x160 = p[25]*x159;
  auto x161 = x106*x72;
  auto x162 = p[14]*x149;
  auto x163 = p[24]*x151;
  auto x164 = x140*x31;
  auto x165 = p[13]*x149;
  auto x166 = p[23]*x151;
  auto x167 = x141*x37;
  auto x168 = p[16]*x157;
  auto x169 = p[26]*x159;
  auto x170 = x122*x79;
  auto x171 = p[18]*x157;
  auto x172 = p[28]*x159;
  auto x173 = x122*x68*x72;
  auto x174 = p[21]*x157;
  auto x175 = p[31]*x159;
  auto x176 = x115*x72;
  auto x177 = x176*x75;
  auto x178 = p[17]*x157;
  auto x179 = p[27]*x159;
  auto x180 = x176*x70;
  auto x181 = p[19]*x157;
  auto x182 = p[29]*x159;
  auto x183 = 0.23329051492939901*x116*x80;
  auto x184 = x153*(-x150 - x152) + x156*(-x154 - x155) + x161*(-x158 - x160) + x164*(x162 + x163) + x167*(-x165 - x166) + x170*(-x168 - x169) + x173*(-x171 - x172) + x177*(x174 + x175) + x180*(-x178 - x179) + x183*(-x181 - x182);

 sum[0]=std::sqrt(f[0])*x54*x60 + 0.66666666666666596*f[0]*x0 + 2.3036290248239504*std::pow(f[0], 2.2915026221291801)*x12*x13 + 1.1539066761044912*std::pow(f[0], 3.9462219947249002)*x78*x80 + x1*x10*x11 + x10*x6*x9 + 2.3541944712588774*x12*x14*x15 + 5.21110255092798*x14*x16*x17 + x22*x33 + x22*x34*x38 + 1.0498073171822955*x39*x44*x45 + x4*x5*(0.78248017483209897*x2 - 0.26082672494403297) + 1.049807317182297*x45*x46*x48 + x5*x8*(0.78248017483209897*x7 - 0.26082672494403297) + x59*x61*x65 + x61*x67*x69 + x70*x73*x74 + x74*x75*x77 + 0.68391798958577998*std::pow(p[0], 2)*std::pow(p[1], 8)/(std::pow(f[0], 5)*std::pow(p[3], 4));

 sum[1]=-2.7105908374491401*p[11]*x16*x82*x85 + x100*x105 + x100*x111 + x100*x113 + x100*x98 + x102*x103*x53 + x103*x134*x68 + x107*x110 + x109*x121*x123*x62 - x11*x24*x81 + x110*x125 + x115*x77*(-2*x24 + 10*x85) - x118*x60*x95 + x12*x38*x94 + x120*x62*x63 + x120*x66*x68 + x126*x128 + x128*x129 - x13*x34*x84 + x130*x132 + x131*x135 + x132*x133 - x15*x82*x83 + x32*x94*(x1 - 2*std::cos(x23)) - x4*x85*x86 - x67*x95*x96 + x87*x89 + x87*x91;

 sum[2]=-p[8]*x84*x87 - 2.7105908374491481*p[9]*x12*x136*x82 + x102*x143*x54 + x105*x138 + x107*x144 - 0.69987154478819702*x108*x117*x142*x69 + x111*x138 + x113*x138 + x115*x73*(-10*x136 - 2*x30) - x121*x142*x96 + x123*x144*x65 + x124*x145 + x125*x144 + x130*x146 + x130*x148 + x133*x146 + x133*x148 + x134*x143*x69 + x135*x147 - x136*x8*x86 + x138*x98 + x139*x89 + x139*x91 + x140*x16*x29 + x141*x35*(x6 + 2*std::cos(x36)) + x145*x65 - 1.0*x17*x31*x82 - x30*x81*x9;

 sum[3]=x184;

 sum[4]=x153*(x150 + x152) + x156*(x154 + x155) + x161*(x158 + x160) + x164*(-x162 - x163) + x167*(x165 + x166) + x170*(x168 + x169) + x173*(x171 + x172) + x177*(-x174 - x175) + x180*(x178 + x179) + x183*(x181 + x182);

 sum[5]=x184;
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot
  auto x0 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x1 = 0.66666666666666596*x0;
  auto x2 = std::pow(p[0], 2)*std::pow(p[1], 8)/std::pow(p[3], 4);
  auto x3 = std::cos(f[2]);
  auto x4 = p[2]*x0;
  auto x5 = 0.80814204956324998*x4;
  auto x6 = p[5]*x5;
  auto x7 = x3*x6;
  auto x8 = std::cos(f[1]);
  auto x9 = p[7]*x5;
  auto x10 = x8*x9;
  auto x11 = 14.377223398316216*std::pow(f[0], 2.3245553203367599);
  auto x12 = std::pow(x8, 2);
  auto x13 = p[10]*x4;
  auto x14 = x13*(0.78248017483209897*x12 - 0.26082672494403297);
  auto x15 = std::pow(x3, 2);
  auto x16 = p[6]*x4;
  auto x17 = x16*(0.78248017483209897*x15 - 0.26082672494403297);
  auto x18 = std::pow(f[0], 3.21110255092798);
  auto x19 = x4*x8;
  auto x20 = p[9]*(3*x15 - 1);
  auto x21 = x19*x20;
  auto x22 = x3*x4;
  auto x23 = p[11]*(1.35529541872457*x12 - 0.45176513957485798);
  auto x24 = x22*x23;
  auto x25 = p[8]*x3;
  auto x26 = x19*x25;
  auto x27 = std::pow(f[0], 2.5);
  auto x28 = std::tan((1.0/2.0)*f[1]);
  auto x29 = 1.0/x28;
  auto x30 = std::tan((1.0/2.0)*f[2]);
  auto x31 = x29*x30;
  auto x32 = 1.5*f[3] - 1.5*f[4] + 1.5*f[5];
  auto x33 = std::cos(x32);
  auto x34 = p[20]*x33;
  auto x35 = std::sin(x32);
  auto x36 = p[30]*x35;
  auto x37 = x34 - x36;
  auto x38 = std::sin(f[2]);
  auto x39 = std::sin(f[1]);
  auto x40 = x39*x4;
  auto x41 = x38*x40;
  auto x42 = x37*x41;
  auto x43 = x31*x42;
  auto x44 = std::pow(f[0], 3.0);
  auto x45 = x31*x44;
  auto x46 = 9.0353027914971591*x45;
  auto x47 = 2*f[1];
  auto x48 = std::sin(x47);
  auto x49 = x39 - x48;
  auto x50 = 1.0*f[3] - 1.0*f[4] + 1.0*f[5];
  auto x51 = std::sin(x50);
  auto x52 = p[24]*x51;
  auto x53 = std::cos(x50);
  auto x54 = p[14]*x53;
  auto x55 = x52 - x54;
  auto x56 = x38*x4;
  auto x57 = x55*x56;
  auto x58 = x49*x57;
  auto x59 = p[13]*x53;
  auto x60 = p[23]*x51;
  auto x61 = x59 - x60;
  auto x62 = 2*f[2];
  auto x63 = std::sin(x62);
  auto x64 = x38 + x63;
  auto x65 = x61*x64;
  auto x66 = x40*x65;
  auto x67 = std::pow(f[0], 1.0);
  auto x68 = x38*x67;
  auto x69 = p[12]*x53;
  auto x70 = p[22]*x51;
  auto x71 = x69 - x70;
  auto x72 = x31*x71;
  auto x73 = x40*x72;
  auto x74 = std::sqrt(f[0]);
  auto x75 = std::sqrt(x39);
  auto x76 = std::pow(x28, -1.0/2.0);
  auto x77 = std::sqrt(x30);
  auto x78 = x76*x77;
  auto x79 = x4*x78;
  auto x80 = x75*x79;
  auto x81 = std::sqrt(x38);
  auto x82 = 0.5*f[3] - 0.5*f[4] + 0.5*f[5];
  auto x83 = std::cos(x82);
  auto x84 = p[15]*x83;
  auto x85 = std::sin(x82);
  auto x86 = p[25]*x85;
  auto x87 = x84 - x86;
  auto x88 = x81*x87;
  auto x89 = 0.34993577239409851*x88;
  auto x90 = 11.087532445765751*std::pow(f[0], 3.7620873481300103);
  auto x91 = std::pow(x38, 2);
  auto x92 = 2*x3;
  auto x93 = -5*x91 + x92 + 4;
  auto x94 = p[17]*x83;
  auto x95 = p[27]*x85;
  auto x96 = x94 - x95;
  auto x97 = x80*x81;
  auto x98 = x96*x97;
  auto x99 = x93*x98;
  auto x100 = 2*x8;
  auto x101 = std::pow(x39, 2);
  auto x102 = x100 + 5*x101 - 4;
  auto x103 = p[31]*x85;
  auto x104 = p[21]*x83;
  auto x105 = x103 - x104;
  auto x106 = x105*x97;
  auto x107 = x102*x106;
  auto x108 = std::pow(f[0], 1.5);
  auto x109 = 2.8868228391288908*x108;
  auto x110 = p[16]*x83;
  auto x111 = p[26]*x85;
  auto x112 = x110 - x111;
  auto x113 = x112*x80;
  auto x114 = 3*x3 + 1;
  auto x115 = x114*x81;
  auto x116 = x113*x115;
  auto x117 = p[18]*x83;
  auto x118 = p[28]*x85;
  auto x119 = x117 - x118;
  auto x120 = x119*x81;
  auto x121 = 3*x8 - 1;
  auto x122 = x121*x80;
  auto x123 = x120*x122;
  auto x124 = p[19]*x83;
  auto x125 = p[29]*x85;
  auto x126 = x124 - x125;
  auto x127 = x121*x126;
  auto x128 = x114*x97;
  auto x129 = x127*x128;
  auto x130 = std::pow(x39, 3.0/2.0);
  auto x131 = x130*x79;
  auto x132 = x120*x131;
  auto x133 = 3.4641874069546694*x27;
  auto x134 = x126*x131;
  auto x135 = x115*x134;
  auto x136 = std::pow(f[0], 3.9462219947249002);
  auto x137 = 3.461720028313473*x136;
  auto x138 = 1.0/f[0];
  auto x139 = std::pow(f[0], 2.0);
  auto x140 = 0.40407102478162499*x139;
  auto x141 = p[7]*x140;
  auto x142 = x20*x40;
  auto x143 = std::pow(f[0], 5.21110255092798);
  auto x144 = 0.45176513957485798*x143;
  auto x145 = x25*x40;
  auto x146 = 0.69987154478819702*std::pow(f[0], 3.2915026221291801);
  auto x147 = x39*x8;
  auto x148 = 1.5649603496641979*std::pow(f[0], 4.3245553203367599);
  auto x149 = x13*x148;
  auto x150 = p[11]*x22;
  auto x151 = 2.7105908374491401*x143;
  auto x152 = x150*x151;
  auto x153 = x56*x8;
  auto x154 = 0.34993577239409901*x71;
  auto x155 = x154*x45;
  auto x156 = std::pow(f[0], 4.5);
  auto x157 = 0.23329051492939901*x156;
  auto x158 = x31*x37;
  auto x159 = x153*x158;
  auto x160 = std::pow(f[0], 5.0);
  auto x161 = 0.45176513957485798*x160;
  auto x162 = x161*x31;
  auto x163 = x8 - 2*std::cos(x47);
  auto x164 = x163*x57;
  auto x165 = x19*x65;
  auto x166 = std::pow(f[0], 3.5);
  auto x167 = 0.98976783055847695*x166;
  auto x168 = x154*x44;
  auto x169 = std::pow(x28, 2);
  auto x170 = -1.0/2.0*x169 - 1.0/2.0;
  auto x171 = x170/x169;
  auto x172 = x171*x30;
  auto x173 = x172*x41;
  auto x174 = 0.23329051492939901*x108;
  auto x175 = x174*x78;
  auto x176 = 1.0/x75;
  auto x177 = x176*x19;
  auto x178 = x177*x88;
  auto x179 = x172*x42;
  auto x180 = std::pow(x28, -3.0/2.0);
  auto x181 = -1.0/4.0*x169 - 1.0/4.0;
  auto x182 = x180*x181*x77;
  auto x183 = 0.46658102985879801*x108;
  auto x184 = x4*x75;
  auto x185 = x184*x88;
  auto x186 = x183*x185;
  auto x187 = x161*x172;
  auto x188 = 10*x147 - 2*x39;
  auto x189 = x106*x188;
  auto x190 = std::pow(f[0], 5.7620873481300103);
  auto x191 = 0.40407102478162499*x190;
  auto x192 = std::pow(f[0], 4.9462219947249002);
  auto x193 = 0.69987154478819702*x192;
  auto x194 = x112*x115;
  auto x195 = 0.16496130509307949*x166;
  auto x196 = x177*x78;
  auto x197 = x195*x196;
  auto x198 = x120*x121;
  auto x199 = 0.32992261018615898*x166;
  auto x200 = x182*x184;
  auto x201 = x199*x200;
  auto x202 = x93*x96;
  auto x203 = 0.20203551239081249*x190;
  auto x204 = x196*x81;
  auto x205 = x203*x204;
  auto x206 = x102*x105;
  auto x207 = x200*x81;
  auto x208 = x191*x202;
  auto x209 = x206*x207;
  auto x210 = 0.1166452574646995*x192;
  auto x211 = x115*x127;
  auto x212 = x196*x211;
  auto x213 = 0.23329051492939901*x192;
  auto x214 = x200*x211;
  auto x215 = -x132*x167 - x135*x193 - x141*x40 - x142*x144 - x145*x146 - x147*x149 - x147*x152 + x153*x155 + x157*x159 + x157*x179 + x162*x164 + x162*x165 + x168*x173 + x175*x178 + x182*x186 + x187*x58 + x187*x66 + x189*x191 + x191*x209 + x194*x197 + x194*x201 + x197*x198 + x198*x201 + x202*x205 + x205*x206 + x207*x208 + x210*x212 + x213*x214;
  auto x216 = std::pow(f[0], 3.3245553203367599);
  auto x217 = 6.7677576062563833*x216;
  auto x218 = 2.3036290248239504*std::pow(f[0], 2.2915026221291801);
  auto x219 = std::pow(f[0], 4.21110255092798);
  auto x220 = 2.3541944712588774*x219;
  auto x221 = 2.2588256978742898*std::pow(f[0], 4.0);
  auto x222 = x221*x31;
  auto x223 = 1.0498073171822955*x166;
  auto x224 = 1.049807317182297*x139;
  auto x225 = x224*x72;
  auto x226 = x172*x221;
  auto x227 = std::pow(f[0], 4.7620873481300103);
  auto x228 = 2.3282925396401293*x227;
  auto x229 = x224*x71;
  auto x230 = x74*x78;
  auto x231 = 0.69987154478819702*x74;
  auto x232 = x185*x231;
  auto x233 = 0.57736456782577816*x27;
  auto x234 = x196*x233;
  auto x235 = 1.1547291356515563*x27;
  auto x236 = x200*x235;
  auto x237 = 1.1641462698200646*x227;
  auto x238 = x204*x237;
  auto x239 = x202*x228;
  auto x240 = 0.57695333805224558*x136;
  auto x241 = 1.1539066761044912*x136;
  auto x242 = -x13*x147*x217 - x138*x215 - x142*x220 - x145*x218 - 14.125166827553224*x147*x150*x219 + x153*x225 + x159*x223 + x164*x222 + x165*x222 + x173*x229 + x177*x230*x89 + x179*x223 + x182*x232 + x189*x228 + x194*x234 + x194*x236 + x198*x234 + x198*x236 + x202*x238 + x206*x238 + x207*x239 + x209*x228 + x212*x240 + x214*x241 + x226*x58 + x226*x66 - x39*x67*x9;
  auto x243 = std::pow(x38, 3.0/2.0);
  auto x244 = x113*x243;
  auto x245 = x243*x80;
  auto x246 = x127*x245;
  auto x247 = p[5]*x140;
  auto x248 = x23*x56;
  auto x249 = 1.0*x143;
  auto x250 = p[8]*x146;
  auto x251 = x3*x38;
  auto x252 = x148*x16;
  auto x253 = p[9]*x19;
  auto x254 = 2.7105908374491481*x143;
  auto x255 = x253*x254;
  auto x256 = std::pow(x30, 2);
  auto x257 = (1.0/2.0)*x256 + 1.0/2.0;
  auto x258 = x257*x29;
  auto x259 = x258*x41;
  auto x260 = x3*x40;
  auto x261 = x258*x42;
  auto x262 = x158*x260;
  auto x263 = x161*x258;
  auto x264 = x162*x49;
  auto x265 = x22*x55;
  auto x266 = x162*x40;
  auto x267 = x3 + 2*std::cos(x62);
  auto x268 = x267*x61;
  auto x269 = 1.0/x81;
  auto x270 = x22*x269;
  auto x271 = x270*x75;
  auto x272 = x271*x87;
  auto x273 = (1.0/4.0)*x256;
  auto x274 = x273 + 1.0/4.0;
  auto x275 = x274/x77;
  auto x276 = x275*x76;
  auto x277 = -10*x251 - 2*x38;
  auto x278 = x277*x98;
  auto x279 = x192*x243;
  auto x280 = x126*x279;
  auto x281 = 0.69987154478819702*x122;
  auto x282 = x271*x78;
  auto x283 = x195*x282;
  auto x284 = x112*x114;
  auto x285 = x119*x121;
  auto x286 = x184*x276;
  auto x287 = x194*x286;
  auto x288 = x121*x199;
  auto x289 = x120*x276;
  auto x290 = x184*x289;
  auto x291 = x203*x282;
  auto x292 = x286*x81;
  auto x293 = x206*x292;
  auto x294 = x210*x282;
  auto x295 = x114*x127;
  auto x296 = x211*x286;
  auto x297 = -x153*x250 + x155*x260 + x157*x261 + x157*x262 - x167*x244 + x168*x259 + x175*x272 + x186*x276 + x191*x278 + x191*x293 + x199*x287 + x202*x291 + x206*x291 + x208*x292 + x213*x296 - x247*x56 - x248*x249 - x251*x252 - x251*x255 + x263*x58 + x263*x66 + x264*x265 + x266*x268 - x280*x281 + x283*x284 + x283*x285 + x288*x290 + x294*x295;
  auto x298 = 5.21110255092798*x219;
  auto x299 = x221*x258;
  auto x300 = x222*x49;
  auto x301 = x222*x40;
  auto x302 = x233*x282;
  auto x303 = x237*x282;
  auto x304 = -p[8]*x153*x218 + x121*x235*x290 - x138*x297 - x16*x217*x251 + x202*x303 + x206*x303 - 14.125166827553265*x219*x251*x253 + x223*x261 + x223*x262 + x224*x3*x73 + x228*x278 + x228*x293 + x229*x259 + 0.34993577239409851*x230*x272 + x232*x276 + x235*x287 + x239*x292 + x240*x282*x295 + x241*x296 - x248*x298 + x265*x300 + x268*x301 + x284*x302 + x285*x302 + x299*x58 + x299*x66 - x6*x68;
  auto x305 = 1.0*x51;
  auto x306 = p[12]*x305;
  auto x307 = 1.0*x53;
  auto x308 = p[22]*x307;
  auto x309 = -x306 - x308;
  auto x310 = x309*x41;
  auto x311 = 0.34993577239409901*x45;
  auto x312 = x31*x41;
  auto x313 = 1.5*p[20]*x35;
  auto x314 = 1.5*p[30]*x33;
  auto x315 = -x313 - x314;
  auto x316 = x157*x315;
  auto x317 = 0.5*x85;
  auto x318 = p[15]*x317;
  auto x319 = 0.5*x83;
  auto x320 = p[25]*x319;
  auto x321 = -x318 - x320;
  auto x322 = x321*x97;
  auto x323 = p[14]*x305;
  auto x324 = p[24]*x307;
  auto x325 = x323 + x324;
  auto x326 = x325*x56;
  auto x327 = p[13]*x305;
  auto x328 = p[23]*x307;
  auto x329 = -x327 - x328;
  auto x330 = x329*x64;
  auto x331 = p[16]*x317;
  auto x332 = p[26]*x319;
  auto x333 = -x331 - x332;
  auto x334 = x128*x333;
  auto x335 = p[18]*x317;
  auto x336 = p[28]*x319;
  auto x337 = -x335 - x336;
  auto x338 = x121*x97;
  auto x339 = x337*x338;
  auto x340 = x191*x97;
  auto x341 = p[21]*x317;
  auto x342 = p[31]*x319;
  auto x343 = x341 + x342;
  auto x344 = x102*x343;
  auto x345 = p[17]*x317;
  auto x346 = p[27]*x319;
  auto x347 = -x345 - x346;
  auto x348 = x347*x93;
  auto x349 = p[19]*x317;
  auto x350 = p[29]*x319;
  auto x351 = -x349 - x350;
  auto x352 = x121*x128;
  auto x353 = x351*x352;
  auto x354 = x183*x322 + x199*x334 + x199*x339 + x213*x353 + x264*x326 + x266*x330 + x310*x311 + x312*x316 + x340*x344 + x340*x348;
  auto x355 = x223*x312;
  auto x356 = x228*x97;
  auto x357 = -x138*x354 + x224*x31*x310 + x231*x322 + x235*x334 + x235*x339 + x241*x353 + x300*x326 + x301*x330 + x315*x355 + x344*x356 + x348*x356;
  auto x358 = x306 + x308;
  auto x359 = x311*x358;
  auto x360 = x313 + x314;
  auto x361 = x157*x360;
  auto x362 = x318 + x320;
  auto x363 = x362*x97;
  auto x364 = -x323 - x324;
  auto x365 = x364*x56;
  auto x366 = x327 + x328;
  auto x367 = x366*x64;
  auto x368 = x331 + x332;
  auto x369 = x128*x368;
  auto x370 = x335 + x336;
  auto x371 = x338*x370;
  auto x372 = x345 + x346;
  auto x373 = x372*x93;
  auto x374 = -x341 - x342;
  auto x375 = x102*x374;
  auto x376 = x349 + x350;
  auto x377 = x352*x376;
  auto x378 = x183*x363 + x199*x369 + x199*x371 + x213*x377 + x264*x365 + x266*x367 + x312*x361 + x340*x373 + x340*x375 + x359*x41;
  auto x379 = -x138*x378 + x224*x312*x358 + x231*x363 + x235*x369 + x235*x371 + x241*x377 + x300*x365 + x301*x367 + x355*x360 + x356*x373 + x356*x375;
  auto x380 = 3.464187406954669*x27;
  auto x381 = 3.4617200283134735*x136;
  auto x382 = 0.69987154478819802*x44*x71;
  auto x383 = x153*x172;
  auto x384 = 0.46658102985879801*x156*x37;
  auto x385 = 0.90353027914971595*x160;
  auto x386 = x172*x385;
  auto x387 = x169 + 1;
  auto x388 = 0.1749678861970495*x41*x45*x71;
  auto x389 = 0.1166452574646995*x108;
  auto x390 = x12*x79/x130;
  auto x391 = x80*x88;
  auto x392 = x389*x391;
  auto x393 = 0.1166452574646995*x156*x43;
  auto x394 = 0.22588256978742899*x160*x31;
  auto x395 = x387*x394;
  auto x396 = x166*x78;
  auto x397 = 1.9795356611169539*x120;
  auto x398 = x19*x75;
  auto x399 = x130*x4;
  auto x400 = x182*x399;
  auto x401 = x168*x41;
  auto x402 = x170*x30*(-x169 - 1)/std::pow(x28, 3);
  auto x403 = x157*x42;
  auto x404 = x181*x186;
  auto x405 = x77*(-3.0/4.0*x169 - 3.0/4.0)/std::pow(x28, 5.0/2.0);
  auto x406 = x161*x402;
  auto x407 = x105*x188;
  auto x408 = x191*x407;
  auto x409 = x75*x81;
  auto x410 = x190*x409*x5;
  auto x411 = 0.082480652546539746*x166;
  auto x412 = x390*x411;
  auto x413 = x387*x411;
  auto x414 = x115*x126;
  auto x415 = x192*x78;
  auto x416 = 1.399743089576394*x415;
  auto x417 = x202*x81;
  auto x418 = 0.10101775619540625*x190;
  auto x419 = x390*x418;
  auto x420 = x206*x81;
  auto x421 = x387*x418;
  auto x422 = x194*x199;
  auto x423 = x181*x184;
  auto x424 = x405*x423;
  auto x425 = x198*x199;
  auto x426 = x177*x182;
  auto x427 = x208*x81;
  auto x428 = x191*x420;
  auto x429 = 0.058322628732349752*x192;
  auto x430 = x129*x429;
  auto x431 = x211*x213;
  auto x432 = 4.3245553203367599*x216;
  auto x433 = f[0]*(f[0]*x1 + x10*x67 + x107*x228 + x116*x235 + x123*x235 + x129*x241 + x14*x432 + x17*x432 + x21*x220 + x218*x26 + x222*x58 + x222*x66 + x223*x43 + x225*x41 + x228*x99 + x231*x391 + x24*x298 + x67*x7 + 0.68391798958577998*x2/std::pow(f[0], 5));
  auto x434 = -x107*x203 - x116*x195 - x123*x195 - x129*x210 - x146*x26 - x155*x41 - x157*x43 - x174*x391 - x203*x99 + (1.0/6.0)*x433;
  auto x435 = x153*x258;
  auto x436 = x19*x3;
  auto x437 = x157*x37;
  auto x438 = x162*x163;
  auto x439 = x162*x19;
  auto x440 = x171*x257;
  auto x441 = x172*x260;
  auto x442 = x161*x440;
  auto x443 = x265*x49;
  auto x444 = x268*x40;
  auto x445 = 0.49488391527923847*x396;
  auto x446 = x112*x243;
  auto x447 = x130*x270;
  auto x448 = x389*x87;
  auto x449 = x176*x269*x436;
  auto x450 = x449*x78;
  auto x451 = x180*x275;
  auto x452 = x277*x96;
  auto x453 = x191*x207;
  auto x454 = x121*x280;
  auto x455 = x114*x126;
  auto x456 = x112*x396;
  auto x457 = 0.082480652546539746*x449;
  auto x458 = x182*x271;
  auto x459 = x195*x458;
  auto x460 = x177*x195;
  auto x461 = x423*x451;
  auto x462 = x418*x450;
  auto x463 = x203*x458;
  auto x464 = x177*x276;
  auto x465 = x203*x464;
  auto x466 = p[11]*x147*x151*x56 + p[9]*x251*x254*x40 + x114*x456*x457 - x119*x445*x447 + x121*x289*x460 + 2.0996146343645909*x134*x279 + x155*x436 + x157*x158*x436 + x164*x263 + x165*x263 - x167*x200*x446 - x167*x289*x399 + x168*x435 + x168*x441 + x174*x178*x276 + x174*x182*x272 - x177*x445*x446 + x187*x443 + x187*x444 - x193*x276*x399*x414 + x194*x276*x460 - 0.34993577239409851*x196*x454 - 0.69987154478819702*x200*x454 + x202*x462 + x202*x463 + x205*x452 + x206*x462 + x206*x463 + x210*x211*x464 + x210*x295*x458 + x250*x41 + x265*x438 + x268*x439 + x284*x459 + x285*x396*x457 + x285*x459 + x291*x407 + x292*x408 + 0.058322628732349752*x295*x415*x449 + x401*x440 + x403*x440 + x404*x451 - 0.34993577239409851*x415*x447*x455 + x417*x465 + x420*x465 + x422*x461 + x425*x461 + x427*x461 + x428*x461 + x431*x461 + x435*x437 + x437*x441 + x442*x58 + x442*x66 + x448*x450 + x452*x453;
  auto x467 = 1.0/x39;
  auto x468 = x354/std::tan(f[1]);
  auto x469 = x309*x311;
  auto x470 = x153*x31;
  auto x471 = x131*x81;
  auto x472 = x167*x471;
  auto x473 = 0.34993577239409901*x44;
  auto x474 = x321*x81;
  auto x475 = x175*x177;
  auto x476 = x183*x200;
  auto x477 = x187*x49;
  auto x478 = x187*x40;
  auto x479 = x188*x340;
  auto x480 = x114*x351;
  auto x481 = x193*x471;
  auto x482 = x115*x333;
  auto x483 = x121*x81;
  auto x484 = x337*x483;
  auto x485 = x115*x121;
  auto x486 = x351*x485;
  auto x487 = x196*x210;
  auto x488 = x200*x213;
  auto x489 = x153*x469 + x172*x310*x473 + x173*x316 + x197*x482 + x197*x484 + x201*x482 + x201*x484 + x205*x344 + x205*x348 + x316*x470 + x326*x438 + x326*x477 + x330*x439 + x330*x478 - x337*x472 + x343*x479 + x344*x453 + x348*x453 + x474*x475 + x474*x476 - x480*x481 + x486*x487 + x486*x488;
  auto x490 = -x354*((1.0/6.0)*x39 - 2.0/3.0*x467) - 2.0/3.0*x468 + x489;
  auto x491 = (1.0/3.0)*x468;
  auto x492 = (1.0/3.0)*x354;
  auto x493 = x467*x492;
  auto x494 = x362*x81;
  auto x495 = x114*x376;
  auto x496 = x115*x368;
  auto x497 = x370*x483;
  auto x498 = x376*x485;
  auto x499 = x153*x359 + x173*x358*x473 + x173*x361 + x197*x496 + x197*x497 + x201*x496 + x201*x497 + x205*x373 + x205*x375 - x3*x491 + x3*x493 + x361*x470 + x365*x438 + x365*x477 + x367*x439 + x367*x478 - x370*x472 + x373*x453 + x374*x479 + x375*x453 + x475*x494 + x476*x494 - x481*x495 + x487*x498 + x488*x498;
  auto x500 = x489 - x491 + x493;
  auto x501 = x258*x260;
  auto x502 = x258*x385;
  auto x503 = x256 + 1;
  auto x504 = x394*x503;
  auto x505 = x15/x243;
  auto x506 = x505*x80;
  auto x507 = x22*x409;
  auto x508 = x411*x503;
  auto x509 = x274*x76*(-x273 - 1.0/4.0)/std::pow(x30, 3.0/2.0);
  auto x510 = x418*x503;
  auto x511 = x411*x505;
  auto x512 = x418*x506;
  auto x513 = x184*x509;
  auto x514 = x271*x276;
  auto x515 = x199*x514;
  auto x516 = 1.0/std::tan(f[2]);
  auto x517 = x492*x516;
  auto x518 = 1.0/x38;
  auto x519 = (1.0/3.0)*x378*x518;
  auto x520 = x259*x473;
  auto x521 = x260*x31;
  auto x522 = x263*x49;
  auto x523 = x263*x40;
  auto x524 = x22*x264;
  auto x525 = x266*x267;
  auto x526 = x167*x245;
  auto x527 = x175*x271;
  auto x528 = x183*x286;
  auto x529 = x277*x340;
  auto x530 = x279*x281;
  auto x531 = x114*x283;
  auto x532 = x121*x283;
  auto x533 = x199*x286;
  auto x534 = x288*x292;
  auto x535 = x191*x292;
  auto x536 = x121*x294;
  auto x537 = x213*x286;
  auto x538 = x259*x316 + x260*x469 + x291*x344 + x291*x348 + x309*x520 + x316*x521 + x321*x527 + x325*x524 + x326*x522 + x329*x525 + x330*x523 - x333*x526 + x333*x531 + x337*x532 + x337*x534 + x344*x535 + x347*x529 + x348*x535 - x351*x530 + x474*x528 + x480*x536 + x482*x533 + x486*x537;
  auto x539 = -x517*x8 + x519*x8 + x538;
  auto x540 = x259*x361 + x260*x359 + x291*x373 + x291*x375 - x354*((1.0/6.0)*x38 - 2.0/3.0*x518) + x358*x520 + x361*x521 + x362*x527 + x364*x524 + x365*x522 + x366*x525 + x367*x523 - x368*x526 + x368*x531 + x370*x532 + x370*x534 + x372*x529 + x373*x535 + x375*x535 - x376*x530 - 2.0/3.0*x378*x516 + x494*x528 + x495*x536 + x496*x533 + x498*x537;
  auto x541 = -x517 + x519 + x538;
  auto x542 = (1.0/18.0)*x433;
  auto x543 = 1.0*x70;
  auto x544 = 1.0*x69;
  auto x545 = x311*x41;
  auto x546 = 2.25*x36;
  auto x547 = 2.25*x34;
  auto x548 = x157*x312;
  auto x549 = 0.25*x86;
  auto x550 = 0.25*x84;
  auto x551 = x183*x97;
  auto x552 = 1.0*x54;
  auto x553 = 1.0*x52;
  auto x554 = x264*x56;
  auto x555 = 1.0*x60;
  auto x556 = 1.0*x59;
  auto x557 = x266*x64;
  auto x558 = 0.25*x111;
  auto x559 = 0.25*x110;
  auto x560 = x128*x199;
  auto x561 = 0.25*x118;
  auto x562 = 0.25*x117;
  auto x563 = x199*x338;
  auto x564 = 0.25*x104;
  auto x565 = 0.25*x103;
  auto x566 = x102*x340;
  auto x567 = 0.25*x95;
  auto x568 = 0.25*x94;
  auto x569 = x340*x93;
  auto x570 = 0.25*x125;
  auto x571 = 0.25*x124;
  auto x572 = x213*x352;
  auto x573 = x545*(x543 - x544) + x548*(x546 - x547) + x551*(x549 - x550) + x554*(x552 - x553) + x557*(x555 - x556) + x560*(x558 - x559) + x563*(x561 - x562) + x566*(x564 - x565) + x569*(x567 - x568) + x572*(x570 - x571);
  auto x574 = (1.0/3.0)*x297*x38;
  auto x575 = (1.0/3.0)*x215*x39;
  auto x576 = (1.0/9.0)*x433;
  auto x577 = x576*x8;
  auto x578 = x545*(-x543 + x544) + x548*(-x546 + x547) + x551*(-x549 + x550) + x554*(-x552 + x553) + x557*(-x555 + x556) + x560*(-x558 + x559) + x563*(-x561 + x562) + x566*(-x564 + x565) + x569*(-x567 + x568) + x572*(-x570 + x571);
  auto x579 = -x3*x575 + x3*x577 - x574*x8 + x578;
  auto x580 = x573 - x575 + x577;
  auto x581 = x3*x576 - x574 + x578;

 sum[0]=5.2787719507969681*std::pow(f[0], 1.2915026221291801)*x26 + 4.5535719051034445*std::pow(f[0], 2.9462219947249002)*x129 + x1 + x10 + x107*x90 + x109*x116 + x109*x123 + x11*x14 + x11*x17 + 9.9137543432988053*x18*x21 + 21.944487245360122*x18*x24 + 3.6743256101380339*x27*x43 + x46*x58 + x46*x66 + 2.0996146343645941*x68*x73 + x7 + x90*x99 + x80*x89/x74 - 3.4195899479289*x2/std::pow(f[0], 6);

 sum[6]=-x132*x380 - x135*x381 + x242;

 sum[12]=-x244*x380 - x246*x381 + x304;

 sum[18]=x357;

 sum[24]=x379;

 sum[30]=x357;

 sum[1]=-x132*x133 - x135*x137 + x242;

 sum[7]=x101*x149 + x101*x152 + x106*x191*(-x100 - 10*x101 + 10*x12) - x107*x421 - x116*x413 - x12*x149 - x12*x152 - x123*x413 - x141*x19 - x144*x21 + x162*x57*(-x39 + 4*x48) - x162*x66 + x164*x386 + x165*x386 - x166*x397*x400 + x178*x182*x183 + x182*x407*x410 - 1.399743089576394*x192*x400*x414 - x194*x412 - x198*x412 + x204*x408 - x211*x390*x429 + x382*x383 + x383*x384 - x387*x388 - x387*x392 - x387*x393 - x387*x430 - x389*x390*x88 - x395*x58 - x395*x66 - x396*x397*x398 - x398*x414*x416 + x401*x402 + x402*x403 + x404*x405 + x406*x58 + x406*x66 - x417*x419 - x419*x420 - x421*x99 + x422*x424 + x422*x426 + x424*x425 + x424*x427 + x424*x428 + x424*x431 + x425*x426 + x426*x427 + x426*x428 + x426*x431 + x434;

 sum[13]=x466;

 sum[19]=x490;

 sum[25]=x499;

 sum[31]=x500;

 sum[2]=-x133*x244 - x137*x246 + x304;

 sum[8]=x466;

 sum[14]=x107*x510 - x113*x114*x511 + x116*x508 - x119*x122*x511 - x122*x429*x455*x505 + x123*x508 - x127*x416*x507 - x15*x252 - x15*x255 - x162*x58 - 1.9795356611169539*x166*x286*x446 + x183*x272*x276 + x186*x509 + x191*x206*x514 + x191*x282*x452 + x191*x98*(-10*x15 + 10*x91 - x92) - x202*x512 - x206*x512 + x208*x514 + x213*x295*x514 - x22*x247 - x24*x249 + x252*x91 + x255*x91 + x266*x61*(-x38 - 4*x63) + x276*x410*x452 + x284*x515 + x285*x515 - 1.399743089576394*x286*x454 + x382*x501 + x384*x501 + x388*x503 + x392*x503 + x393*x503 + x422*x513 + x425*x513 + x427*x513 + x428*x513 + x430*x503 + x431*x513 + x434 + x443*x502 + x444*x502 - x448*x506 - 1.9795356611169539*x456*x507 + x504*x58 + x504*x66 + x510*x99;

 sum[20]=x539;

 sum[26]=x540;

 sum[32]=x541;

 sum[3]=x357;

 sum[9]=x490;

 sum[15]=x539;

 sum[21]=(1.0/6.0)*x215*x48 + x542*(x101 + 2) + x573;

 sum[27]=x579;

 sum[33]=x580;

 sum[4]=x379;

 sum[10]=x499;

 sum[16]=x540;

 sum[22]=x579;

 sum[28]=(1.0/6.0)*x297*x63 + x542*(x91 + 2) + x573;

 sum[34]=x581;

 sum[5]=x357;

 sum[11]=x500;

 sum[17]=x541;

 sum[23]=x580;

 sum[29]=x581;

 sum[35]=x573 + x576;
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  auto x0 = std::pow(p[0], 2)*std::pow(p[1], 8)/std::pow(p[3], 4);
  auto x1 = 33.420651142226113*std::pow(f[0], 1.3245553203367599);
  auto x2 = std::cos(f[1]);
  auto x3 = std::pow(x2, 2);
  auto x4 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x5 = p[2]*x4;
  auto x6 = p[10]*x5;
  auto x7 = x6*(0.78248017483209897*x3 - 0.26082672494403297);
  auto x8 = std::cos(f[2]);
  auto x9 = std::pow(x8, 2);
  auto x10 = p[6]*x5;
  auto x11 = x10*(0.78248017483209897*x9 - 0.26082672494403297);
  auto x12 = std::pow(f[0], 2.21110255092798);
  auto x13 = x5*x8;
  auto x14 = p[11]*(1.35529541872457*x3 - 0.45176513957485798);
  auto x15 = x13*x14;
  auto x16 = x13*x2;
  auto x17 = p[8]*x16;
  auto x18 = p[9]*x5;
  auto x19 = x18*(3*x9 - 1);
  auto x20 = x19*x2;
  auto x21 = std::sin(f[1]);
  auto x22 = std::sin(f[2]);
  auto x23 = x21*x22;
  auto x24 = x23*x5;
  auto x25 = 1.0*f[3] - 1.0*f[4] + 1.0*f[5];
  auto x26 = std::cos(x25);
  auto x27 = p[12]*x26;
  auto x28 = std::sin(x25);
  auto x29 = p[22]*x28;
  auto x30 = x27 - x29;
  auto x31 = std::tan((1.0/2.0)*f[1]);
  auto x32 = 1.0/x31;
  auto x33 = std::tan((1.0/2.0)*f[2]);
  auto x34 = x32*x33;
  auto x35 = x30*x34;
  auto x36 = 2.0996146343645941*x35;
  auto x37 = x24*x36;
  auto x38 = std::pow(f[0], 2.0);
  auto x39 = x34*x38;
  auto x40 = 27.105908374491477*x39;
  auto x41 = 2*f[1];
  auto x42 = std::sin(x41);
  auto x43 = x21 - x42;
  auto x44 = p[24]*x28;
  auto x45 = p[14]*x26;
  auto x46 = x44 - x45;
  auto x47 = x22*x5;
  auto x48 = x46*x47;
  auto x49 = x43*x48;
  auto x50 = 2*f[2];
  auto x51 = std::sin(x50);
  auto x52 = x22 + x51;
  auto x53 = p[13]*x26;
  auto x54 = p[23]*x28;
  auto x55 = x5*(x53 - x54);
  auto x56 = x52*x55;
  auto x57 = x21*x56;
  auto x58 = std::pow(f[0], 1.5);
  auto x59 = 1.5*f[3] - 1.5*f[4] + 1.5*f[5];
  auto x60 = std::cos(x59);
  auto x61 = p[20]*x60;
  auto x62 = std::sin(x59);
  auto x63 = p[30]*x62;
  auto x64 = x61 - x63;
  auto x65 = x34*x64;
  auto x66 = x24*x65;
  auto x67 = std::sqrt(x31);
  auto x68 = 1.0/x67;
  auto x69 = std::sqrt(x33);
  auto x70 = x68*x69;
  auto x71 = std::sqrt(x22);
  auto x72 = x5*x71;
  auto x73 = x70*x72;
  auto x74 = std::sqrt(x21);
  auto x75 = 0.5*f[3] - 0.5*f[4] + 0.5*f[5];
  auto x76 = std::cos(x75);
  auto x77 = p[15]*x76;
  auto x78 = std::sin(x75);
  auto x79 = p[25]*x78;
  auto x80 = x77 - x79;
  auto x81 = x74*x80;
  auto x82 = 0.17496788619704925*x81;
  auto x83 = x73*x82;
  auto x84 = std::sqrt(f[0]);
  auto x85 = 4.330234258693336*x84;
  auto x86 = p[16]*x76;
  auto x87 = p[26]*x78;
  auto x88 = x86 - x87;
  auto x89 = x74*x88;
  auto x90 = 3*x8 + 1;
  auto x91 = x73*x90;
  auto x92 = x89*x91;
  auto x93 = p[18]*x76;
  auto x94 = p[28]*x78;
  auto x95 = x93 - x94;
  auto x96 = 3*x2 - 1;
  auto x97 = x73*x74;
  auto x98 = x96*x97;
  auto x99 = x95*x98;
  auto x100 = 41.71226553619632*std::pow(f[0], 2.7620873481300103);
  auto x101 = std::pow(x22, 2);
  auto x102 = 2*x8;
  auto x103 = -5*x101 + x102 + 4;
  auto x104 = p[17]*x76;
  auto x105 = p[27]*x78;
  auto x106 = x104 - x105;
  auto x107 = x106*x97;
  auto x108 = x103*x107;
  auto x109 = 2*x2;
  auto x110 = std::pow(x21, 2);
  auto x111 = x109 + 5*x110 - 4;
  auto x112 = p[31]*x78;
  auto x113 = p[21]*x76;
  auto x114 = x112 - x113;
  auto x115 = x114*x97;
  auto x116 = x111*x115;
  auto x117 = p[19]*x76;
  auto x118 = p[29]*x78;
  auto x119 = x117 - x118;
  auto x120 = x90*x98;
  auto x121 = 1.0/f[0];
  auto x122 = p[7]*x5;
  auto x123 = 0.40407102478162499*x38;
  auto x124 = x122*x123;
  auto x125 = x124*x21;
  auto x126 = x19*x21;
  auto x127 = std::pow(f[0], 5.21110255092798);
  auto x128 = 0.45176513957485798*x127;
  auto x129 = x126*x128;
  auto x130 = std::pow(f[0], 3.2915026221291801);
  auto x131 = 0.69987154478819702*x130;
  auto x132 = x21*x8;
  auto x133 = p[8]*x5;
  auto x134 = x132*x133;
  auto x135 = x131*x134;
  auto x136 = std::pow(f[0], 4.3245553203367599);
  auto x137 = 1.5649603496641979*x136;
  auto x138 = x2*x21;
  auto x139 = x138*x6;
  auto x140 = 2.7105908374491401*x127;
  auto x141 = p[11]*x13;
  auto x142 = x138*x141;
  auto x143 = x140*x142;
  auto x144 = std::pow(f[0], 3.0);
  auto x145 = 0.34993577239409901*x144;
  auto x146 = x2*x22;
  auto x147 = x146*x5;
  auto x148 = x147*x35;
  auto x149 = x145*x148;
  auto x150 = x147*x65;
  auto x151 = std::pow(f[0], 4.5);
  auto x152 = 0.23329051492939901*x151;
  auto x153 = x150*x152;
  auto x154 = std::pow(f[0], 5.0);
  auto x155 = 0.45176513957485798*x154;
  auto x156 = std::cos(x41);
  auto x157 = 2*x156;
  auto x158 = -x157 + x2;
  auto x159 = x34*x48;
  auto x160 = x158*x159;
  auto x161 = x155*x160;
  auto x162 = x155*x34;
  auto x163 = x2*x56;
  auto x164 = x162*x163;
  auto x165 = std::pow(f[0], 3.5);
  auto x166 = 0.98976783055847695*x165;
  auto x167 = std::pow(x21, 3.0/2.0);
  auto x168 = x167*x73;
  auto x169 = x168*x95;
  auto x170 = std::pow(x31, 2);
  auto x171 = 1.0/x170;
  auto x172 = -1.0/2.0*x170 - 1.0/2.0;
  auto x173 = x171*x172;
  auto x174 = x173*x33;
  auto x175 = x145*x24;
  auto x176 = x175*x30;
  auto x177 = x174*x176;
  auto x178 = x73*x80;
  auto x179 = 0.23329051492939901*x58;
  auto x180 = 1.0/x74;
  auto x181 = x180*x2;
  auto x182 = x179*x181;
  auto x183 = x152*x24;
  auto x184 = x174*x64;
  auto x185 = x183*x184;
  auto x186 = 0.46658102985879801*x58;
  auto x187 = x72*x81;
  auto x188 = x186*x187;
  auto x189 = std::pow(x31, -3.0/2.0);
  auto x190 = -1.0/4.0*x170 - 1.0/4.0;
  auto x191 = x190*x69;
  auto x192 = x189*x191;
  auto x193 = x188*x192;
  auto x194 = x155*x174;
  auto x195 = x194*x49;
  auto x196 = 2*x21;
  auto x197 = 10*x138 - x196;
  auto x198 = x115*x197;
  auto x199 = std::pow(f[0], 5.7620873481300103);
  auto x200 = 0.40407102478162499*x199;
  auto x201 = x168*x90;
  auto x202 = std::pow(f[0], 4.9462219947249002);
  auto x203 = x119*x202;
  auto x204 = 0.69987154478819702*x203;
  auto x205 = 0.16496130509307949*x165;
  auto x206 = x181*x88;
  auto x207 = x206*x91;
  auto x208 = x181*x73;
  auto x209 = x205*x95;
  auto x210 = x209*x96;
  auto x211 = 0.32992261018615898*x165;
  auto x212 = x211*x72;
  auto x213 = x192*x89;
  auto x214 = x213*x90;
  auto x215 = x212*x214;
  auto x216 = x192*x72;
  auto x217 = x216*x74;
  auto x218 = x211*x95;
  auto x219 = x218*x96;
  auto x220 = x217*x219;
  auto x221 = 0.20203551239081249*x199;
  auto x222 = x208*x221;
  auto x223 = x103*x106;
  auto x224 = x111*x114;
  auto x225 = x200*x217;
  auto x226 = x223*x225;
  auto x227 = x224*x225;
  auto x228 = 0.1166452574646995*x203;
  auto x229 = x181*x91;
  auto x230 = x229*x96;
  auto x231 = 0.23329051492939901*x203;
  auto x232 = x217*x96;
  auto x233 = x232*x90;
  auto x234 = x231*x233;
  auto x235 = -x125 - x129 - x135 - x137*x139 - x143 + x149 + x153 + x161 + x164 - x166*x169 + x177 + x178*x182 + x185 + x193 + x194*x57 + x195 + x198*x200 - x201*x204 + x205*x207 + x208*x210 + x215 + x220 + x222*x223 + x222*x224 + x226 + x227 + x228*x230 + x234;
  auto x236 = -x121*x235;
  auto x237 = std::pow(f[0], 2.5);
  auto x238 = 3.4641874069546694*x237;
  auto x239 = std::pow(f[0], 3.9462219947249002);
  auto x240 = x119*x239;
  auto x241 = 3.461720028313473*x240;
  auto x242 = std::pow(f[0], 1.0);
  auto x243 = 0.80814204956324998*x5;
  auto x244 = p[7]*x243;
  auto x245 = x21*x244;
  auto x246 = std::pow(f[0], 3.3245553203367599);
  auto x247 = 6.7677576062563833*x246;
  auto x248 = std::pow(f[0], 2.2915026221291801);
  auto x249 = 2.3036290248239504*x248;
  auto x250 = std::pow(f[0], 4.21110255092798);
  auto x251 = 2.3541944712588774*x250;
  auto x252 = 14.125166827553224*x250;
  auto x253 = std::pow(f[0], 4.0);
  auto x254 = 2.2588256978742898*x253;
  auto x255 = x163*x34;
  auto x256 = 1.0498073171822955*x165;
  auto x257 = 1.049807317182297*x38;
  auto x258 = x174*x254;
  auto x259 = std::pow(f[0], 4.7620873481300103);
  auto x260 = 2.3282925396401293*x259;
  auto x261 = x24*x256;
  auto x262 = x24*x30;
  auto x263 = x257*x262;
  auto x264 = x178*x181;
  auto x265 = 0.34993577239409851*x84;
  auto x266 = 0.69987154478819702*x84;
  auto x267 = x187*x192;
  auto x268 = 0.57736456782577816*x237;
  auto x269 = x95*x96;
  auto x270 = x208*x269;
  auto x271 = 1.1547291356515563*x237;
  auto x272 = x214*x72;
  auto x273 = x217*x269;
  auto x274 = 1.1641462698200646*x259;
  auto x275 = x208*x274;
  auto x276 = x217*x260;
  auto x277 = 0.57695333805224558*x240;
  auto x278 = 1.1539066761044912*x240;
  auto x279 = -x126*x251 - x134*x249 - x139*x247 - x142*x252 + x148*x257 + x150*x256 + x160*x254 + x174*x263 + x184*x261 + x198*x260 + x207*x268 + x223*x275 + x223*x276 + x224*x275 + x224*x276 + x230*x277 + x233*x278 - x242*x245 + x254*x255 + x258*x49 + x258*x57 + x264*x265 + x266*x267 + x268*x270 + x271*x272 + x271*x273;
  auto x280 = -x169*x238 - x201*x241 + x279;
  auto x281 = x236 + x280;
  auto x282 = x121*x281;
  auto x283 = 8.6604685173866738*x58;
  auto x284 = std::pow(f[0], 2.9462219947249002);
  auto x285 = x119*x284;
  auto x286 = 13.660715715310332*x285;
  auto x287 = std::pow(f[0], 3.21110255092798);
  auto x288 = 9.9137543432988053*x287;
  auto x289 = 5.2787719507969681*std::pow(f[0], 1.2915026221291801);
  auto x290 = std::pow(f[0], 2.3245553203367599);
  auto x291 = 22.499784556629233*x290;
  auto x292 = 3.6743256101380339*x237;
  auto x293 = 9.0353027914971591*x144;
  auto x294 = x242*x36;
  auto x295 = x24*x64;
  auto x296 = x174*x295;
  auto x297 = x174*x293;
  auto x298 = 2.0996146343645941*x242;
  auto x299 = x24*x298;
  auto x300 = x174*x30;
  auto x301 = std::pow(f[0], 3.7620873481300103);
  auto x302 = 11.087532445765751*x301;
  auto x303 = 1.0/x84;
  auto x304 = 0.17496788619704925*x264;
  auto x305 = 0.34993577239409851*x303;
  auto x306 = 5.5437662228828755*x301;
  auto x307 = x208*x223;
  auto x308 = x224*x306;
  auto x309 = x217*x302;
  auto x310 = 1.4434114195644454*x58;
  auto x311 = 2.8868228391288908*x58;
  auto x312 = 2.2767859525517222*x285;
  auto x313 = 4.5535719051034445*x285;
  auto x314 = -x126*x288 - x134*x289 - x139*x291 - 59.482526059792669*x142*x287 + x147*x294 + x150*x292 + x160*x293 + x198*x302 + x207*x310 + x208*x308 + x223*x309 + x224*x309 + x230*x312 + x233*x313 - x245 + x255*x293 + x267*x305 + x270*x310 + x272*x311 + x273*x311 + x292*x296 + x297*x49 + x297*x57 + x299*x300 + x303*x304 + x306*x307;
  auto x315 = p[5]*x5;
  auto x316 = x123*x315;
  auto x317 = x22*x316;
  auto x318 = 1.0*x127;
  auto x319 = x14*x47;
  auto x320 = x318*x319;
  auto x321 = x133*x146;
  auto x322 = x131*x321;
  auto x323 = x22*x8;
  auto x324 = x10*x137;
  auto x325 = 2.7105908374491481*x127;
  auto x326 = x18*x2;
  auto x327 = x323*x326;
  auto x328 = x325*x327;
  auto x329 = std::pow(x33, 2);
  auto x330 = (1.0/2.0)*x329 + 1.0/2.0;
  auto x331 = x32*x330;
  auto x332 = x30*x331;
  auto x333 = x175*x332;
  auto x334 = x132*x5;
  auto x335 = x145*x35;
  auto x336 = x334*x335;
  auto x337 = x331*x64;
  auto x338 = x183*x337;
  auto x339 = x334*x65;
  auto x340 = x152*x339;
  auto x341 = x155*x331;
  auto x342 = x341*x57;
  auto x343 = x162*x43;
  auto x344 = x13*x46;
  auto x345 = x343*x344;
  auto x346 = x162*x21;
  auto x347 = std::cos(x50);
  auto x348 = 2*x347;
  auto x349 = x348 + x8;
  auto x350 = x349*x55;
  auto x351 = x346*x350;
  auto x352 = std::pow(x22, 3.0/2.0);
  auto x353 = x5*x70;
  auto x354 = x352*x353;
  auto x355 = x354*x89;
  auto x356 = 1.0/x71;
  auto x357 = x179*x356;
  auto x358 = x13*x70;
  auto x359 = x358*x81;
  auto x360 = 1.0/x69;
  auto x361 = (1.0/4.0)*x329;
  auto x362 = x361 + 1.0/4.0;
  auto x363 = x360*x362;
  auto x364 = x363*x68;
  auto x365 = 2*x22;
  auto x366 = -10*x323 - x365;
  auto x367 = x107*x366;
  auto x368 = x354*x74;
  auto x369 = x368*x96;
  auto x370 = x89*x90;
  auto x371 = x356*x358;
  auto x372 = x205*x371;
  auto x373 = x371*x74;
  auto x374 = x211*x89;
  auto x375 = x364*x72;
  auto x376 = x375*x90;
  auto x377 = x364*x95;
  auto x378 = x74*x96;
  auto x379 = x377*x378;
  auto x380 = x221*x373;
  auto x381 = x375*x74;
  auto x382 = x200*x223;
  auto x383 = x224*x381;
  auto x384 = x373*x96;
  auto x385 = x384*x90;
  auto x386 = x231*x378;
  auto x387 = -x166*x355 + x188*x364 + x200*x367 + x200*x383 - x204*x369 + x210*x373 + x212*x379 + x223*x380 + x224*x380 + x228*x385 - x317 - x320 - x322 - x323*x324 - x328 + x333 + x336 + x338 + x340 + x341*x49 + x342 + x345 + x351 + x357*x359 + x370*x372 + x374*x376 + x376*x386 + x381*x382;
  auto x388 = -x121*x387;
  auto x389 = p[5]*x243;
  auto x390 = x22*x389;
  auto x391 = x10*x247;
  auto x392 = 5.21110255092798*x250;
  auto x393 = 14.125166827553265*x250;
  auto x394 = x254*x331;
  auto x395 = x254*x34;
  auto x396 = x344*x43;
  auto x397 = x21*x350;
  auto x398 = x24*x337;
  auto x399 = x24*x332;
  auto x400 = x257*x35;
  auto x401 = x265*x356;
  auto x402 = x187*x364;
  auto x403 = x370*x371;
  auto x404 = x269*x373;
  auto x405 = x370*x375;
  auto x406 = x379*x72;
  auto x407 = x274*x373;
  auto x408 = x223*x381;
  auto x409 = x376*x378;
  auto x410 = x223*x407 + x224*x407 - x242*x390 - x249*x321 + x256*x339 + x256*x398 + x257*x399 + x260*x367 + x260*x383 + x260*x408 + x266*x402 + x268*x403 + x268*x404 + x271*x405 + x271*x406 + x277*x385 + x278*x409 - x319*x392 - x323*x391 - x327*x393 + x334*x400 + x359*x401 + x394*x49 + x394*x57 + x395*x396 + x395*x397;
  auto x411 = -x238*x355 - x241*x369 + x410;
  auto x412 = x388 + x411;
  auto x413 = x121*x412;
  auto x414 = x287*x327;
  auto x415 = 21.944487245360122*x287;
  auto x416 = x10*x323;
  auto x417 = x293*x331;
  auto x418 = x293*x34;
  auto x419 = x371*x82;
  auto x420 = x223*x373;
  auto x421 = -x289*x321 - x291*x416 + x292*x339 + x292*x398 + x294*x334 + x298*x399 + x302*x367 + x302*x383 + x302*x408 + x303*x419 + x305*x402 + x306*x420 + x308*x373 + x310*x403 + x310*x404 + x311*x405 + x311*x406 + x312*x385 + x313*x409 - x319*x415 - x390 + x396*x418 + x397*x418 + x417*x49 + x417*x57;
  auto x422 = -59.482526059792839*x414 + x421;
  auto x423 = p[20]*x62;
  auto x424 = 1.5*x423;
  auto x425 = p[30]*x60;
  auto x426 = 1.5*x425;
  auto x427 = -x424 - x426;
  auto x428 = x34*x427;
  auto x429 = x183*x428;
  auto x430 = 0.5*x78;
  auto x431 = p[15]*x430;
  auto x432 = 0.5*x76;
  auto x433 = p[25]*x432;
  auto x434 = -x431 - x433;
  auto x435 = x434*x97;
  auto x436 = p[16]*x430;
  auto x437 = p[26]*x432;
  auto x438 = -x436 - x437;
  auto x439 = x90*x97;
  auto x440 = x438*x439;
  auto x441 = p[18]*x430;
  auto x442 = p[28]*x432;
  auto x443 = -x441 - x442;
  auto x444 = x443*x98;
  auto x445 = x200*x97;
  auto x446 = p[21]*x430;
  auto x447 = p[31]*x432;
  auto x448 = x446 + x447;
  auto x449 = x111*x448;
  auto x450 = p[17]*x430;
  auto x451 = p[27]*x432;
  auto x452 = -x450 - x451;
  auto x453 = x103*x452;
  auto x454 = p[19]*x430;
  auto x455 = p[29]*x432;
  auto x456 = -x454 - x455;
  auto x457 = x202*x456;
  auto x458 = 0.23329051492939901*x120;
  auto x459 = 1.0*x28;
  auto x460 = p[12]*x459;
  auto x461 = 1.0*x26;
  auto x462 = p[22]*x461;
  auto x463 = -x460 - x462;
  auto x464 = x34*x463;
  auto x465 = x175*x464;
  auto x466 = p[14]*x459;
  auto x467 = p[24]*x461;
  auto x468 = x466 + x467;
  auto x469 = x468*x47;
  auto x470 = x343*x469;
  auto x471 = p[13]*x459;
  auto x472 = p[23]*x461;
  auto x473 = -x471 - x472;
  auto x474 = x5*x52;
  auto x475 = x473*x474;
  auto x476 = x346*x475;
  auto x477 = x465 + x470 + x476;
  auto x478 = x186*x435 + x211*x440 + x211*x444 + x429 + x445*x449 + x445*x453 + x457*x458 + x477;
  auto x479 = x43*x469;
  auto x480 = x21*x475;
  auto x481 = x24*x428;
  auto x482 = 1.049807317182297*x39;
  auto x483 = x24*x463;
  auto x484 = x260*x97;
  auto x485 = x239*x456;
  auto x486 = 1.1539066761044912*x120;
  auto x487 = x256*x481 + x266*x435 + x271*x440 + x271*x444 + x395*x479 + x395*x480 + x449*x484 + x453*x484 + x482*x483 + x485*x486;
  auto x488 = -x121*x478 + x487;
  auto x489 = x121*x488;
  auto x490 = x302*x97;
  auto x491 = 4.5535719051034445*x120;
  auto x492 = x284*x491;
  auto x493 = x292*x481 + x299*x464 + x305*x435 + x311*x440 + x311*x444 + x418*x479 + x418*x480 + x449*x490 + x453*x490 + x456*x492;
  auto x494 = -2*x489 + x493;
  auto x495 = x424 + x426;
  auto x496 = x34*x495;
  auto x497 = x183*x496;
  auto x498 = x431 + x433;
  auto x499 = x498*x97;
  auto x500 = x436 + x437;
  auto x501 = x439*x500;
  auto x502 = x441 + x442;
  auto x503 = x502*x98;
  auto x504 = x450 + x451;
  auto x505 = x103*x504;
  auto x506 = -x446 - x447;
  auto x507 = x111*x506;
  auto x508 = x454 + x455;
  auto x509 = x202*x508;
  auto x510 = x460 + x462;
  auto x511 = x34*x510;
  auto x512 = x175*x511;
  auto x513 = -x466 - x467;
  auto x514 = x47*x513;
  auto x515 = x343*x514;
  auto x516 = x471 + x472;
  auto x517 = x474*x516;
  auto x518 = x346*x517;
  auto x519 = x512 + x515 + x518;
  auto x520 = x186*x499 + x211*x501 + x211*x503 + x445*x505 + x445*x507 + x458*x509 + x497 + x519;
  auto x521 = x43*x514;
  auto x522 = x21*x517;
  auto x523 = x24*x496;
  auto x524 = x482*x510;
  auto x525 = x239*x508;
  auto x526 = x24*x524 + x256*x523 + x266*x499 + x271*x501 + x271*x503 + x395*x521 + x395*x522 + x484*x505 + x484*x507 + x486*x525;
  auto x527 = -x121*x520 + x526;
  auto x528 = x121*x527;
  auto x529 = x292*x523 + x299*x511 + x305*x499 + x311*x501 + x311*x503 + x418*x521 + x418*x522 + x490*x505 + x490*x507 + x492*x508;
  auto x530 = std::pow(f[0], -2);
  auto x531 = 8.660468517386672*x58;
  auto x532 = 13.660715715310333*x285;
  auto x533 = -x121*x280 - x169*x531 - x201*x532 + x235*x530 - x282 + x314;
  auto x534 = x2*x97;
  auto x535 = x237*x95;
  auto x536 = 6.9283748139093388*x535;
  auto x537 = x167*x72;
  auto x538 = x192*x537;
  auto x539 = x2*x439;
  auto x540 = 6.9234400566269461*x240;
  auto x541 = x110*x6;
  auto x542 = x3*x6;
  auto x543 = x110*x141;
  auto x544 = x141*x3;
  auto x545 = -x21 + 4*x42;
  auto x546 = x159*x545;
  auto x547 = 0.69987154478819802*x144;
  auto x548 = x147*x300;
  auto x549 = 0.46658102985879801*x151;
  auto x550 = x147*x184;
  auto x551 = x158*x48;
  auto x552 = 0.90353027914971595*x154;
  auto x553 = x174*x552;
  auto x554 = -x109 - 10*x110 + 10*x3;
  auto x555 = x115*x554;
  auto x556 = x170 + 1;
  auto x557 = x24*x35;
  auto x558 = 0.1749678861970495*x144;
  auto x559 = x557*x558;
  auto x560 = 0.1166452574646995*x58;
  auto x561 = x3/x167;
  auto x562 = x178*x561;
  auto x563 = x73*x81;
  auto x564 = x563*x58;
  auto x565 = 0.1166452574646995*x564;
  auto x566 = 0.1166452574646995*x151;
  auto x567 = x566*x66;
  auto x568 = 0.22588256978742899*x556;
  auto x569 = x154*x49;
  auto x570 = x34*x569;
  auto x571 = x154*x57;
  auto x572 = x34*x571;
  auto x573 = x165*x95;
  auto x574 = 1.9795356611169539*x573;
  auto x575 = -x170 - 1;
  auto x576 = x172*x575/std::pow(x31, 3);
  auto x577 = x33*x576;
  auto x578 = x577*x64;
  auto x579 = std::pow(x31, -5.0/2.0);
  auto x580 = -3.0/4.0*x170 - 3.0/4.0;
  auto x581 = x191*x580;
  auto x582 = x579*x581;
  auto x583 = x186*x216;
  auto x584 = x181*x80;
  auto x585 = x155*x577;
  auto x586 = x114*x197;
  auto x587 = x200*x208;
  auto x588 = x199*x74;
  auto x589 = x586*x588;
  auto x590 = x192*x71;
  auto x591 = x243*x590;
  auto x592 = 0.082480652546539746*x165;
  auto x593 = x561*x88;
  auto x594 = x593*x91;
  auto x595 = x561*x96;
  auto x596 = x573*x73;
  auto x597 = 0.082480652546539746*x596;
  auto x598 = x556*x592;
  auto x599 = 1.399743089576394*x203;
  auto x600 = x538*x90;
  auto x601 = 0.10101775619540625*x199;
  auto x602 = x561*x73;
  auto x603 = x223*x602;
  auto x604 = x224*x602;
  auto x605 = x556*x601;
  auto x606 = x582*x72;
  auto x607 = x374*x90;
  auto x608 = x606*x74;
  auto x609 = x216*x90;
  auto x610 = x211*x609;
  auto x611 = x181*x216;
  auto x612 = x200*x224;
  auto x613 = x595*x91;
  auto x614 = 0.058322628732349752*x203;
  auto x615 = x120*x203;
  auto x616 = 0.058322628732349752*x556;
  auto x617 = x386*x90;
  auto x618 = x181*x96;
  auto x619 = x609*x618;
  auto x620 = 0.1166452574646995*x120;
  auto x621 = -x108*x221 - x116*x221 - x131*x17 - x152*x66 - x175*x35 - x179*x563 - x203*x620 - x205*x92 - x205*x99;
  auto x622 = -x108*x605 - x116*x605 - x124*x2 - x128*x20 + x137*x541 - x137*x542 + x140*x543 - x140*x544 + x155*x546 - x162*x57 + x163*x553 + x176*x577 + x183*x578 + x188*x582 + x200*x555 + x206*x610 + x219*x608 + x219*x611 + x231*x619 + x382*x608 + x382*x611 + x49*x585 - x534*x574 - x538*x574 - x539*x599 + x547*x548 + x549*x550 + x551*x553 - x556*x559 - x556*x565 - x556*x567 - x560*x562 - x568*x570 - x568*x572 + x57*x585 + x583*x584 + x586*x587 + x589*x591 - x592*x594 - x595*x597 - x598*x92 - x598*x99 - x599*x600 - x601*x603 - x601*x604 + x606*x607 + x606*x617 + x608*x612 + x611*x612 - x613*x614 - x615*x616 + x621;
  auto x623 = 0.66666666666666596*x4;
  auto x624 = x0/std::pow(f[0], 5);
  auto x625 = 4.3245553203367599*x246;
  auto x626 = x389*x8;
  auto x627 = x242*x626;
  auto x628 = x2*x244;
  auto x629 = x242*x628;
  auto x630 = x17*x249;
  auto x631 = x20*x251;
  auto x632 = x15*x392;
  auto x633 = x395*x49;
  auto x634 = x395*x57;
  auto x635 = x256*x66;
  auto x636 = x24*x400;
  auto x637 = f[0]*x623 + x108*x260 + x11*x625 + x116*x260 + x240*x486 + x266*x563 + x271*x92 + x271*x99 + 0.68391798958577998*x624 + x625*x7 + x627 + x629 + x630 + x631 + x632 + x633 + x634 + x635 + x636;
  auto x638 = (1.0/6.0)*f[0];
  auto x639 = x637*x638;
  auto x640 = x622 + x639;
  auto x641 = x121*x640;
  auto x642 = 14.377223398316216*x290;
  auto x643 = x108*x302 + x11*x642 + x116*x302 + x15*x415 + x17*x289 + x20*x288 + x242*x37 + x285*x491 + x292*x66 + x305*x563 + x311*x92 + x311*x99 + x418*x49 + x418*x57 + x623 + x626 + x628 + x642*x7 - 3.4195899479289*x0/std::pow(f[0], 6);
  auto x644 = x638*x643;
  auto x645 = x563*x84;
  auto x646 = x120*x240;
  auto x647 = -x108*x274 - x116*x274 - x268*x92 - x268*x99 - x630 - x635 - x636 + x644 - 0.34993577239409851*x645 - 0.57695333805224558*x646;
  auto x648 = 4.5176513957485795*x253;
  auto x649 = x174*x648;
  auto x650 = 2.0996146343645909*x165;
  auto x651 = 2.0996146343645941*x38;
  auto x652 = x253*x34;
  auto x653 = 1.1294128489371449*x652;
  auto x654 = x556*x653;
  auto x655 = x165*x66;
  auto x656 = 0.52490365859114774*x655;
  auto x657 = x38*x557;
  auto x658 = 0.52490365859114851*x657;
  auto x659 = 0.17496788619704925*x84;
  auto x660 = x83*x84;
  auto x661 = x254*x577;
  auto x662 = x260*x586;
  auto x663 = 4.6565850792802586*x259;
  auto x664 = x217*x586;
  auto x665 = x187*x266;
  auto x666 = x216*x266;
  auto x667 = 0.28868228391288908*x237;
  auto x668 = 0.28868228391288908*x535;
  auto x669 = x595*x73;
  auto x670 = x556*x667;
  auto x671 = 0.58207313491003232*x259;
  auto x672 = x556*x671;
  auto x673 = x271*x370;
  auto x674 = x269*x271;
  auto x675 = x271*x609;
  auto x676 = x260*x608;
  auto x677 = x260*x611;
  auto x678 = 0.28847666902612279*x240;
  auto x679 = 0.28847666902612279*x646;
  auto x680 = x378*x90;
  auto x681 = x278*x680;
  auto x682 = -x108*x672 - x116*x672 + x163*x649 + x206*x675 + x208*x662 + x223*x676 + x223*x677 + x224*x676 + x224*x677 + x247*x541 - x247*x542 + x252*x543 - x252*x544 + x254*x546 + x260*x555 + x261*x578 + x263*x577 + x278*x619 - x49*x654 + x49*x661 + x548*x651 + x550*x650 + x551*x649 - x556*x656 - x556*x658 - x556*x660 - x556*x679 - x562*x659 - x57*x654 + x57*x661 + x582*x665 + x584*x666 - x594*x667 - x603*x671 - x604*x671 + x606*x673 + x606*x681 + x608*x674 + x611*x674 - x613*x678 + x663*x664 - x668*x669 - x670*x92 - x670*x99;
  auto x683 = -x121*x622 - x540*x600 - x629 - x631 - x634 - x641 + x647 + x682;
  auto x684 = x167*x354;
  auto x685 = x240*x684;
  auto x686 = x206*x354;
  auto x687 = x237*x686;
  auto x688 = x352*x5;
  auto x689 = x213*x688;
  auto x690 = x181*x354;
  auto x691 = x690*x96;
  auto x692 = 1.7308600141567365*x240;
  auto x693 = x192*x688;
  auto x694 = x378*x693;
  auto x695 = x133*x23;
  auto x696 = x18*x21;
  auto x697 = x325*x696;
  auto x698 = p[11]*x47;
  auto x699 = x140*x698;
  auto x700 = x147*x332;
  auto x701 = x147*x337;
  auto x702 = x16*x65;
  auto x703 = 2.0996146343645909*x203;
  auto x704 = x158*x344;
  auto x705 = x2*x350;
  auto x706 = x173*x330;
  auto x707 = x145*x334;
  auto x708 = x64*x706;
  auto x709 = x152*x334;
  auto x710 = x155*x706;
  auto x711 = 0.49488391527923847*x165;
  auto x712 = x167*x371;
  auto x713 = 0.49488391527923847*x573;
  auto x714 = x166*x688;
  auto x715 = x377*x537;
  auto x716 = x16*x356*x70;
  auto x717 = x180*x80;
  auto x718 = x716*x717;
  auto x719 = x13*x81;
  auto x720 = x192*x357;
  auto x721 = x375*x80;
  auto x722 = x190*x363;
  auto x723 = x189*x722;
  auto x724 = x106*x366;
  auto x725 = x200*x381;
  auto x726 = 0.34993577239409851*x203;
  auto x727 = x712*x90;
  auto x728 = x364*x537;
  auto x729 = x728*x90;
  auto x730 = x180*x88;
  auto x731 = x716*x730*x90;
  auto x732 = 0.082480652546539746*x573;
  auto x733 = x180*x716;
  auto x734 = x733*x96;
  auto x735 = x13*x356;
  auto x736 = x214*x735;
  auto x737 = x192*x735;
  auto x738 = x737*x74;
  auto x739 = x206*x376;
  auto x740 = x205*x72;
  auto x741 = x377*x618;
  auto x742 = x72*x723;
  auto x743 = x74*x742;
  auto x744 = x223*x733;
  auto x745 = x224*x601;
  auto x746 = x221*x738;
  auto x747 = x181*x375;
  auto x748 = x221*x223;
  auto x749 = x224*x747;
  auto x750 = x734*x90;
  auto x751 = x680*x737;
  auto x752 = x376*x618;
  auto x753 = x131*x695 + x138*x699 + x145*x700 + x152*x701 + x152*x702 + x16*x335 + x162*x704 + x162*x705 + x163*x341 - x166*x715 + x176*x706 + x182*x721 + x183*x708 + x184*x709 + x188*x723 + x194*x396 + x194*x397 - x204*x694 - x204*x729 + x205*x736 + x205*x739 + x210*x738 - x213*x714 + x219*x743 + x221*x749 + x222*x724 + x223*x746 + x224*x746 + x225*x724 + x228*x751 + x228*x752 + x300*x707 + x323*x697 + x341*x551 + x380*x586 + x382*x743 + x49*x710 + x560*x718 + x57*x710 + x586*x725 + x592*x731 + x601*x744 + x607*x742 + x612*x743 + x614*x750 + x617*x742 + x684*x703 - x686*x711 - x691*x726 - x712*x713 + x719*x720 - x726*x727 + x732*x734 + x733*x745 + x740*x741 + x747*x748;
  auto x754 = x323*x696;
  auto x755 = x138*x698;
  auto x756 = x254*x706;
  auto x757 = x184*x334;
  auto x758 = x300*x334;
  auto x759 = x192*x719;
  auto x760 = x181*x265;
  auto x761 = x274*x738;
  auto x762 = x223*x747;
  auto x763 = x260*x743;
  auto x764 = -2*x121*x753 + x16*x400 + x163*x394 + x223*x761 + x223*x763 + x224*x671*x733 + x224*x761 + x224*x763 + x249*x695 + x252*x755 + x256*x701 + x256*x702 + x256*x757 + x257*x700 + x257*x758 + x258*x396 + x258*x397 + x261*x708 + x263*x706 + x268*x269*x738 + x268*x72*x741 + x268*x736 + x268*x739 + x274*x749 + x274*x762 + x275*x724 + x276*x724 + x277*x751 + x277*x752 + x381*x662 + x393*x754 + x394*x551 + x395*x704 + x395*x705 + x401*x759 + x407*x586 + x49*x756 + x57*x756 + x659*x718 + x665*x723 + x667*x731 + x668*x734 + x671*x744 + x673*x742 + x674*x743 + x678*x750 + x681*x742 + x721*x760;
  auto x765 = -x238*x689 - x241*x694 - 1.7320937034773347*x687 - x691*x692 + x764;
  auto x766 = x535*x712;
  auto x767 = -x238*x715 - x241*x729 - x692*x727 - 1.7320937034773347*x766;
  auto x768 = 10.385160084940418*x685 + x765 + x767;
  auto x769 = x145*x147;
  auto x770 = x147*x428;
  auto x771 = x158*x469;
  auto x772 = x2*x475;
  auto x773 = x168*x443;
  auto x774 = x174*x175;
  auto x775 = x434*x73;
  auto x776 = x174*x427;
  auto x777 = x434*x74;
  auto x778 = x197*x448;
  auto x779 = 0.69987154478819702*x201;
  auto x780 = x229*x438;
  auto x781 = x443*x96;
  auto x782 = x208*x781;
  auto x783 = x438*x74;
  auto x784 = x211*x783;
  auto x785 = x217*x781;
  auto x786 = 0.1166452574646995*x230;
  auto x787 = 0.23329051492939901*x233;
  auto x788 = x152*x770 + x162*x771 + x162*x772 - x166*x773 + x182*x775 + x183*x776 + x194*x479 + x194*x480 + x205*x780 + x205*x782 + x211*x785 + x222*x449 + x222*x453 + x225*x449 + x225*x453 + x445*x778 - x457*x779 + x457*x786 + x457*x787 + x463*x774 + x464*x769 + x583*x777 + x609*x784;
  auto x789 = -x121*x788;
  auto x790 = 3.461720028313473*x239;
  auto x791 = x201*x790;
  auto x792 = x463*x482;
  auto x793 = x174*x257;
  auto x794 = 0.57695333805224558*x230;
  auto x795 = 1.1539066761044912*x233;
  auto x796 = x147*x792 + x256*x770 + x258*x479 + x258*x480 + x261*x776 + x268*x780 + x268*x782 + x271*x785 + x275*x449 + x275*x453 + x276*x449 + x276*x453 + x395*x771 + x395*x772 + x483*x793 + x484*x778 + x485*x794 + x485*x795 + x666*x777 + x675*x783 + x760*x775;
  auto x797 = -x238*x773 - x456*x791 + x789 + x796;
  auto x798 = 1.0/x21;
  auto x799 = (2.0/3.0)*x798;
  auto x800 = (1.0/6.0)*x21 - x799;
  auto x801 = std::tan(f[1]);
  auto x802 = 1.0/x801;
  auto x803 = (2.0/3.0)*x802;
  auto x804 = -x478*x800 - x478*x803 + x788;
  auto x805 = x121*x804;
  auto x806 = -x488*x800 - x488*x803 - x805;
  auto x807 = x797 + x806;
  auto x808 = x168*x502;
  auto x809 = x158*x514;
  auto x810 = x2*x517;
  auto x811 = x147*x496;
  auto x812 = x197*x506;
  auto x813 = x174*x495;
  auto x814 = x24*x510;
  auto x815 = x498*x73;
  auto x816 = x498*x74;
  auto x817 = x229*x500;
  auto x818 = x502*x96;
  auto x819 = x208*x818;
  auto x820 = x500*x74;
  auto x821 = x217*x818;
  auto x822 = x147*x524 + x256*x811 + x258*x521 + x258*x522 + x261*x813 + x268*x817 + x268*x819 + x271*x821 + x275*x505 + x275*x507 + x276*x505 + x276*x507 + x395*x809 + x395*x810 + x484*x812 + x525*x794 + x525*x795 + x666*x816 + x675*x820 + x760*x815 + x793*x814;
  auto x823 = x211*x820;
  auto x824 = x152*x811 + x162*x809 + x162*x810 - x166*x808 + x182*x815 + x183*x813 + x194*x521 + x194*x522 + x205*x817 + x205*x819 + x211*x821 + x222*x505 + x222*x507 + x225*x505 + x225*x507 + x445*x812 - x509*x779 + x509*x786 + x509*x787 + x510*x774 + x511*x769 + x583*x816 + x609*x823;
  auto x825 = (1.0/3.0)*x478;
  auto x826 = x802*x825;
  auto x827 = x798*x825;
  auto x828 = -x8*x826 + x8*x827 + x824;
  auto x829 = x121*x828;
  auto x830 = (1.0/3.0)*x488;
  auto x831 = x802*x830;
  auto x832 = x798*x830;
  auto x833 = -x121*x824 - x8*x831 + x8*x832 - x829;
  auto x834 = -x238*x808 - x508*x791 + x822 + x833;
  auto x835 = x788 - x826 + x827;
  auto x836 = x121*x835;
  auto x837 = -x831 + x832 - x836;
  auto x838 = x797 + x837;
  auto x839 = -x121*x411 - x355*x531 - x369*x532 + x387*x530 - x413;
  auto x840 = x358*x71;
  auto x841 = x237*x89;
  auto x842 = 6.9283748139093388*x841;
  auto x843 = x364*x688;
  auto x844 = x378*x540;
  auto x845 = x101*x326;
  auto x846 = x326*x9;
  auto x847 = x332*x334;
  auto x848 = x334*x337;
  auto x849 = -x22 - 4*x51;
  auto x850 = x55*x849;
  auto x851 = x331*x396;
  auto x852 = x331*x397;
  auto x853 = x329 + 1;
  auto x854 = 0.22588256978742899*x853;
  auto x855 = 10*x101 - x102 - 10*x9;
  auto x856 = x107*x855;
  auto x857 = 1.0/x352;
  auto x858 = x353*x9;
  auto x859 = x857*x858;
  auto x860 = x560*x859;
  auto x861 = 1.9795356611169539*x165;
  auto x862 = x861*x89;
  auto x863 = x592*x853;
  auto x864 = -x361 - 1.0/4.0;
  auto x865 = x68*x864;
  auto x866 = std::pow(x33, 3.0/2.0);
  auto x867 = x362/x866;
  auto x868 = x865*x867;
  auto x869 = x356*x364;
  auto x870 = x719*x869;
  auto x871 = x601*x853;
  auto x872 = x373*x724;
  auto x873 = x588*x724;
  auto x874 = x364*x71;
  auto x875 = x243*x874;
  auto x876 = x859*x90;
  auto x877 = x165*x89;
  auto x878 = 0.082480652546539746*x877;
  auto x879 = x573*x859;
  auto x880 = 0.082480652546539746*x378;
  auto x881 = x378*x599;
  auto x882 = x223*x859;
  auto x883 = 0.10101775619540625*x588;
  auto x884 = x224*x859;
  auto x885 = x72*x868;
  auto x886 = x74*x885;
  auto x887 = x13*x869;
  auto x888 = x379*x735;
  auto x889 = 0.058322628732349752*x853;
  auto x890 = x74*x887;
  auto x891 = x378*x876;
  auto x892 = x887*x90;
  auto x893 = x101*x324 + x108*x871 + x116*x871 - x15*x318 - x162*x49 + x186*x870 + x188*x868 + x200*x856 + x200*x872 + x211*x888 + x219*x886 - x316*x8 - x324*x9 + x325*x845 - x325*x846 + x346*x850 + x382*x886 + x382*x890 + x386*x892 + x547*x847 + x549*x848 + x552*x851 + x552*x852 + x559*x853 + x565*x853 + x567*x853 + x570*x854 + x572*x854 + x607*x885 + x607*x887 + x612*x886 + x612*x890 - x614*x891 + x615*x889 + x617*x885 + x621 - x81*x860 - x840*x862 - x840*x881 - x843*x862 - x843*x881 + x863*x92 + x863*x99 + x873*x875 - x876*x878 - x879*x880 - x882*x883 - x883*x884;
  auto x894 = x639 + x893;
  auto x895 = x121*x894;
  auto x896 = x21*x395;
  auto x897 = x653*x853;
  auto x898 = x667*x853;
  auto x899 = x671*x853;
  auto x900 = x381*x724;
  auto x901 = x187*x68;
  auto x902 = x864*x867;
  auto x903 = x901*x902;
  auto x904 = x378*x859;
  auto x905 = x671*x74;
  auto x906 = x260*x886;
  auto x907 = x260*x890;
  auto x908 = x101*x391 + x108*x899 + x116*x899 + x223*x906 + x223*x907 + x224*x906 + x224*x907 + x260*x856 + x260*x872 + x266*x870 + x266*x903 + x271*x888 - x391*x9 + x393*x845 - x393*x846 + x49*x897 + x57*x897 + x648*x851 + x648*x852 + x650*x848 + x651*x847 + x656*x853 + x658*x853 + x660*x853 + x663*x900 - x668*x904 + x673*x885 + x673*x887 + x674*x886 - x678*x891 + x679*x853 + x681*x885 + x681*x887 - x82*x84*x859 - 0.28868228391288908*x841*x876 + x850*x896 - x882*x905 - x884*x905 + x898*x92 + x898*x99;
  auto x909 = -x121*x893 - x627 - x632 - x633 + x647 - x843*x844 - x895 + x908;
  auto x910 = x175*x331;
  auto x911 = x331*x427;
  auto x912 = x334*x428;
  auto x913 = x13*x468;
  auto x914 = x349*x5;
  auto x915 = x473*x914;
  auto x916 = x368*x438;
  auto x917 = x358*x777;
  auto x918 = x375*x777;
  auto x919 = x366*x452;
  auto x920 = 0.69987154478819702*x369;
  auto x921 = x783*x90;
  auto x922 = x378*x443;
  auto x923 = x375*x922;
  auto x924 = 0.1166452574646995*x385;
  auto x925 = 0.23329051492939901*x409;
  auto x926 = x152*x912 - x166*x916 + x183*x911 + x186*x918 + x211*x923 + x341*x479 + x341*x480 + x343*x913 + x346*x915 + x357*x917 + x372*x921 + x372*x922 + x376*x784 + x380*x449 + x380*x453 + x445*x919 + x449*x725 + x453*x725 - x457*x920 + x457*x924 + x457*x925 + x463*x910 + x464*x707;
  auto x927 = -x121*x926;
  auto x928 = x369*x790;
  auto x929 = x395*x43;
  auto x930 = x257*x331;
  auto x931 = x268*x371;
  auto x932 = x268*x373;
  auto x933 = x271*x376;
  auto x934 = x260*x381;
  auto x935 = 0.57695333805224558*x385;
  auto x936 = 1.1539066761044912*x409;
  auto x937 = x256*x912 + x261*x911 + x266*x918 + x271*x923 + x334*x792 + x394*x479 + x394*x480 + x401*x917 + x407*x449 + x407*x453 + x449*x934 + x453*x934 + x483*x930 + x484*x919 + x485*x935 + x485*x936 + x781*x932 + x783*x933 + x896*x915 + x913*x929 + x921*x931;
  auto x938 = -x238*x916 - x456*x928 + x927 + x937;
  auto x939 = std::tan(f[2]);
  auto x940 = 1.0/x939;
  auto x941 = x825*x940;
  auto x942 = 1.0/x22;
  auto x943 = (1.0/3.0)*x942;
  auto x944 = x520*x943;
  auto x945 = -x2*x941 + x2*x944 + x926;
  auto x946 = x121*x945;
  auto x947 = x830*x940;
  auto x948 = x527*x943;
  auto x949 = -x2*x947 + x2*x948 - x946;
  auto x950 = x938 + x949;
  auto x951 = x368*x500;
  auto x952 = x13*x513;
  auto x953 = x516*x914;
  auto x954 = x331*x495;
  auto x955 = x334*x496;
  auto x956 = x366*x504;
  auto x957 = x358*x816;
  auto x958 = x375*x816;
  auto x959 = x820*x90;
  auto x960 = x378*x502;
  auto x961 = x375*x960;
  auto x962 = x256*x955 + x261*x954 + x266*x958 + x271*x961 + x334*x524 + x394*x521 + x394*x522 + x401*x957 + x407*x505 + x407*x507 + x484*x956 + x505*x934 + x507*x934 + x525*x935 + x525*x936 + x814*x930 + x818*x932 + x820*x933 + x896*x953 + x929*x952 + x931*x959;
  auto x963 = x152*x955 - x166*x951 + x183*x954 + x186*x958 + x211*x961 + x341*x521 + x341*x522 + x343*x952 + x346*x953 + x357*x957 + x372*x959 + x372*x960 + x376*x823 + x380*x505 + x380*x507 + x445*x956 + x505*x725 + x507*x725 - x509*x920 + x509*x924 + x509*x925 + x510*x910 + x511*x707;
  auto x964 = (2.0/3.0)*x942;
  auto x965 = (1.0/6.0)*x22 - x964;
  auto x966 = (2.0/3.0)*x940;
  auto x967 = -x478*x965 - x520*x966 + x963;
  auto x968 = x121*x967;
  auto x969 = -x121*x963 - x488*x965 - x527*x966 - x968;
  auto x970 = -x238*x951 - x508*x928 + x962 + x969;
  auto x971 = x926 - x941 + x944;
  auto x972 = x121*x971;
  auto x973 = -x947 + x948 - x972;
  auto x974 = x938 + x973;
  auto x975 = -x121*x487 + x478*x530 - x489 + x493;
  auto x976 = (1.0/6.0)*x42;
  auto x977 = x110 + 2;
  auto x978 = (1.0/18.0)*f[0];
  auto x979 = x977*x978;
  auto x980 = 1.0*x29;
  auto x981 = 1.0*x27;
  auto x982 = x980 - x981;
  auto x983 = x175*x34;
  auto x984 = 2.25*x63;
  auto x985 = 2.25*x61;
  auto x986 = x984 - x985;
  auto x987 = x34*x986;
  auto x988 = 0.25*x79;
  auto x989 = 0.25*x77;
  auto x990 = x988 - x989;
  auto x991 = x97*x990;
  auto x992 = 1.0*x45;
  auto x993 = 1.0*x44;
  auto x994 = x992 - x993;
  auto x995 = x47*x994;
  auto x996 = 1.0*x54;
  auto x997 = 1.0*x53;
  auto x998 = x996 - x997;
  auto x999 = x474*x998;
  auto x1000 = 0.25*x87;
  auto x1001 = 0.25*x86;
  auto x1002 = x1000 - x1001;
  auto x1003 = x1002*x439;
  auto x1004 = 0.25*x94;
  auto x1005 = 0.25*x93;
  auto x1006 = x1004 - x1005;
  auto x1007 = x1006*x98;
  auto x1008 = 0.25*x113;
  auto x1009 = 0.25*x112;
  auto x1010 = x1008 - x1009;
  auto x1011 = x1010*x111;
  auto x1012 = 0.25*x105;
  auto x1013 = 0.25*x104;
  auto x1014 = x1012 - x1013;
  auto x1015 = x1014*x103;
  auto x1016 = 0.25*x118;
  auto x1017 = 0.25*x117;
  auto x1018 = x1016 - x1017;
  auto x1019 = x202*x458;
  auto x1020 = x1003*x211 + x1007*x211 + x1011*x445 + x1015*x445 + x1018*x1019 + x183*x987 + x186*x991 + x343*x995 + x346*x999 + x982*x983;
  auto x1021 = x1020 + x235*x976 + x637*x979;
  auto x1022 = x1021*x121;
  auto x1023 = x643*x979;
  auto x1024 = x24*x482;
  auto x1025 = x239*x486;
  auto x1026 = x1003*x271 + x1007*x271 + x1011*x484 + x1015*x484 + x1018*x1025 + x1024*x982 + x261*x987 + x266*x991 + x896*x999 + x929*x995;
  auto x1027 = -x1020*x121 + x1026;
  auto x1028 = -x1022 + x1023 + x1027 + x281*x976;
  auto x1029 = -x980 + x981;
  auto x1030 = -x984 + x985;
  auto x1031 = x1030*x34;
  auto x1032 = -x988 + x989;
  auto x1033 = x1032*x97;
  auto x1034 = -x992 + x993;
  auto x1035 = x1034*x47;
  auto x1036 = -x996 + x997;
  auto x1037 = x1036*x474;
  auto x1038 = -x1000 + x1001;
  auto x1039 = x1038*x439;
  auto x1040 = -x1004 + x1005;
  auto x1041 = x1040*x98;
  auto x1042 = -x1012 + x1013;
  auto x1043 = x103*x1042;
  auto x1044 = -x1008 + x1009;
  auto x1045 = x1044*x111;
  auto x1046 = -x1016 + x1017;
  auto x1047 = x1019*x1046 + x1029*x983 + x1031*x183 + x1033*x186 + x1035*x343 + x1037*x346 + x1039*x211 + x1041*x211 + x1043*x445 + x1045*x445;
  auto x1048 = -x1047*x121;
  auto x1049 = (1.0/3.0)*x22;
  auto x1050 = x1049*x387;
  auto x1051 = (1.0/3.0)*x21;
  auto x1052 = x1051*x235;
  auto x1053 = (1.0/9.0)*f[0];
  auto x1054 = x1053*x637;
  auto x1055 = x1054*x2;
  auto x1056 = x1047 - x1050*x2 - x1052*x8 + x1055*x8;
  auto x1057 = x1056*x121;
  auto x1058 = x1049*x412;
  auto x1059 = x1051*x281;
  auto x1060 = x1053*x643;
  auto x1061 = x1060*x2;
  auto x1062 = x1024*x1029 + x1025*x1046 + x1031*x261 + x1033*x266 + x1035*x929 + x1037*x896 + x1039*x271 + x1041*x271 + x1043*x484 + x1045*x484;
  auto x1063 = x1061*x8 + x1062;
  auto x1064 = x1048 - x1057 - x1058*x2 - x1059*x8 + x1063;
  auto x1065 = x1020 - x1052 + x1055;
  auto x1066 = x1065*x121;
  auto x1067 = x1027 - x1059 + x1061 - x1066;
  auto x1068 = -x121*x526 + x520*x530 - x528 + x529;
  auto x1069 = (1.0/6.0)*x51;
  auto x1070 = x101 + 2;
  auto x1071 = x1070*x978;
  auto x1072 = x1020 + x1069*x387 + x1071*x637;
  auto x1073 = x1072*x121;
  auto x1074 = x1071*x643;
  auto x1075 = x1027 + x1069*x412 - x1073 + x1074;
  auto x1076 = x1047 - x1050 + x1054*x8;
  auto x1077 = x1076*x121;
  auto x1078 = x1060*x8 + x1062;
  auto x1079 = x1048 - x1058 - x1077 + x1078;
  auto x1080 = x1020 + x1054;
  auto x1081 = x1080*x121;
  auto x1082 = x1027 + x1060 - x1081;
  auto x1083 = 6.928374813909338*x535;
  auto x1084 = 6.9234400566269469*x240;
  auto x1085 = -x1083*x534 - x1083*x538 - x1084*x539;
  auto x1086 = 10.38516008494042*x685;
  auto x1087 = 3.464187406954669*x237;
  auto x1088 = 1.7308600141567367*x240;
  auto x1089 = 3.4617200283134735*x240;
  auto x1090 = x1086 - x1087*x715 - x1088*x727 - x1089*x729 - 1.7320937034773345*x766;
  auto x1091 = 3.4617200283134735*x201;
  auto x1092 = -x1087*x773 - x1091*x485 + x796;
  auto x1093 = x1092 + x789;
  auto x1094 = -x1087*x808 - x1091*x525 + x822;
  auto x1095 = 0.134690341593875*x242;
  auto x1096 = x315*x8;
  auto x1097 = x122*x2;
  auto x1098 = 0.67345170796937492*x242;
  auto x1099 = x15*x250;
  auto x1100 = x20*x250;
  auto x1101 = 0.37647094964571493*x652;
  auto x1102 = 1.882354748228575*x652;
  auto x1103 = f[0]*x4;
  auto x1104 = 0.72075922005612658*x246;
  auto x1105 = x17*x248;
  auto x1106 = 0.38490971188385215*x237;
  auto x1107 = 0.7760975132133765*x259;
  auto x1108 = -x108*x1107 + x11*x1104 + 0.11111111111111099*x1103 + x1104*x7 - 1.9196908540199586*x1105 - x1106*x92 - x1106*x99 - x1107*x116 + 0.11398633159763*x624 + x644 - 0.23329051492939901*x645 - 0.38463555870149707*x646 - 0.8748394309852463*x655 - 0.87483943098524752*x657;
  auto x1109 = -x1087*x169 - x1091*x240 + x279;
  auto x1110 = x1109 + x236;
  auto x1111 = (1.0/3.0)*f[0];
  auto x1112 = 6.2598413986567918*x136;
  auto x1113 = x165*x169;
  auto x1114 = x33*x556;
  auto x1115 = x262*x558;
  auto x1116 = x1114*x1115;
  auto x1117 = x295*x566;
  auto x1118 = x1114*x1117;
  auto x1119 = 0.22588256978742899*x1114;
  auto x1120 = x1119*x569;
  auto x1121 = x1119*x571;
  auto x1122 = x201*x203;
  auto x1123 = x48*x545;
  auto x1124 = 1.355295418724574*x174;
  auto x1125 = x264*x560;
  auto x1126 = std::pow(x2, 3)/std::pow(x21, 5.0/2.0);
  auto x1127 = x154*x160;
  auto x1128 = 0.67764770936228702*x556;
  auto x1129 = x154*x255;
  auto x1130 = x180*x3;
  auto x1131 = x67*x69;
  auto x1132 = x556*x560;
  auto x1133 = x1132*x187;
  auto x1134 = 0.60610653717243745*x199;
  auto x1135 = 0.52490365859114851*x144;
  auto x1136 = x174*x262;
  auto x1137 = 1.049807317182297*x144;
  auto x1138 = 0.34993577239409851*x151;
  auto x1139 = 0.69987154478819702*x58;
  auto x1140 = 0.69987154478819702*x151;
  auto x1141 = x207*x592;
  auto x1142 = x597*x618;
  auto x1143 = 1.355295418724574*x154*x577;
  auto x1144 = 0.15152663429310936*x223;
  auto x1145 = x1126*x199*x73;
  auto x1146 = 0.15152663429310936*x224;
  auto x1147 = x114*x554;
  auto x1148 = 1.2122130743448749*x199;
  auto x1149 = 0.12372097881980962*x1126;
  auto x1150 = x147*x30;
  auto x1151 = x307*x601;
  auto x1152 = x208*x745;
  auto x1153 = x878*x90;
  auto x1154 = x556*x72;
  auto x1155 = x1131*x1154;
  auto x1156 = x573*x880;
  auto x1157 = x1114*x173;
  auto x1158 = 0.67764770936228702*x1157;
  auto x1159 = x171*x575;
  auto x1160 = x537*x582;
  auto x1161 = x189*x580*x69;
  auto x1162 = 0.52490365859114774*x203;
  auto x1163 = 0.30305326858621873*x199;
  auto x1164 = x1134*x217;
  auto x1165 = x556*x58;
  auto x1166 = 0.34993577239409851*x58;
  auto x1167 = x216*x561;
  auto x1168 = x1155*x883;
  auto x1169 = x172*x33*x575*(-3.0/2.0*x170 - 3.0/2.0)/std::pow(x31, 4);
  auto x1170 = x230*x614;
  auto x1171 = x183*x64;
  auto x1172 = x581*(-5.0/4.0*x170 - 5.0/4.0)/std::pow(x31, 7.0/2.0);
  auto x1173 = x1169*x155;
  auto x1174 = 0.087483943098524627*x203;
  auto x1175 = x1154*x1161;
  auto x1176 = x614*x680;
  auto x1177 = 4.1992292687291819*x203;
  auto x1178 = x1177*x2;
  auto x1179 = x199*x556;
  auto x1180 = x1163*x1167;
  auto x1181 = 0.12372097881980962*x556;
  auto x1182 = x165*x207;
  auto x1183 = x596*x618;
  auto x1184 = 0.24744195763961924*x165;
  auto x1185 = x1184*x593;
  auto x1186 = x573*x595;
  auto x1187 = x1175*x883;
  auto x1188 = x1172*x72;
  auto x1189 = x1188*x74;
  auto x1190 = x1134*x223;
  auto x1191 = x181*x606;
  auto x1192 = x1134*x224;
  auto x1193 = x206*x90;
  auto x1194 = x606*x618;
  auto x1195 = 0.17496788619704925*x203;
  auto x1196 = x1195*x595;
  auto x1197 = x135 - x149 - x153;
  auto x1198 = x1109*x638 + x1197;
  auto x1199 = 3.4617200283134735*x369;
  auto x1200 = -x1087*x355 - x1199*x240 + x410;
  auto x1201 = x322 - x336 - x340;
  auto x1202 = x1200*x638 + x1201;
  auto x1203 = x162*x545;
  auto x1204 = x165*x355;
  auto x1205 = x547*x706;
  auto x1206 = x16*x547;
  auto x1207 = x147*x549;
  auto x1208 = x16*x549;
  auto x1209 = x1177*x167;
  auto x1210 = x552*x706;
  auto x1211 = x203*x369;
  auto x1212 = x556*x558;
  auto x1213 = x334*x35;
  auto x1214 = x356*x359;
  auto x1215 = x1214*x560;
  auto x1216 = x556*x566;
  auto x1217 = x331*x568;
  auto x1218 = x34*x396;
  auto x1219 = x154*x568;
  auto x1220 = x34*x397;
  auto x1221 = x330*x576;
  auto x1222 = x1221*x155;
  auto x1223 = x1153*x371;
  auto x1224 = x384*x732;
  auto x1225 = 0.058322628732349752*x58;
  auto x1226 = x371*x561;
  auto x1227 = x560*x561;
  auto x1228 = x166*x95;
  auto x1229 = x716*x74;
  auto x1230 = x192*x714;
  auto x1231 = x167*x735;
  auto x1232 = x1231*x192;
  auto x1233 = x537*x723;
  auto x1234 = x601*x602;
  auto x1235 = x420*x601;
  auto x1236 = x373*x745;
  auto x1237 = x16*x717;
  auto x1238 = x579*x580*x722;
  auto x1239 = x186*x72;
  auto x1240 = x1239*x723;
  auto x1241 = x221*x733;
  auto x1242 = x200*x724;
  auto x1243 = x200*x586;
  auto x1244 = x243*x71*x723;
  auto x1245 = 0.041240326273269873*x371;
  auto x1246 = x877*x90;
  auto x1247 = x1246*x371;
  auto x1248 = 0.041240326273269873*x556;
  auto x1249 = x384*x573;
  auto x1250 = x385*x614;
  auto x1251 = x204*x90;
  auto x1252 = x204*x618;
  auto x1253 = x376*x74;
  auto x1254 = 0.050508878097703123*x199;
  auto x1255 = x1226*x1254;
  auto x1256 = 0.050508878097703123*x1179;
  auto x1257 = x375*x561;
  auto x1258 = x582*x735;
  auto x1259 = x1258*x74;
  auto x1260 = x192*x356;
  auto x1261 = x16*x730;
  auto x1262 = x1261*x205*x90;
  auto x1263 = x16*x180;
  auto x1264 = x1260*x1263;
  auto x1265 = x1238*x72;
  auto x1266 = x1265*x74;
  auto x1267 = x212*x723;
  auto x1268 = x181*x742;
  auto x1269 = x221*x224;
  auto x1270 = 0.029161314366174876*x203;
  auto x1271 = x203*x385;
  auto x1272 = x1263*x96;
  auto x1273 = x228*x90;
  auto x1274 = x1272*x1273;
  auto x1275 = -x110*x699 + x1123*x341 - x1132*x402 + x1147*x380 + x1147*x725 + x1150*x1205 - x1156*x375*x556 - 0.058322628732349752*x1165*x1214 + x1171*x1221 + x1178*x368 + x1185*x354 - x1186*x1245 + x1193*x1267 + x1196*x354 + x1203*x344 + 0.24744195763961924*x1204*x556 + 0.49488391527923847*x1204 + x1206*x300 + x1207*x708 + x1208*x184 + x1209*x693 + x1210*x163 + x1210*x551 + 0.17496788619704925*x1211*x556 + 0.34993577239409851*x1211 - x1212*x1213 - x1212*x399 - x1215 - x1216*x339 - x1216*x398 - x1217*x569 - x1217*x571 - x1218*x1219 - x1219*x1220 + x1221*x176 + x1222*x49 + x1222*x57 - x1223 - x1224 - x1225*x1226*x80 - x1227*x721 - x1228*x1229 - x1228*x1232 - x1229*x1251 - x1230*x206 - x1232*x1251 - x1233*x574 - x1233*x599*x90 - x1234*x724 - x1235 - x1236 + x1237*x720 + x1238*x188 + x1240*x584 + x1241*x586 + x1242*x608 + x1242*x611 + x1243*x738 + x1243*x747 + x1244*x589 - x1245*x165*x593*x90 - x1247*x1248 - x1248*x1249 - x1250 - x1252*x693 - x1253*x2*x599 - x1255*x223 - x1255*x224 - x1256*x224*x373 - x1256*x420 - x1257*x223*x601 - x1257*x745 + x1258*x205*x370 + x1258*x228*x680 + x1259*x1269 + x1259*x210 + x1259*x748 + x1260*x1262 + x1260*x1274 + x1264*x1269 + x1264*x210 + x1264*x748 + x1265*x607 + x1265*x617 + x1266*x219 + x1266*x382 + x1266*x612 + x1268*x219 + x1268*x382 + x1268*x612 - x1270*x371*x595*x90 - 0.029161314366174876*x1271*x556 - x179*x402 - x2*x381*x574 - x204*x378*x582*x688 - x205*x405 - x205*x406 - x221*x367 - x221*x383 - x221*x408 - x228*x409 + x231*x618*x742*x90 + x3*x699 + x30*x577*x707 + x328 - x333 - x338 - x342 - x351 + x357*x582*x719 - x367*x605 - x375*x595*x732 - x376*x556*x878 - x376*x592*x593 - x376*x595*x614 - x383*x605 + x396*x585 + x397*x585 - x408*x605 - x409*x556*x614 + x553*x704 + x553*x705 + x578*x709 - x582*x714*x89;
  auto x1276 = x800*x835;
  auto x1277 = (4.0/3.0)*x802;
  auto x1278 = x221*x97;
  auto x1279 = -x1278*x449 - x1278*x453 - x179*x435 - x205*x440 - x205*x444 - x429 - x457*x620 - x465;
  auto x1280 = x1279 + x487*x638;
  auto x1281 = x463*x547;
  auto x1282 = x147*x174;
  auto x1283 = x445*x554;
  auto x1284 = x1212*x24;
  auto x1285 = x1219*x34;
  auto x1286 = x443*x861;
  auto x1287 = x175*x577;
  auto x1288 = x183*x577;
  auto x1289 = x1239*x582;
  auto x1290 = x181*x583;
  auto x1291 = x588*x591;
  auto x1292 = x438*x592;
  auto x1293 = x561*x91;
  auto x1294 = x592*x669;
  auto x1295 = 1.399743089576394*x457;
  auto x1296 = x605*x97;
  auto x1297 = x212*x582;
  auto x1298 = x181*x610;
  auto x1299 = x211*x611;
  auto x1300 = x200*x608;
  auto x1301 = x200*x611;
  auto x1302 = 0.058322628732349752*x613;
  auto x1303 = x120*x616;
  auto x1304 = 0.23329051492939901*x457;
  auto x1305 = x606*x680;
  auto x1306 = -x1132*x435 + x1203*x469 + x1207*x776 - x1216*x481 - x1227*x775 - x1234*x449 - x1234*x453 + x1281*x1282 + x1283*x448 - x1284*x464 - x1285*x479 - x1285*x480 - x1286*x534 - x1286*x538 + x1287*x463 + x1288*x427 + x1289*x777 + x1290*x434 + x1291*x778 - x1292*x1293 - x1294*x443 - x1295*x539 - x1295*x600 - x1296*x449 - x1296*x453 + x1297*x921 + x1297*x922 + x1298*x438 + x1299*x781 + x1300*x449 + x1300*x453 + x1301*x449 + x1301*x453 - x1302*x457 - x1303*x457 + x1304*x1305 + x1304*x619 - x440*x598 - x444*x598 - x476 + x479*x585 + x480*x585 + x553*x771 + x553*x772 + x587*x778;
  auto x1307 = x1280 + x1306;
  auto x1308 = x803*x835;
  auto x1309 = x799*x804;
  auto x1310 = -x1278*x505 - x1278*x507 - x179*x499 - x205*x501 - x205*x503 - x497 - x509*x620 - x512;
  auto x1311 = x1310 + x526*x638;
  auto x1312 = x510*x547;
  auto x1313 = x502*x861;
  auto x1314 = x500*x592;
  auto x1315 = 1.399743089576394*x509;
  auto x1316 = 0.23329051492939901*x509;
  auto x1317 = -x1132*x499 + x1203*x514 + x1207*x813 - x1216*x523 - x1227*x815 - x1234*x505 - x1234*x507 + x1282*x1312 + x1283*x506 - x1284*x511 - x1285*x521 - x1285*x522 + x1287*x510 + x1288*x495 + x1289*x816 + x1290*x498 + x1291*x812 - x1293*x1314 - x1294*x502 - x1296*x505 - x1296*x507 + x1297*x959 + x1297*x960 + x1298*x500 + x1299*x818 + x1300*x505 + x1300*x507 + x1301*x505 + x1301*x507 - x1302*x509 - x1303*x509 + x1305*x1316 - x1313*x534 - x1313*x538 - x1315*x539 - x1315*x600 + x1316*x619 - x501*x598 - x503*x598 - x518 + x521*x585 + x522*x585 + x553*x809 + x553*x810 + x587*x812;
  auto x1318 = -x1087*x689 - x1088*x691 - x1089*x694 - 1.7320937034773345*x687 + x764;
  auto x1319 = x1090 + x1318;
  auto x1320 = x1200 + x388;
  auto x1321 = x1201 + x1275 + x1320*x638;
  auto x1322 = x162*x2;
  auto x1323 = x331*x552;
  auto x1324 = x558*x853;
  auto x1325 = x566*x853;
  auto x1326 = x334*x549;
  auto x1327 = x21*x850;
  auto x1328 = x1115*x853;
  auto x1329 = x264*x58;
  auto x1330 = x560*x853;
  auto x1331 = x1117*x853;
  auto x1332 = x174*x854;
  auto x1333 = x106*x855;
  auto x1334 = x5*x857*x9;
  auto x1335 = x1334*x81;
  auto x1336 = x70*x71;
  auto x1337 = x537*x868;
  auto x1338 = x688*x723;
  auto x1339 = x217*x221;
  auto x1340 = 0.041240326273269873*x853;
  auto x1341 = x189*x190*x902;
  auto x1342 = x356*x719;
  auto x1343 = x1254*x853;
  auto x1344 = x217*x871;
  auto x1345 = x1334*x192;
  auto x1346 = x1254*x181;
  auto x1347 = x1345*x883;
  auto x1348 = x181*x885;
  auto x1349 = x205*x96;
  auto x1350 = x1341*x72;
  auto x1351 = x219*x74;
  auto x1352 = x723*x735;
  auto x1353 = x203*x230;
  auto x1354 = x203*x889;
  auto x1355 = x1263*x869;
  auto x1356 = x1350*x74;
  auto x1357 = x1352*x74;
  auto x1358 = -x101*x697 - 0.24744195763961924*x1113*x853 + 0.49488391527923847*x1113 - 0.17496788619704925*x1122*x853 + 0.34993577239409851*x1122 - x1125 + x1127*x854 + x1129*x854 - x1141 - x1142 - x1151 - x1152 - x1156*x1345 - x1170 - x1176*x1345 + x1182*x1340 + x1183*x1340 + x1193*x740*x868 + x1195*x167*x876 + x1205*x30*x334 + x1206*x332 + x1208*x337 + x1209*x840 + x1209*x843 + x1210*x396 + x1210*x397 - x1225*x584*x859 - x1228*x1337 - x1231*x166*x377 + x1237*x357*x364 + x1241*x724 + x1242*x738 + x1242*x747 + x1243*x886 + x1243*x890 + x1244*x873 - x1251*x1337 - x1251*x167*x887 - x1252*x843 - x1261*x1336*x166 + x1262*x869 + x1263*x1349*x356*x377 + x1269*x1348 + x1269*x1355 - x1270*x618*x876 - x1272*x1336*x204 + x1273*x618*x885 + x1274*x869 - x13*x213*x71*x861 - x13*x590*x881 + x1322*x850 + x1323*x704 + x1323*x705 + x1324*x148 + x1325*x150 + x1326*x708 + x1327*x194 + x1328*x174 + x1329*x889 + x1330*x267 + x1331*x174 + x1332*x569 + x1332*x571 + x1333*x222 + x1333*x225 - x1334*x214*x592 - x1335*x192*x560 - x1338*x862 - x1338*x881 - x1339*x223 - x1339*x224 + x1341*x188 + x1342*x186*x723 + x1343*x208*x224 + x1343*x307 + x1344*x223 + x1344*x224 - x1346*x882 - x1346*x884 - x1347*x223 - x1347*x224 + x1348*x210 + x1348*x748 + x1350*x1351 + x1350*x607 + x1350*x617 + x1351*x1352 + x1352*x607 + x1352*x617 + 0.029161314366174876*x1353*x853 + x1354*x233 + x1355*x748 + x1356*x382 + x1356*x612 + x1357*x382 + x1357*x612 + x143 - x161 - 0.041240326273269873*x165*x206*x876 + 0.24744195763961924*x167*x879 - x177 - x179*x267 + x182*x80*x885 - x185 - x195 - x198*x221 + x198*x871 - x205*x272 - x206*x364*x714 - x210*x217 - x228*x233 + x232*x732*x853 + x272*x863 - 0.10101775619540625*x589*x859 - 0.041240326273269873*x618*x879 + x697*x9;
  auto x1359 = x1110*x638 + x1197 + x1358;
  auto x1360 = (1.0/3.0)*x940;
  auto x1361 = x1360*x835;
  auto x1362 = x828*x943;
  auto x1363 = x331*x769;
  auto x1364 = x145*x16;
  auto x1365 = x147*x152;
  auto x1366 = x152*x16;
  auto x1367 = 2.0996146343645909*x684;
  auto x1368 = x158*x162;
  auto x1369 = x175*x706;
  auto x1370 = x174*x707;
  auto x1371 = x183*x706;
  auto x1372 = x194*x43;
  auto x1373 = x194*x21;
  auto x1374 = x690*x711;
  auto x1375 = x711*x712;
  auto x1376 = x166*x728;
  auto x1377 = x560*x733;
  auto x1378 = x13*x720;
  auto x1379 = x182*x375;
  auto x1380 = 0.34993577239409851*x457;
  auto x1381 = 0.69987154478819702*x457;
  auto x1382 = x733*x90;
  auto x1383 = x592*x733;
  auto x1384 = x205*x737;
  auto x1385 = x181*x205*x376;
  auto x1386 = x205*x747;
  auto x1387 = x601*x733;
  auto x1388 = x221*x747;
  auto x1389 = x200*x743;
  auto x1390 = 0.058322628732349752*x750;
  auto x1391 = 0.1166452574646995*x457;
  auto x1392 = x680*x742;
  auto x1393 = -x1230*x783 + x1240*x777 + x1267*x921 + x1267*x922 + x1292*x1382 + x1304*x1392 + x1322*x915 + x1363*x463 + x1364*x464 + x1365*x911 + x1366*x428 + x1367*x457 + x1368*x913 + x1369*x463 + x1370*x463 + x1371*x427 + x1372*x913 + x1373*x915 - x1374*x438 - x1375*x443 - x1376*x443 + x1377*x434 + x1378*x777 + x1379*x434 - x1380*x691 - x1380*x727 - x1381*x694 - x1381*x729 + x1383*x781 + x1384*x921 + x1384*x922 + x1385*x438 + x1386*x781 + x1387*x449 + x1387*x453 + x1388*x449 + x1388*x453 + x1389*x449 + x1389*x453 + x1390*x457 + x1391*x751 + x1391*x752 + x222*x919 + x225*x919 + x341*x771 + x341*x772 + x380*x778 + x449*x746 + x453*x746 + x479*x710 + x480*x710 + x709*x776 + x725*x778;
  auto x1394 = -x1361*x2 + x1362*x2 + x1393;
  auto x1395 = -x800*x971 - x803*x945;
  auto x1396 = x1394 + x1395;
  auto x1397 = (1.0/3.0)*x802;
  auto x1398 = x1397*x971;
  auto x1399 = (1.0/3.0)*x798;
  auto x1400 = x1399*x945;
  auto x1401 = x1363*x510;
  auto x1402 = x1364*x511;
  auto x1403 = x1365*x954;
  auto x1404 = x1366*x496;
  auto x1405 = x1367*x509;
  auto x1406 = x341*x809;
  auto x1407 = x341*x810;
  auto x1408 = x1368*x952;
  auto x1409 = x1322*x953;
  auto x1410 = x1369*x510;
  auto x1411 = x1370*x510;
  auto x1412 = x1371*x495;
  auto x1413 = x709*x813;
  auto x1414 = x521*x710;
  auto x1415 = x522*x710;
  auto x1416 = x1372*x952;
  auto x1417 = x1373*x953;
  auto x1418 = -x1374*x500;
  auto x1419 = -x1375*x502;
  auto x1420 = -x1230*x820;
  auto x1421 = -x1376*x502;
  auto x1422 = x1377*x498;
  auto x1423 = x1378*x816;
  auto x1424 = x1379*x498;
  auto x1425 = x1240*x816;
  auto x1426 = x380*x812;
  auto x1427 = x222*x956;
  auto x1428 = x225*x956;
  auto x1429 = x725*x812;
  auto x1430 = 0.34993577239409851*x509;
  auto x1431 = -x1430*x691;
  auto x1432 = -x1430*x727;
  auto x1433 = 0.69987154478819702*x509;
  auto x1434 = -x1433*x694;
  auto x1435 = -x1433*x729;
  auto x1436 = x1314*x1382;
  auto x1437 = x1383*x818;
  auto x1438 = x1384*x959;
  auto x1439 = x1384*x960;
  auto x1440 = x1385*x500;
  auto x1441 = x1386*x818;
  auto x1442 = x1267*x959;
  auto x1443 = x1267*x960;
  auto x1444 = x1387*x505;
  auto x1445 = x1387*x507;
  auto x1446 = x505*x746;
  auto x1447 = x507*x746;
  auto x1448 = x1388*x505;
  auto x1449 = x1388*x507;
  auto x1450 = x1389*x505;
  auto x1451 = x1389*x507;
  auto x1452 = x1390*x509;
  auto x1453 = 0.1166452574646995*x509;
  auto x1454 = x1453*x751;
  auto x1455 = x1453*x752;
  auto x1456 = x1316*x1392;
  auto x1457 = -x1398*x8 + x1400*x8 + x1401 + x1402 + x1403 + x1404 + x1405 + x1406 + x1407 + x1408 + x1409 + x1410 + x1411 + x1412 + x1413 + x1414 + x1415 + x1416 + x1417 + x1418 + x1419 + x1420 + x1421 + x1422 + x1423 + x1424 + x1425 + x1426 + x1427 + x1428 + x1429 + x1431 + x1432 + x1434 + x1435 + x1436 + x1437 + x1438 + x1439 + x1440 + x1441 + x1442 + x1443 + x1444 + x1445 + x1446 + x1447 + x1448 + x1449 + x1450 + x1451 + x1452 + x1454 + x1455 + x1456;
  auto x1458 = -x828*x966 - x835*x965;
  auto x1459 = x1457 + x1458;
  auto x1460 = -x1361 + x1362 + x1393;
  auto x1461 = -x1398 + x1400;
  auto x1462 = x1460 + x1461;
  auto x1463 = x1092 - x487*x800 - x487*x803 - 2*x805;
  auto x1464 = x2/x110;
  auto x1465 = std::pow(x801, 2);
  auto x1466 = (-x1465 - 1)/x1465;
  auto x1467 = x1279 + x488*x638;
  auto x1468 = x1306 + x1467;
  auto x1469 = -x1276 - 2.0/3.0*x1466*x478 + x1468 - x478*((2.0/3.0)*x1464 + (1.0/6.0)*x2) - x788*x800 - x788*x803 - x803*x804;
  auto x1470 = x1394 - x800*x926 - x803*x926;
  auto x1471 = x1065*x800;
  auto x1472 = x34*x769;
  auto x1473 = x166*x168;
  auto x1474 = x182*x73;
  auto x1475 = x174*x183;
  auto x1476 = x583*x74;
  auto x1477 = x197*x445;
  auto x1478 = x202*x779;
  auto x1479 = x205*x229;
  auto x1480 = x1349*x208;
  auto x1481 = x1002*x90;
  auto x1482 = x211*x217;
  auto x1483 = x211*x232;
  auto x1484 = x1018*x202;
  auto x1485 = x1002*x1479 - x1006*x1473 + x1006*x1480 + x1006*x1483 + x1010*x1477 + x1011*x222 + x1011*x225 + x1015*x222 + x1015*x225 - x1018*x1478 + x1322*x999 + x1365*x987 + x1368*x995 + x1372*x995 + x1373*x999 + x1472*x982 + x1474*x990 + x1475*x986 + x1476*x990 + x1481*x1482 + x1484*x786 + x1484*x787 + x774*x982;
  auto x1486 = -x1020*x800 - x1020*x803 + x1485;
  auto x1487 = -x1021*x803 + x1110*x979 - x1471 + x1486 + x640*x976;
  auto x1488 = x1065*x802;
  auto x1489 = (1.0/3.0)*x1488;
  auto x1490 = x1021*x1399;
  auto x1491 = x1038*x90;
  auto x1492 = x1046*x202;
  auto x1493 = x1029*x1472 + x1029*x774 + x1030*x1475 + x1031*x1365 + x1032*x1474 + x1032*x1476 + x1035*x1368 + x1035*x1372 + x1037*x1322 + x1037*x1373 + x1038*x1479 - x1040*x1473 + x1040*x1480 + x1040*x1483 + x1043*x222 + x1043*x225 + x1044*x1477 + x1045*x222 + x1045*x225 - x1046*x1478 + x1482*x1491 + x1492*x786 + x1492*x787;
  auto x1494 = -x1489*x8 + x1490*x8 + x1493;
  auto x1495 = x1049*x753;
  auto x1496 = x1051*x640;
  auto x1497 = x1053*x1110;
  auto x1498 = x1497*x2;
  auto x1499 = -x1495*x2 - x1496*x8 + x1498*x8;
  auto x1500 = -x1047*x800 - x1047*x803 + x1494 + x1499;
  auto x1501 = -x1496 + x1498;
  auto x1502 = x1486 - x1489 + x1490 + x1501;
  auto x1503 = x1397*x487;
  auto x1504 = x1399*x487;
  auto x1505 = x1094 - x1503*x8 + x1504*x8 - 2*x829;
  auto x1506 = x1397*x788;
  auto x1507 = x1397*x835;
  auto x1508 = x1399*x788;
  auto x1509 = x1399*x804;
  auto x1510 = x1464*x825;
  auto x1511 = x1466*x825;
  auto x1512 = x1310 + x527*x638;
  auto x1513 = x1317 - x1506*x8 - x1507*x8 + x1508*x8 + x1509*x8 - x1510*x8 - x1511*x8 + x1512;
  auto x1514 = x1397*x926;
  auto x1515 = x1399*x926;
  auto x1516 = x1401 + x1402 + x1403 + x1404 + x1405 + x1406 + x1407 + x1408 + x1409 + x1410 + x1411 + x1412 + x1413 + x1414 + x1415 + x1416 + x1417 + x1418 + x1419 + x1420 + x1421 + x1422 + x1423 + x1424 + x1425 + x1426 + x1427 + x1428 + x1429 + x1431 + x1432 + x1434 + x1435 + x1436 + x1437 + x1438 + x1439 + x1440 + x1441 + x1442 + x1443 + x1444 + x1445 + x1446 + x1447 + x1448 + x1449 + x1450 + x1451 + x1452 + x1454 + x1455 + x1456 + x1458 - x1514*x8 + x1515*x8 + x22*x826 - x22*x827;
  auto x1517 = x1020*x1397;
  auto x1518 = x1020*x1399;
  auto x1519 = -x1517*x8 + x1518*x8;
  auto x1520 = -x1056*x803 - x1076*x800;
  auto x1521 = x1493 + x1499 + x1519 + x1520;
  auto x1522 = x1047*x1397;
  auto x1523 = x1076*x1397;
  auto x1524 = x1047*x1399;
  auto x1525 = x1056*x1399;
  auto x1526 = x1069*x753 + x1485;
  auto x1527 = x1071*x1110 - x1522*x8 - x1523*x8 + x1524*x8 + x1525*x8 + x1526;
  auto x1528 = -x1495;
  auto x1529 = x1493 + x1497*x8 + x1528;
  auto x1530 = -x1523 + x1525;
  auto x1531 = x1519 + x1529 + x1530;
  auto x1532 = x1092 - x1503 + x1504 - 2*x836;
  auto x1533 = x1468 - x1506 - x1507 + x1508 + x1509 - x1510 - x1511;
  auto x1534 = x1460 - x1514 + x1515;
  auto x1535 = -x1080*x800;
  auto x1536 = x1485 - x1517 + x1518;
  auto x1537 = -2.0/3.0*x1488 + x1501 + x1535 + x1536;
  auto x1538 = x1080*x1397;
  auto x1539 = x1065*x1399;
  auto x1540 = -x1538*x8 + x1539*x8;
  auto x1541 = -x1522 + x1524 + x1529 + x1540;
  auto x1542 = x1497 + x1536 - x1538 + x1539;
  auto x1543 = 6.928374813909338*x841;
  auto x1544 = x1084*x378;
  auto x1545 = -x1543*x840 - x1543*x843 - x1544*x840;
  auto x1546 = -x1087*x916 - x1199*x485 + x937;
  auto x1547 = x1546 + x927;
  auto x1548 = -x1087*x951 - x1199*x525 + x962;
  auto x1549 = 1.355295418724574*x331;
  auto x1550 = x1327*x154;
  auto x1551 = x144*x399;
  auto x1552 = x151*x398;
  auto x1553 = 0.67764770936228702*x154*x853;
  auto x1554 = x32*x329;
  auto x1555 = x1330*x901;
  auto x1556 = x1554*x854;
  auto x1557 = x331*x854;
  auto x1558 = x353*std::pow(x8, 3)/std::pow(x22, 5.0/2.0);
  auto x1559 = x356*x858;
  auto x1560 = x72*x853;
  auto x1561 = x1153*x1560;
  auto x1562 = x68*x866;
  auto x1563 = x1560*x1562;
  auto x1564 = x1558*x588;
  auto x1565 = 0.12372097881980962*x1558;
  auto x1566 = x378*x573;
  auto x1567 = x1563*x883;
  auto x1568 = x688*x868;
  auto x1569 = x375*x378;
  auto x1570 = x13*x874;
  auto x1571 = x360*x865;
  auto x1572 = x1560*x1571;
  auto x1573 = x1354*x680*x72;
  auto x1574 = x362*x865*(-3.0/4.0*x329 - 3.0/4.0)/std::pow(x33, 5.0/2.0);
  auto x1575 = x199*x853;
  auto x1576 = 1.2122130743448749*x873;
  auto x1577 = 0.12372097881980962*x853;
  auto x1578 = x1572*x883;
  auto x1579 = x1334*x364;
  auto x1580 = 0.30305326858621873*x1579*x588;
  auto x1581 = x1574*x72;
  auto x1582 = x735*x868;
  auto x1583 = x1582*x74;
  auto x1584 = x1581*x74;
  auto x1585 = x966*x971;
  auto x1586 = x964*x967;
  auto x1587 = x331*x334;
  auto x1588 = x346*x5*x849;
  auto x1589 = x1323*x43;
  auto x1590 = x1323*x21;
  auto x1591 = x1324*x24;
  auto x1592 = x154*x34*x854;
  auto x1593 = x445*x855;
  auto x1594 = x783*x861;
  auto x1595 = x1239*x868;
  auto x1596 = x186*x887;
  auto x1597 = x871*x97;
  auto x1598 = x200*x373;
  auto x1599 = x588*x875;
  auto x1600 = x592*x876;
  auto x1601 = x592*x904;
  auto x1602 = x1295*x378;
  auto x1603 = x859*x883;
  auto x1604 = x212*x868;
  auto x1605 = x211*x887;
  auto x1606 = x120*x889;
  auto x1607 = x200*x886;
  auto x1608 = x200*x890;
  auto x1609 = 0.058322628732349752*x891;
  auto x1610 = x1304*x680;
  auto x1611 = x1281*x1587 + x1325*x481 + x1326*x911 + x1330*x435 + x1588*x473 + x1589*x913 + x1590*x915 + x1591*x464 + x1592*x479 + x1592*x480 + x1593*x452 - x1594*x840 - x1594*x843 + x1595*x777 + x1596*x777 + x1597*x449 + x1597*x453 + x1598*x919 + x1599*x919 - x1600*x783 - x1601*x443 - x1602*x840 - x1602*x843 - x1603*x449 - x1603*x453 + x1604*x921 + x1604*x922 + x1605*x922 + x1606*x457 + x1607*x449 + x1607*x453 + x1608*x449 + x1608*x453 - x1609*x457 + x1610*x885 + x1610*x887 + x440*x863 + x444*x863 - x470 - x777*x860 + x784*x892;
  auto x1612 = x1280 + x1611;
  auto x1613 = x965*x971;
  auto x1614 = (4.0/3.0)*x940;
  auto x1615 = x820*x861;
  auto x1616 = x1315*x378;
  auto x1617 = x1316*x680;
  auto x1618 = x1312*x1587 + x1325*x523 + x1326*x954 + x1330*x499 + x1588*x516 + x1589*x952 + x1590*x953 + x1591*x511 + x1592*x521 + x1592*x522 + x1593*x504 + x1595*x816 + x1596*x816 + x1597*x505 + x1597*x507 + x1598*x956 + x1599*x956 - x1600*x820 - x1601*x502 - x1603*x505 - x1603*x507 + x1604*x959 + x1604*x960 + x1605*x960 + x1606*x509 + x1607*x505 + x1607*x507 + x1608*x505 + x1608*x507 - x1609*x509 - x1615*x840 - x1615*x843 - x1616*x840 - x1616*x843 + x1617*x885 + x1617*x887 + x501*x863 + x503*x863 - x515 - x816*x860 + x823*x892;
  auto x1619 = x1360*x487;
  auto x1620 = x526*x943;
  auto x1621 = x1546 - x1619*x2 + x1620*x2 - 2*x946;
  auto x1622 = x1360*x788;
  auto x1623 = x824*x943;
  auto x1624 = x1393 + x1395 - x1622*x2 + x1623*x2 + x21*x941 - x21*x944;
  auto x1625 = x1360*x926;
  auto x1626 = x1360*x971;
  auto x1627 = x943*x963;
  auto x1628 = x943*x967;
  auto x1629 = x8/x101;
  auto x1630 = (1.0/3.0)*x1629*x520;
  auto x1631 = std::pow(x939, 2);
  auto x1632 = (-x1631 - 1)/x1631;
  auto x1633 = x1632*x825;
  auto x1634 = x1467 + x1611;
  auto x1635 = -x1625*x2 - x1626*x2 + x1627*x2 + x1628*x2 - x1630*x2 - x1633*x2 + x1634;
  auto x1636 = x1065*x1360;
  auto x1637 = x1056*x943;
  auto x1638 = x34*x707;
  auto x1639 = x183*x331;
  auto x1640 = x341*x43;
  auto x1641 = x21*x341;
  auto x1642 = x13*x343;
  auto x1643 = x346*x914;
  auto x1644 = x166*x368;
  auto x1645 = x74*x990;
  auto x1646 = x357*x358;
  auto x1647 = x186*x375;
  auto x1648 = x366*x445;
  auto x1649 = x202*x920;
  auto x1650 = x372*x74;
  auto x1651 = x372*x378;
  auto x1652 = x1253*x211;
  auto x1653 = x1569*x211;
  auto x1654 = -x1002*x1644 + x1002*x1652 + x1006*x1651 + x1006*x1653 + x1011*x380 + x1011*x725 + x1014*x1648 + x1015*x380 + x1015*x725 - x1018*x1649 + x1481*x1650 + x1484*x924 + x1484*x925 + x1638*x982 + x1639*x986 + x1640*x995 + x1641*x999 + x1642*x994 + x1643*x998 + x1645*x1646 + x1645*x1647 + x709*x987 + x910*x982;
  auto x1655 = x1654 + x753*x976;
  auto x1656 = x1020*x1360;
  auto x1657 = x1047*x943;
  auto x1658 = -x1656*x2 + x1657*x2;
  auto x1659 = x1320*x979 - x1636*x2 + x1637*x2 + x1655 + x1658;
  auto x1660 = x1047*x1360;
  auto x1661 = x1020*x943;
  auto x1662 = x1032*x74;
  auto x1663 = x1029*x1638 + x1029*x910 + x1030*x1639 + x1031*x709 + x1034*x1642 + x1035*x1640 + x1036*x1643 + x1037*x1641 - x1038*x1644 + x1038*x1652 + x1040*x1651 + x1040*x1653 + x1042*x1648 + x1043*x380 + x1043*x725 + x1045*x380 + x1045*x725 - x1046*x1649 + x1491*x1650 + x1492*x924 + x1492*x925 + x1646*x1662 + x1647*x1662;
  auto x1664 = -x1056*x966 - x1065*x965 + x1663;
  auto x1665 = x1051*x753;
  auto x1666 = x1049*x894;
  auto x1667 = x1053*x1320;
  auto x1668 = x1667*x2;
  auto x1669 = -x1665*x8 - x1666*x2 + x1668*x8;
  auto x1670 = -x1660*x2 + x1661*x2 + x1664 + x1669;
  auto x1671 = -x1665;
  auto x1672 = x1668 + x1671;
  auto x1673 = -x1636 + x1637 + x1654;
  auto x1674 = x1658 + x1672 + x1673;
  auto x1675 = x1548 - x487*x965 - x526*x966 - 2*x968;
  auto x1676 = x1457 - x788*x965 - x824*x966;
  auto x1677 = x1512 - x1613 + x1618 - 2.0/3.0*x1632*x520 - x478*((2.0/3.0)*x1629 + (1.0/6.0)*x8) - x926*x965 - x963*x966 - x966*x967;
  auto x1678 = -x1020*x965 - x1047*x966;
  auto x1679 = x1076*x940;
  auto x1680 = (1.0/3.0)*x1679;
  auto x1681 = x1072*x943;
  auto x1682 = -x1680*x2 + x1681*x2;
  auto x1683 = x1663 + x1669 + x1678 + x1682;
  auto x1684 = x1076*x965;
  auto x1685 = -x1020*x966 - x1047*x965 + x1069*x894 + x1071*x1320 - x1072*x966 + x1654 - x1684;
  auto x1686 = x1663 + x1681;
  auto x1687 = -x1666 + x1667*x8;
  auto x1688 = x1678 - x1680 + x1686 + x1687;
  auto x1689 = x1546 - x1619 + x1620 - 2*x972;
  auto x1690 = x1393 + x1461 - x1622 + x1623;
  auto x1691 = -x1625 - x1626 + x1627 + x1628 - x1630 - x1633 + x1634;
  auto x1692 = x1654 - x1656 + x1657;
  auto x1693 = x1080*x1360;
  auto x1694 = x1076*x943;
  auto x1695 = -x1693*x2 + x1694*x2;
  auto x1696 = x1672 + x1692 + x1695;
  auto x1697 = -x1080*x965;
  auto x1698 = -x1660 + x1661 + x1663 - 2.0/3.0*x1679 + x1687 + x1697;
  auto x1699 = x1667 + x1692 - x1693 + x1694;
  auto x1700 = 0.067345170796937498*x38;
  auto x1701 = 0.26082672494403297*x136;
  auto x1702 = 0.075294189929142996*x127;
  auto x1703 = 0.45176513957485664*x127;
  auto x1704 = 0.075294189929142996*x572;
  auto x1705 = 0.075294189929142996*x154;
  auto x1706 = 0.029161314366174917*x144*x557;
  auto x1707 = x562*x58;
  auto x1708 = 0.019440876244116584*x556;
  auto x1709 = x151*x66;
  auto x1710 = 0.037647094964571498*x556;
  auto x1711 = 0.11664525746469967*x144;
  auto x1712 = 0.077763504976466336*x151;
  auto x1713 = 0.15058837985828599*x154;
  auto x1714 = x1713*x174;
  auto x1715 = 0.067345170796937498*x199;
  auto x1716 = 0.013746775424423291*x165;
  auto x1717 = x595*x596;
  auto x1718 = x1716*x556;
  auto x1719 = 0.016836292699234374*x199;
  auto x1720 = x1719*x556;
  auto x1721 = 0.058322628732349835*x144;
  auto x1722 = 0.038881752488233168*x151;
  auto x1723 = 0.077763504976466336*x58;
  auto x1724 = x187*x582;
  auto x1725 = x1723*x584;
  auto x1726 = 0.075294189929142996*x577;
  auto x1727 = x1715*x208;
  auto x1728 = 0.134690341593875*x199;
  auto x1729 = x203*x613;
  auto x1730 = 0.0097204381220582919*x615;
  auto x1731 = 0.054987101697693164*x606;
  auto x1732 = 0.054987101697693164*x165;
  auto x1733 = x206*x609;
  auto x1734 = 0.054987101697693164*x573;
  auto x1735 = x1734*x618;
  auto x1736 = x1715*x608;
  auto x1737 = x1715*x611;
  auto x1738 = 0.038881752488233168*x203;
  auto x1739 = 0.1166452574646995*x130;
  auto x1740 = x1721*x557;
  auto x1741 = 0.038881752488233168*x564;
  auto x1742 = x1722*x66;
  auto x1743 = 0.027493550848846582*x165;
  auto x1744 = x1743*x92;
  auto x1745 = x1743*x99;
  auto x1746 = 0.033672585398468749*x199;
  auto x1747 = x108*x1746;
  auto x1748 = x116*x1746;
  auto x1749 = x17*x1739 + x1740 + x1741 + x1742 + x1744 + x1745 + x1747 + x1748 + 0.019440876244116584*x615;
  auto x1750 = x122*x21;
  auto x1751 = 0.075294189929142996*x174;
  auto x1752 = x1746*x224;
  auto x1753 = x1715*x217;
  auto x1754 = x1053*x488;
  auto x1755 = 3.375*x423;
  auto x1756 = 3.375*x425;
  auto x1757 = x183*x34;
  auto x1758 = 0.125*x78;
  auto x1759 = p[15]*x1758;
  auto x1760 = 0.125*x76;
  auto x1761 = p[25]*x1760;
  auto x1762 = x186*x97;
  auto x1763 = p[16]*x1758;
  auto x1764 = p[26]*x1760;
  auto x1765 = x211*x439;
  auto x1766 = p[18]*x1758;
  auto x1767 = p[28]*x1760;
  auto x1768 = x211*x98;
  auto x1769 = p[17]*x1758;
  auto x1770 = p[27]*x1760;
  auto x1771 = x103*x445;
  auto x1772 = p[21]*x1758;
  auto x1773 = p[31]*x1760;
  auto x1774 = x111*x445;
  auto x1775 = p[19]*x1758;
  auto x1776 = p[29]*x1760;
  auto x1777 = x1019*(x1775 + x1776) + x1757*(x1755 + x1756) + x1762*(x1759 + x1761) + x1765*(x1763 + x1764) + x1768*(x1766 + x1767) + x1771*(x1769 + x1770) + x1774*(-x1772 - x1773) + x519;
  auto x1778 = x1777 + x487*x979 + x788*x976;
  auto x1779 = (2.0/3.0)*x21;
  auto x1780 = x1779*x804;
  auto x1781 = (2.0/3.0)*x146;
  auto x1782 = f[0]*x488;
  auto x1783 = (2.0/9.0)*x1782;
  auto x1784 = x1783*x2;
  auto x1785 = x1019*(-x1775 - x1776) + x1757*(-x1755 - x1756) + x1762*(-x1759 - x1761) + x1765*(-x1763 - x1764) + x1768*(-x1766 - x1767) + x1771*(-x1769 - x1770) + x1774*(x1772 + x1773) + x477;
  auto x1786 = 0.074074074074073987*x1103;
  auto x1787 = 0.075990887731753332*x624;
  auto x1788 = 0.089793561062583321*x242;
  auto x1789 = x1096*x1788;
  auto x1790 = x1097*x1788;
  auto x1791 = 0.48050614670408442*x246;
  auto x1792 = x1791*x7;
  auto x1793 = x11*x1791;
  auto x1794 = 0.26157716347320858*x1100;
  auto x1795 = 0.25595878053599447*x1105;
  auto x1796 = 0.57901139454755335*x1099;
  auto x1797 = 0.11664525746469967*x657;
  auto x1798 = 0.077763504976466336*x645;
  auto x1799 = 0.2509806330971433*x652;
  auto x1800 = x1799*x49;
  auto x1801 = x1799*x57;
  auto x1802 = 0.11664525746469949*x655;
  auto x1803 = 0.12830323729461737*x237;
  auto x1804 = x1803*x92;
  auto x1805 = x1803*x99;
  auto x1806 = 0.25869917107112544*x259;
  auto x1807 = x108*x1806;
  auto x1808 = x116*x1806;
  auto x1809 = 0.12821185290049902*x646;
  auto x1810 = -x1786 - x1787 - x1789 - x1790 - x1792 - x1793 - x1794 - x1795 - x1796 - x1797 - x1798 - x1800 - x1801 - x1802 - x1804 - x1805 - x1807 - x1808 - x1809;
  auto x1811 = x1810*x2;
  auto x1812 = x1049*x411;
  auto x1813 = x1051*x280;
  auto x1814 = -2*x1057 + x1063 - x1811*x8 - x1812*x2 - x1813*x8;
  auto x1815 = 0.134690341593875*x38;
  auto x1816 = x22*x315;
  auto x1817 = 0.23329051492939901*x130;
  auto x1818 = 0.52165344988806595*x136;
  auto x1819 = x10*x1818;
  auto x1820 = 0.33333333333333331*x127;
  auto x1821 = 0.90353027914971595*x127;
  auto x1822 = 0.15058837985828599*x331;
  auto x1823 = x1713*x34;
  auto x1824 = 0.15552700995293267*x58;
  auto x1825 = x376*x877;
  auto x1826 = x1566*x375;
  auto x1827 = x1715*x373;
  auto x1828 = 0.077763504976466336*x203;
  auto x1829 = -0.32992261018615898*x1204 - 0.23329051492939901*x1211 + x1213*x1711 + x1214*x1723 + x1218*x1713 + 0.054987101697693164*x1247 + 0.11664525746469967*x1551 + 0.077763504976466336*x1552 + x1712*x339 + x1715*x420 + x1728*x367 + x1728*x383 + x1728*x408 + x1734*x384 + x1738*x385 - x1815*x1816 - x1817*x321 - x1819*x323 - x1820*x319 - x1821*x327 + x1822*x569 + x1822*x571 + x1823*x397 + x1824*x402 + 0.10997420339538633*x1825 + 0.10997420339538633*x1826 + x1827*x224 + x1828*x409;
  auto x1830 = 0.15058837985828599*x127;
  auto x1831 = 0.90353027914971329*x127;
  auto x1832 = 0.15058837985828599*x174;
  auto x1833 = 0.10997420339538633*x165;
  auto x1834 = 0.10997420339538633*x573;
  auto x1835 = x1728*x217;
  auto x1836 = -0.32992261018615898*x1113 - 0.23329051492939901*x1122 + x1136*x1711 + 0.054987101697693164*x1183 - x126*x1830 + 0.077763504976466336*x1329 - x134*x1817 + 0.038881752488233168*x1353 - x139*x1818 - x142*x1831 + x148*x1711 + x150*x1712 + x160*x1713 + x1712*x296 + x1713*x255 + x1715*x307 + x1727*x224 + x1728*x198 + x1732*x207 - x1750*x1815 + x1824*x267 + x1828*x233 + x1832*x569 + x1832*x571 + x1833*x272 + x1834*x232 + x1835*x223 + x1835*x224;
  auto x1837 = x1836*x2;
  auto x1838 = x154*x1822;
  auto x1839 = x1711*x262;
  auto x1840 = x1712*x295;
  auto x1841 = 0.15058837985828599*x706;
  auto x1842 = 0.038881752488233168*x58;
  auto x1843 = x1728*x586;
  auto x1844 = 0.054987101697693164*x1566;
  auto x1845 = 0.10997420339538633*x742;
  auto x1846 = x1715*x738;
  auto x1847 = x1728*x743;
  auto x1848 = x1246*x1845 + x1392*x1828 + x1566*x1845 + x16*x1711*x35 + x163*x1838 + x1711*x700 + x1711*x758 + x1712*x701 + x1712*x702 + x1712*x757 + x1714*x396 + x1714*x397 + x1715*x749 + x1715*x762 + x1723*x356*x759 + x1725*x375 + x1727*x724 + x1732*x736 + x1732*x739 + x1735*x375 + x1738*x751 + x1738*x752 + x1743*x731 + x1746*x744 + x1752*x733 + x1817*x695 + x1821*x754 + x1823*x704 + x1823*x705 + x1824*x187*x723 + x1827*x586 + x1831*x755 + x1835*x724 + x1838*x551 + x1839*x706 + x1840*x706 + x1841*x569 + x1841*x571 + x1842*x718 + x1843*x381 + x1844*x737 + x1846*x223 + x1846*x224 + x1847*x223 + x1847*x224 + 0.69987154478819691*x203*x684 + 0.019440876244116584*x203*x750 - x205*x686 - x209*x712 - x211*x689 - x211*x715 - x228*x691 - x228*x727 - x231*x694 - x231*x729 + 0.027493550848846582*x573*x734;
  auto x1849 = 0.23329051492939934*x144;
  auto x1850 = 0.15552700995293267*x151;
  auto x1851 = 0.30117675971657198*x154;
  auto x1852 = x174*x1851;
  auto x1853 = 0.075294189929142996*x34;
  auto x1854 = x1853*x569;
  auto x1855 = 0.65984522037231796*x573;
  auto x1856 = 0.15058837985828599*x577;
  auto x1857 = 0.26938068318774999*x199;
  auto x1858 = 0.46658102985879801*x203;
  auto x1859 = 0.10997420339538633*x606;
  auto x1860 = x1728*x608;
  auto x1861 = x1728*x611;
  auto x1862 = -x108*x1715 - x116*x1715 - x17*x1817 - x1711*x557 - x1712*x66 - x1732*x92 - x1732*x99 - 0.077763504976466336*x564 - 0.038881752488233168*x615;
  auto x1863 = x21*(-x1097*x1815 + x1246*x1859 + x1305*x1828 + x1566*x1859 + x163*x1852 - x1704*x556 - 0.038881752488233168*x1707 - x1708*x615 + x1713*x546 - 0.027493550848846582*x1717 + x1724*x1824 + x1728*x555 - 0.019440876244116584*x1729 + x1733*x1833 - x1740*x556 - x1741*x556 - x1742*x556 - x1743*x594 - x1744*x556 - x1745*x556 - x1746*x603 - x1746*x604 - x1747*x556 - x1748*x556 + x1818*x541 - x1818*x542 + x1824*x216*x584 + x1828*x619 - x1830*x20 + x1831*x543 - x1831*x544 + x1834*x216*x618 + x1839*x577 + x1840*x577 + x1843*x208 + x1849*x548 + x1850*x550 + x1852*x551 - x1854*x556 - x1855*x534 - x1855*x538 + x1856*x569 + x1856*x571 + x1857*x664 - x1858*x539 - x1858*x600 + x1860*x223 + x1860*x224 + x1861*x223 + x1861*x224 + x1862 - 0.15058837985828599*x572);
  auto x1864 = x1053*x1109;
  auto x1865 = x1864*x2;
  auto x1866 = -x1054*x132 - x146*x1848 + x1494 + x1520 + x1829*x23 - x1837*x8 - x1863*x8 + x1865*x8;
  auto x1867 = x1829*x8;
  auto x1868 = 0.15058837985828599*x34;
  auto x1869 = x81*x859;
  auto x1870 = 0.65984522037231796*x877;
  auto x1871 = x876*x877;
  auto x1872 = x378*x879;
  auto x1873 = x1858*x378;
  auto x1874 = x1746*x74;
  auto x1875 = 0.10997420339538633*x885;
  auto x1876 = 0.10997420339538633*x887;
  auto x1877 = 0.019440876244116584*x853;
  auto x1878 = x1728*x886;
  auto x1879 = x1728*x890;
  auto x1880 = x203*x891;
  auto x1881 = x1828*x680;
  auto x1882 = x22*(x101*x1819 - x1096*x1815 + x1246*x1875 + x1246*x1876 - x15*x1820 + x1550*x1868 + x1566*x1875 + x1566*x1876 + x1704*x853 + x1728*x856 + x1728*x872 + x1740*x853 + x1741*x853 + x1742*x853 + x1744*x853 + x1745*x853 + x1747*x853 + x1748*x853 - x1819*x9 + x1821*x845 - x1821*x846 + x1824*x870 + x1824*x903 - x1842*x1869 + x1849*x847 + x1850*x848 + x1851*x851 + x1851*x852 + x1854*x853 + x1857*x900 + x1862 - x1868*x569 - x1870*x840 - x1870*x843 - 0.027493550848846582*x1871 - 0.027493550848846582*x1872 - x1873*x840 - x1873*x843 - x1874*x882 - x1874*x884 + x1877*x615 + x1878*x223 + x1878*x224 + x1879*x223 + x1879*x224 - 0.019440876244116584*x1880 + x1881*x885 + x1881*x887);
  auto x1883 = x1053*x1200;
  auto x1884 = x1883*x2;
  auto x1885 = -x1054*x146 - x132*x1848 + x1664 + x1682 + x1836*x23 - x1867*x2 - x1882*x2 + x1884*x8;
  auto x1886 = x1754*x8;
  auto x1887 = x1886*x2;
  auto x1888 = x1051*x804;
  auto x1889 = x1049*x945;
  auto x1890 = x1049*x926;
  auto x1891 = x1051*x788;
  auto x1892 = x1053*x487;
  auto x1893 = x1892*x2;
  auto x1894 = -x1890*x2 - x1891*x8 + x1893*x8;
  auto x1895 = x1785 + x1887 - x1888*x8 - x1889*x2 + x1894 + x527*x979 + x828*x976;
  auto x1896 = x1049*x963;
  auto x1897 = x1051*x824;
  auto x1898 = x1049*x967;
  auto x1899 = x1051*x828;
  auto x1900 = x1053*x526;
  auto x1901 = x1900*x2;
  auto x1902 = x1053*x527;
  auto x1903 = x1902*x2;
  auto x1904 = (1.0/18.0)*x1782;
  auto x1905 = x1070*x1904 + x1777;
  auto x1906 = x1069*x945 - x1896*x2 - x1897*x8 - x1898*x2 - x1899*x8 + x1901*x8 + x1903*x8 + x1905;
  auto x1907 = x1785 + x1886;
  auto x1908 = -x1889 + x1907;
  auto x1909 = -x1899 + x1903;
  auto x1910 = x1894 + x1908 + x1909;
  auto x1911 = x1026 + x1061 - 2*x1066 - x1811 - x1813;
  auto x1912 = -x1054*x21 + x1485 - x1488 + x1490 + x1535 - x1837 - x1863 + x1865;
  auto x1913 = x1671 + x1673 + x1695 + x1884;
  auto x1914 = x1754*x2 + x1777 - x1891 + x1893;
  auto x1915 = -x1888 + x1904*x977 + x1914 + x835*x976;
  auto x1916 = x1051*x835;
  auto x1917 = x1049*x971;
  auto x1918 = x1887 - x1916*x8 - x1917*x2;
  auto x1919 = -x1897 + x1901 + x1908 + x1918;
  auto x1920 = x1754 + x1914 - x1916;
  auto x1921 = x10*x1701;
  auto x1922 = 0.16666666666666666*x127;
  auto x1923 = 0.037647094964571498*x853;
  auto x1924 = 0.016836292699234374*x588;
  auto x1925 = x1716*x853;
  auto x1926 = x1719*x853;
  auto x1927 = 0.054987101697693164*x1246;
  auto x1928 = x1715*x886;
  auto x1929 = x1715*x890;
  auto x1930 = x1738*x680;
  auto x1931 = 0.075294189929142996*x331;
  auto x1932 = (2.0/9.0)*f[0]*x527*x8;
  auto x1933 = x1069*x926 + x1071*x487 + x1777;
  auto x1934 = (2.0/3.0)*x22;
  auto x1935 = -2*x1077 + x1078 - x1810*x8 - x1812;
  auto x1936 = x1493 + x1528 + x1530 + x1540 + x1864*x8;
  auto x1937 = -x1054*x22 - x1679 + x1686 + x1697 - x1867 - x1882 + x1883*x8;
  auto x1938 = -x1890 + x1892*x8;
  auto x1939 = x1785 + x1909 + x1918 + x1938;
  auto x1940 = x1069*x971 - x1896 - x1898 + x1900*x8 + x1902*x8 + x1905;
  auto x1941 = x1902 + x1907 - x1917 + x1938;
  auto x1942 = x1777 + x1892;

 sum[0]=6.8175478160762513*std::pow(f[0], 0.2915026221291801)*x17 + 13.415833701377133*std::pow(f[0], 1.9462219947249002)*x119*x120 + x1*x11 + x1*x7 + x100*x108 + x100*x116 + 70.465998972382408*x12*x15 + 31.834081861040136*x12*x20 + x37 + x40*x49 + x40*x57 + 9.1858140253450848*x58*x66 + x85*x92 + x85*x99 - x83/x58 + 20.5175396875734*x0/std::pow(f[0], 7);

 sum[36]=x533;

 sum[72]=-59.482526059792832*x414 + x421 + x839;

 sum[108]=x975;

 sum[144]=x1068;

 sum[180]=x975;

 sum[6]=x533;

 sum[42]=-x1084*x600 + x1085 + x1095*x1096 - x1097*x1098 + 0.86851709182132997*x1099 - 1.9618287260490646*x1100 + x1101*x49 - x1102*x57 + x1108 - 2*x641 + x682;

 sum[78]=x1319;

 sum[114]=x1463;

 sum[150]=x1505;

 sum[186]=x1532;

 sum[12]=x422 + x839;

 sum[48]=x1319;

 sum[84]=x1095*x1097 - x1096*x1098 - 4.3425854591066502*x1099 + 0.39236574520981288*x1100 + x1101*x57 - x1102*x49 + x1108 - x1544*x843 + x1545 - 2*x895 + x908;

 sum[120]=x1621;

 sum[156]=x1675;

 sum[192]=x1689;

 sum[18]=x975;

 sum[54]=x1463;

 sum[90]=x1621;

 sum[126]=-2*x1022 + x1023 + x1026 + x280*x976 - x637*(-1.0/18.0*x110 - 1.0/9.0);

 sum[162]=x1814;

 sum[198]=x1911;

 sum[24]=x1068;

 sum[60]=x1505;

 sum[96]=x1675;

 sum[132]=x1814;

 sum[168]=x1026 + x1069*x411 - 2*x1073 + x1074 - x637*(-1.0/18.0*x101 - 1.0/9.0);

 sum[204]=x1935;

 sum[30]=x975;

 sum[66]=x1532;

 sum[102]=x1689;

 sum[138]=x1911;

 sum[174]=x1935;

 sum[210]=x1026 + x1060 - 2*x1081 + x1786 + x1787 + x1789 + x1790 + x1792 + x1793 + x1794 + x1795 + x1796 + x1797 + x1798 + x1800 + x1801 + x1802 + x1804 + x1805 + x1807 + x1808 + x1809;

 sum[1]=-x169*x283 - x201*x286 - 2*x282 + x314;

 sum[37]=x1085 + x683;

 sum[73]=x1086 + x1318 + x767;

 sum[109]=x807;

 sum[145]=x834;

 sum[181]=x838;

 sum[7]=-x534*x536 - x536*x538 - x539*x540 + x683;

 sum[43]=x1110*x1111 + x1112*x139 + 0.74232587291885777*x1113*x556 + 2.4744195763961923*x1113 - x1116*x1159 - x1116 - x1118*x1159 - x1118 - x1120*x1159 - x1120 - x1121*x1159 - x1121 + 0.52490365859114774*x1122*x556 + 1.7496788619704926*x1122 + x1123*x1124*x154 - x1124*x571 + x1125 + x1126*x1174*x91*x96 + 0.17496788619704925*x1126*x178*x58 - x1127*x1128 - x1128*x1129 - x1130*x1162*x91 - 0.74232587291885777*x1130*x596 - x1131*x1133 - x1133*x1161 + x1134*x1147*x208 - x1134*x198 - x1135*x1157*x262 - x1135*x148*x556 - x1136*x1137 + x1137*x1150*x577 - x1138*x1157*x295 - x1138*x150*x556 - x1139*x267 + x1139*x584*x606 + x1140*x147*x578 - x1140*x296 + x1141 + x1142 + x1143*x163 + x1143*x551 + x1144*x1145 + x1145*x1146 - x1146*x1179*x208 + x1147*x1148*x217 + x1148*x586*x611 + x1149*x165*x88*x91 + x1149*x596*x96 + x115*x200*(-40*x138 + x196) + x1151 + x1152 - x1153*x1155 - x1153*x1175 - x1155*x1156 - x1155*x1176 - x1156*x1175 - x1158*x569 - x1158*x571 - 2.9693034916754311*x1160*x573 - x1160*x703*x90 - x1163*x198*x556 - x1163*x586*x602 - x1164*x223 - x1164*x224 - x1165*x304 - x1166*x1167*x80 - x1168*x223 - x1168*x224 + x1169*x1171 + x1169*x176 + x1170 + x1172*x188 + x1173*x49 + x1173*x57 - x1174*x230*x556 - x1175*x1176 - x1178*x217*x90 - 0.15152663429310936*x1179*x307 - x1180*x223 - x1180*x224 - x1181*x1182 - x1181*x1183 - x1185*x609 - 0.24744195763961924*x1186*x216 - x1187*x223 - x1187*x224 + x1188*x607 + x1188*x617 + x1189*x219 + x1189*x382 + x1189*x612 + x1190*x1191 + x1191*x1192 + x1193*x606*x711 + x1194*x713 + x1194*x726*x90 - x1196*x609 + x1198 + x125 + 10.84236334979656*x127*x142 + x129 + x155*x159*(8*x156 - x2) - x164 - x193*x556 - 5.9386069833508612*x2*x217*x573 - x215*x556 - x220*x556 - x226*x556 - x227*x556 - x232*x713 - x233*x726 - x234*x556 - x272*x711 + 1.2122130743448749*x589*x606;

 sum[79]=x1321;

 sum[115]=x1469;

 sum[151]=x1513;

 sum[187]=x1533;

 sum[13]=x768;

 sum[49]=x1321;

 sum[85]=x1198 + x1358;

 sum[121]=x1624;

 sum[157]=x1676;

 sum[193]=x1690;

 sum[19]=x807;

 sum[55]=x1469;

 sum[91]=x1624;

 sum[127]=-x1021*x1277 + x1054*x138 + x1109*x979 - 2*x1471 + x1485 - x157*(0.16496130509307949*x1113 + 0.1166452574646995*x1122 - x1136*x1721 - 0.027493550848846582*x1183 + x126*x1702 - 0.038881752488233168*x1329 + x134*x1739 - 0.019440876244116584*x1353 + x139*x1701 + x142*x1703 - x148*x1721 - x150*x1722 - x160*x1705 + x1700*x1750 - x1705*x255 - x1715*x198 - x1722*x296 - x1723*x267 - x1732*x272 - x1734*x232 - x1738*x233 - x1743*x207 - x1746*x307 - x1751*x569 - x1751*x571 - x1752*x208 - x1753*x223 - x1753*x224) - x42*(x108*x1720 + x1097*x1700 + x116*x1720 - x1246*x1731 - x1305*x1738 - x1566*x1731 - x163*x1714 - x1701*x541 + x1701*x542 + x1702*x20 - x1703*x543 + x1703*x544 + x1704 - x1705*x546 + x1706*x556 + 0.019440876244116584*x1707 + x1708*x1709 + x1708*x564 + x1710*x570 + x1710*x572 - x1711*x548 - x1712*x550 - x1714*x551 - x1715*x555 + x1716*x594 + 0.013746775424423291*x1717 + x1718*x92 + x1718*x99 + x1719*x603 + x1719*x604 - x1721*x262*x577 - x1722*x295*x577 - x1723*x1724 - x1725*x216 - x1726*x569 - x1726*x571 - x1727*x586 - x1728*x664 + 0.0097204381220582919*x1729 + x1730*x556 - x1732*x1733 - x1735*x216 - x1736*x223 - x1736*x224 - x1737*x223 - x1737*x224 - x1738*x619 + x1749 + x218*x534 + x218*x538 + x231*x539 + x231*x600);

 sum[163]=x1866;

 sum[199]=x1912;

 sum[25]=x834;

 sum[61]=x1513;

 sum[97]=x1676;

 sum[133]=x1866;

 sum[169]=x1056*x799*x8 + x1071*x1109 - x1076*x8*x803 + x1526;

 sum[205]=x1936;

 sum[31]=x838;

 sum[67]=x1533;

 sum[103]=x1690;

 sum[139]=x1912;

 sum[175]=x1936;

 sum[211]=x1065*x799 - x1080*x803 + x1485 + x1864;

 sum[2]=-x283*x355 - x286*x369 - 2*x413 + x422;

 sum[38]=x1090 + x765;

 sum[74]=x1545 + x909;

 sum[110]=x950;

 sum[146]=x970;

 sum[182]=x974;

 sum[8]=x768;

 sum[44]=x1202 + x1275;

 sum[80]=x1359;

 sum[116]=x1470;

 sum[152]=x1516;

 sum[188]=x1534;

 sum[14]=-x840*x842 - x840*x844 - x842*x843 + x909;

 sum[50]=x1359;

 sum[86]=x107*x200*(40*x323 + x365) + x1111*x1320 + x1112*x416 + x1134*x1333*x373 - x1134*x367 - x1134*x383 - x1134*x408 + x1135*x1213*x853 + x1138*x339*x853 + x1139*x1342*x868 - x1139*x402 + x1144*x1564 + x1146*x1564 + x1146*x1575*x373 + x1148*x1333*x381 + x1156*x1563 + x1156*x1572 - x1162*x1559*x378 + x1163*x367*x853 - x1166*x1335*x364 + x1174*x1558*x680 + x1174*x385*x853 - x1177*x1570*x378 - x1184*x1579*x370 + x1190*x1583 + x1192*x1583 - x1195*x1579*x680 + x1202 - 0.74232587291885777*x1204*x853 + 2.4744195763961923*x1204 - 0.52490365859114774*x1211*x853 + 1.7496788619704926*x1211 + x1215 + x1218*x1553 + x1220*x1553 + x1223 + x1224 + x1235 + x1236 + x1246*x1565 + x1247*x1577 + x1249*x1577 + x1250 + 10.842363349796592*x127*x327 + x1324*x399 + x1325*x398 + x1328*x1554 + x1331*x1554 + x1351*x1581 + x1549*x1550 - x1549*x569 - 1.049807317182297*x1551 - 0.69987154478819702*x1552 + x1555*x360*x864 + x1555*x866 + x1556*x569 + x1556*x571 + x1557*x569 + x1557*x571 + x1558*x58*x82 - 0.74232587291885777*x1559*x877 + x1561*x1562 + x1561*x1571 + x1562*x1573 + x1565*x1566 - 0.24744195763961924*x1566*x1579 + x1567*x223 + x1567*x224 - x1568*x378*x703 - 2.9693034916754311*x1568*x877 - x1569*x713 - 5.9386069833508612*x1570*x877 + x1571*x1573 + x1574*x188 + 0.15152663429310936*x1575*x420 + x1576*x885 + x1576*x887 + x1578*x223 + x1578*x224 - x1580*x223 - x1580*x224 + x1581*x607 + x1581*x617 + x1582*x370*x711 + x1582*x378*x713 + x1582*x680*x726 + x1584*x382 + x1584*x612 + x317 + x320 - x345 + x346*x55*(-8*x347 - x8) - x405*x711 - x409*x726 + x419*x58*x853 - 0.30305326858621873*x859*x873;

 sum[122]=x1635;

 sum[158]=x1677;

 sum[194]=x1691;

 sum[20]=x950;

 sum[56]=x1470;

 sum[92]=x1635;

 sum[128]=x1056*x2*x964 - x1065*x2*x966 + x1200*x979 + x1655;

 sum[164]=x1885;

 sum[200]=x1913;

 sum[26]=x970;

 sum[62]=x1516;

 sum[98]=x1677;

 sum[134]=x1885;

 sum[170]=x1054*x323 + x1071*x1200 - x1072*x1614 + x1654 - 2*x1684 - x348*(0.16496130509307949*x1204 + 0.1166452574646995*x1211 - x1213*x1721 - x1214*x1842 - x1218*x1705 - x1220*x1705 - 0.027493550848846582*x1247 - 0.027493550848846582*x1249 - 0.019440876244116584*x1271 + x128*x327 - 0.058322628732349835*x1551 - 0.038881752488233168*x1552 + x1700*x1816 - x1715*x367 - x1715*x383 - x1715*x408 - x1722*x339 - x1723*x402 - x1738*x409 + x1739*x321 - x1746*x420 - x1752*x373 - 0.054987101697693164*x1825 - 0.054987101697693164*x1826 + x1921*x323 + x1922*x319 - x1931*x569 - x1931*x571) - x51*(-x101*x1921 - x108*x1926 + x1096*x1700 - x116*x1926 - x128*x845 + x128*x846 + x15*x1922 - x1550*x1853 - x1706*x853 - x1709*x1877 - x1711*x847 - x1712*x848 - x1713*x851 - x1713*x852 - x1715*x856 - x1715*x872 - x1723*x870 - x1723*x903 - x1728*x900 - x1730*x853 + x1749 - x1844*x885 - x1844*x887 + x1854 + 0.019440876244116584*x1869*x58 + 0.013746775424423291*x1871 + 0.013746775424423291*x1872 - x1877*x564 + 0.0097204381220582919*x1880 + x1921*x9 - x1923*x570 - x1923*x572 + x1924*x882 + x1924*x884 - x1925*x92 - x1925*x99 - x1927*x885 - x1927*x887 - x1928*x223 - x1928*x224 - x1929*x223 - x1929*x224 - x1930*x885 - x1930*x887 + x374*x840 + x374*x843 + x386*x840 + x386*x843);

 sum[206]=x1937;

 sum[32]=x974;

 sum[68]=x1534;

 sum[104]=x1691;

 sum[140]=x1913;

 sum[176]=x1937;

 sum[212]=x1076*x964 - x1080*x966 + x1654 + x1883;

 sum[3]=x494;

 sum[39]=x1093 + x806;

 sum[75]=x1547 + x949;

 sum[111]=x1028;

 sum[147]=x1064;

 sum[183]=x1067;

 sum[9]=x807;

 sum[45]=-2*x1276 - x1277*x804 + x1307;

 sum[81]=x1396;

 sum[117]=x1487;

 sum[153]=x1521;

 sum[189]=x1537;

 sum[15]=x950;

 sum[51]=x1396;

 sum[87]=-x1585*x2 + x1586*x2 + x1612;

 sum[123]=x1659;

 sum[159]=x1683;

 sum[195]=x1696;

 sum[21]=x1028;

 sum[57]=x1487;

 sum[93]=x1659;

 sum[129]=x1754*x977 + x1778 + (1.0/3.0)*x42*x804;

 sum[165]=x1895;

 sum[201]=x1915;

 sum[27]=x1064;

 sum[63]=x1521;

 sum[99]=x1683;

 sum[135]=x1895;

 sum[171]=-x1779*x8*x828 - x1781*x967 + x1932*x2 + x1933;

 sum[207]=x1939;

 sum[33]=x1067;

 sum[69]=x1537;

 sum[105]=x1696;

 sum[141]=x1915;

 sum[177]=x1939;

 sum[213]=-x1779*x835 + x1784 + x1942;

 sum[4]=-2*x528 + x529;

 sum[40]=x1094 + x833;

 sum[76]=x1548 + x969;

 sum[112]=x1064;

 sum[148]=x1075;

 sum[184]=x1079;

 sum[10]=x834;

 sum[46]=-x1308*x8 + x1309*x8 + x1311 + x1317;

 sum[82]=x1459;

 sum[118]=x1500;

 sum[154]=x1527;

 sum[190]=x1541;

 sum[16]=x970;

 sum[52]=x1459;

 sum[88]=x1311 - 2*x1613 - x1614*x967 + x1618;

 sum[124]=x1670;

 sum[160]=x1685;

 sum[196]=x1698;

 sum[22]=x1064;

 sum[58]=x1500;

 sum[94]=x1670;

 sum[130]=-x1780*x8 - x1781*x945 + x1784*x8 + x1785 + x526*x979 + x824*x976;

 sum[166]=x1906;

 sum[202]=x1919;

 sum[28]=x1075;

 sum[64]=x1527;

 sum[100]=x1685;

 sum[136]=x1906;

 sum[172]=x1069*x963 + x1070*x1902 + x1071*x526 + x1785 + (1.0/3.0)*x51*x967;

 sum[208]=x1940;

 sum[34]=x1079;

 sum[70]=x1541;

 sum[106]=x1698;

 sum[142]=x1919;

 sum[178]=x1940;

 sum[214]=x1783*x8 + x1785 + x1900 - x1934*x971;

 sum[5]=x494;

 sum[41]=x1093 + x837;

 sum[77]=x1547 + x973;

 sum[113]=x1067;

 sum[149]=x1079;

 sum[185]=x1082;

 sum[11]=x838;

 sum[47]=x1307 - x1308 + x1309;

 sum[83]=x1462;

 sum[119]=x1502;

 sum[155]=x1531;

 sum[191]=x1542;

 sum[17]=x974;

 sum[53]=x1462;

 sum[89]=-x1585 + x1586 + x1612;

 sum[125]=x1674;

 sum[161]=x1688;

 sum[197]=x1699;

 sum[23]=x1067;

 sum[59]=x1502;

 sum[95]=x1674;

 sum[131]=x1778 - x1780 + x1784;

 sum[167]=x1910;

 sum[203]=x1920;

 sum[29]=x1079;

 sum[65]=x1531;

 sum[101]=x1688;

 sum[137]=x1910;

 sum[173]=x1932 + x1933 - x1934*x967;

 sum[209]=x1941;

 sum[35]=x1082;

 sum[71]=x1542;

 sum[107]=x1699;

 sum[143]=x1920;

 sum[179]=x1941;

 sum[215]=x1783 + x1942;
       
        return sum;
	}
    
    int getnF()
    {
        return nF;
    }
    
    int getnP()
    {
        return nP;
    }

};
#endif