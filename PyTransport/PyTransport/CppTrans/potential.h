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
// Potential file rewriten at Fri Feb  7 19:17:22 2020

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
nP=9;

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
  auto x3 = std::pow(f[0], 2)*x2;
  auto x4 = p[2]/std::pow(M_PI, 3.0/2.0);
  auto x5 = (9.0/4.0)*x3*x4;
  auto x6 = (1.0/2.0)*f[3] - 1.0/2.0*f[4] + (1.0/2.0)*f[5];
  sum=(3.0/2.0)*std::sqrt(3)*std::pow(f[0], 3.0/2.0)*x2*x4*(p[7]*std::cos(x6) - p[8]*std::sin(x6))*std::sqrt(std::sin(f[1]))*std::sqrt(std::sin(f[2]))*std::sqrt(std::tan((1.0/2.0)*f[2]))/std::sqrt(std::tan((1.0/2.0)*f[1])) + p[5]*x5*std::cos(f[2]) + p[6]*x5*std::cos(f[1]) + x1*(1 - 27.0/32.0*x0/(std::pow(M_PI, 2)*std::pow(f[0], 4)*std::pow(p[3], 4))) + (1.0/3.0)*x3;
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot
  auto x0 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x1 = f[0]*x0;
  auto x2 = std::cos(f[2]);
  auto x3 = p[2]/std::pow(M_PI, 3.0/2.0);
  auto x4 = p[5]*x3;
  auto x5 = (9.0/2.0)*x1;
  auto x6 = std::cos(f[1]);
  auto x7 = p[6]*x3;
  auto x8 = (1.0/2.0)*f[3] - 1.0/2.0*f[4] + (1.0/2.0)*f[5];
  auto x9 = std::cos(x8);
  auto x10 = std::sin(x8);
  auto x11 = p[7]*x9 - p[8]*x10;
  auto x12 = (9.0/4.0)*x0;
  auto x13 = std::sqrt(3);
  auto x14 = std::sin(f[1]);
  auto x15 = std::sqrt(x14);
  auto x16 = std::sin(f[2]);
  auto x17 = std::sqrt(x16);
  auto x18 = std::tan((1.0/2.0)*f[1]);
  auto x19 = std::pow(x18, -1.0/2.0);
  auto x20 = std::tan((1.0/2.0)*f[2]);
  auto x21 = std::sqrt(x20);
  auto x22 = x13*x15*x17*x19*x21*x3;
  auto x23 = std::pow(f[0], 2)*x12;
  auto x24 = x11*x13*x17*x21*x3;
  auto x25 = std::pow(f[0], 3.0/2.0)*x0;
  auto x26 = (3.0/2.0)*x25;
  auto x27 = x15*x26;
  auto x28 = (3.0/4.0)*x19*x25;
  auto x29 = x11*x13*x3;
  auto x30 = (1.0/2.0)*p[7]*x10;
  auto x31 = (1.0/2.0)*p[8]*x9;
  auto x32 = x22*x26;
  auto x33 = x32*(-x30 - x31);

 sum[0]=std::sqrt(f[0])*x11*x12*x22 + (2.0/3.0)*x1 + x2*x4*x5 + x5*x6*x7 + (27.0/4.0)*std::pow(p[0], 2)*std::pow(p[1], 8)/(std::pow(M_PI, 2)*std::pow(f[0], 5)*std::pow(p[3], 4));

 sum[1]=-x14*x23*x7 + x24*x27*(-1.0/4.0*std::pow(x18, 2) - 1.0/4.0)/std::pow(x18, 3.0/2.0) + x24*x28*x6/x15;

 sum[2]=x15*x2*x21*x28*x29/x17 - x16*x23*x4 + x17*x19*x27*x29*((1.0/4.0)*std::pow(x20, 2) + 1.0/4.0)/x21;

 sum[3]=x33;

 sum[4]=x32*(x30 + x31);

 sum[5]=x33;
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot
  auto x0 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x1 = (2.0/3.0)*x0;
  auto x2 = std::pow(p[0], 2)*std::pow(p[1], 8)/(std::pow(M_PI, 2)*std::pow(p[3], 4));
  auto x3 = std::cos(f[2]);
  auto x4 = p[2]*x0/std::pow(M_PI, 3.0/2.0);
  auto x5 = (9.0/2.0)*x4;
  auto x6 = p[5]*x5;
  auto x7 = x3*x6;
  auto x8 = std::cos(f[1]);
  auto x9 = p[6]*x5;
  auto x10 = x8*x9;
  auto x11 = std::sqrt(f[0]);
  auto x12 = std::sqrt(3);
  auto x13 = std::sin(f[1]);
  auto x14 = std::sqrt(x13);
  auto x15 = std::sin(f[2]);
  auto x16 = std::sqrt(x15);
  auto x17 = std::tan((1.0/2.0)*f[1]);
  auto x18 = std::pow(x17, -1.0/2.0);
  auto x19 = std::tan((1.0/2.0)*f[2]);
  auto x20 = std::sqrt(x19);
  auto x21 = x12*x14*x16*x18*x20;
  auto x22 = (1.0/2.0)*f[3] - 1.0/2.0*f[4] + (1.0/2.0)*f[5];
  auto x23 = std::cos(x22);
  auto x24 = p[7]*x23;
  auto x25 = std::sin(x22);
  auto x26 = p[8]*x25;
  auto x27 = x24 - x26;
  auto x28 = x27*x4;
  auto x29 = x21*x28;
  auto x30 = 1.0/f[0];
  auto x31 = (9.0/4.0)*x4;
  auto x32 = std::pow(f[0], 2)*x31;
  auto x33 = p[6]*x32;
  auto x34 = std::pow(f[0], 3.0/2.0);
  auto x35 = x28*x34;
  auto x36 = (3.0/2.0)*x35;
  auto x37 = std::pow(x17, -3.0/2.0);
  auto x38 = std::pow(x17, 2);
  auto x39 = -1.0/4.0*x38 - 1.0/4.0;
  auto x40 = x12*x20;
  auto x41 = x16*x40;
  auto x42 = x14*x37*x39*x41;
  auto x43 = (3.0/4.0)*x34;
  auto x44 = x28*x43;
  auto x45 = 1.0/x14;
  auto x46 = x18*x45;
  auto x47 = x41*x46;
  auto x48 = x47*x8;
  auto x49 = -x13*x33 + x36*x42 + x44*x48;
  auto x50 = x11*x31;
  auto x51 = x27*x50;
  auto x52 = (9.0/8.0)*x11*x28;
  auto x53 = -f[0]*x13*x9 - x30*x49 + x42*x51 + x48*x52;
  auto x54 = p[5]*x32;
  auto x55 = x14*x18;
  auto x56 = 1.0/x20;
  auto x57 = std::pow(x19, 2);
  auto x58 = (1.0/4.0)*x57;
  auto x59 = x58 + 1.0/4.0;
  auto x60 = x12*x16*x56*x59;
  auto x61 = x55*x60;
  auto x62 = 1.0/x16;
  auto x63 = x40*x62;
  auto x64 = x3*x63;
  auto x65 = x55*x64;
  auto x66 = -x15*x54 + x36*x61 + x44*x65;
  auto x67 = -f[0]*x15*x6 - x30*x66 + x51*x61 + x52*x65;
  auto x68 = (1.0/2.0)*p[7]*x25;
  auto x69 = (1.0/2.0)*p[8]*x23;
  auto x70 = -x68 - x69;
  auto x71 = x21*x70;
  auto x72 = x50*x71;
  auto x73 = std::pow(x13, 2);
  auto x74 = std::pow(x8, 2);
  auto x75 = (1.0/6.0)*x73 + (1.0/9.0)*x74;
  auto x76 = 1.0/x73;
  auto x77 = x30*x76;
  auto x78 = 6*x75*x77;
  auto x79 = (2.0/3.0)*x74*x77;
  auto x80 = -x79;
  auto x81 = x34*x4;
  auto x82 = (3.0/2.0)*x81;
  auto x83 = x71*x82;
  auto x84 = std::pow(x15, 2);
  auto x85 = 1.0/x84;
  auto x86 = 2*x76 + 2*x85 - 1;
  auto x87 = (1.0/3.0)*x30*x86;
  auto x88 = std::pow(x3, 2);
  auto x89 = x30*x85;
  auto x90 = (2.0/3.0)*x88*x89;
  auto x91 = x72 - x83*(x78 + x80) - x83*(-x78*x8 + x8*x87 - x8*x90);
  auto x92 = x68 + x69;
  auto x93 = x21*x50;
  auto x94 = (1.0/6.0)*x84 + (1.0/9.0)*x88;
  auto x95 = 6*x89*x94;
  auto x96 = -x90;
  auto x97 = x21*x82;
  auto x98 = -x83*(-x3*x79 + x3*x87 - x3*x95) + x92*x93 - x92*x97*(x95 + x96);
  auto x99 = x72 - x83*(x80 + x87 + x96);
  auto x100 = (3.0/8.0)*x35;
  auto x101 = (3.0/8.0)*x29*x34;
  auto x102 = x14*x39;
  auto x103 = x36*x41;
  auto x104 = f[0]*(f[0]*x1 + f[0]*x10 + f[0]*x7 + x27*x93 + (27.0/4.0)*x2/std::pow(f[0], 5));
  auto x105 = (1.0/6.0)*x104 - x29*x43;
  auto x106 = x102*x37;
  auto x107 = x46*x8;
  auto x108 = x100*x107*x64 + x106*x36*x60 + x106*x44*x64 + x107*x44*x60;
  auto x109 = (1.0/3.0)*x74;
  auto x110 = (1.0/6.0)*x86;
  auto x111 = (1.0/3.0)*x88;
  auto x112 = (1.0/4.0)*x81;
  auto x113 = x112*x70;
  auto x114 = x113*x48;
  auto x115 = x42*x82;
  auto x116 = x115*x70;
  auto x117 = -x114 + x116 - x83*(-x109/x13 - x110*x13 + x111*x13*x85);
  auto x118 = (1.0/2.0)*x81;
  auto x119 = x118*x70;
  auto x120 = x119*x47;
  auto x121 = x120*x3;
  auto x122 = x4*x43;
  auto x123 = x115*x92 - x121*x8 + x121 + x122*x48*x92;
  auto x124 = x114 + x116 + x120;
  auto x125 = x12*x36*x55*x59;
  auto x126 = x55*x92;
  auto x127 = x126*x63;
  auto x128 = x118*x127;
  auto x129 = x61*x70*x82;
  auto x130 = -x119*x65*x8 + x122*x65*x70 + x128*x8 + x129;
  auto x131 = -x112*x127*x3 + x126*x60*x82 - x83*(x109*x15*x76 - x110*x15 - x111/x15);
  auto x132 = x113*x65 + x128 + x129;
  auto x133 = (1.0/3.0)*x13*x49;
  auto x134 = (1.0/4.0)*x24;
  auto x135 = (1.0/4.0)*x26;
  auto x136 = x97*(-x134 + x135);
  auto x137 = (1.0/3.0)*x15*x66;
  auto x138 = (1.0/9.0)*x104;
  auto x139 = x138*x8;
  auto x140 = x97*(x134 - x135);
  auto x141 = -x133*x3 - x137*x8 + x139*x3 + x140;
  auto x142 = -x133 + x136 + x139;
  auto x143 = -x137 + x138*x3 + x140;

 sum[0]=x1 + x10 + x7 + (9.0/8.0)*x29/x11 - 135.0/4.0*x2/std::pow(f[0], 6);

 sum[6]=x53;

 sum[12]=x67;

 sum[18]=x91;

 sum[24]=x98;

 sum[30]=x99;

 sum[1]=x53;

 sum[7]=-x100*x18*x41*x74/std::pow(x13, 3.0/2.0) - x101*(x38 + 1) + x102*x103*(-3.0/4.0*x38 - 3.0/4.0)/std::pow(x17, 5.0/2.0) + x103*x37*x39*x45*x8 + x105 - x33*x8;

 sum[13]=x108;

 sum[19]=x117;

 sum[25]=x123;

 sum[31]=x124;

 sum[2]=x67;

 sum[8]=x108;

 sum[14]=-x100*x40*x55*x88/std::pow(x15, 3.0/2.0) + x101*(x57 + 1) + x105 + x125*x16*(-x58 - 1.0/4.0)/std::pow(x19, 3.0/2.0) + x125*x3*x56*x62 - x3*x54;

 sum[20]=x130;

 sum[26]=x131;

 sum[32]=x132;

 sum[3]=x91;

 sum[9]=x117;

 sum[15]=x130;

 sum[21]=x104*x75 + x133*x8 + x136;

 sum[27]=x141;

 sum[33]=x142;

 sum[4]=x98;

 sum[10]=x123;

 sum[16]=x131;

 sum[22]=x141;

 sum[28]=x104*x94 + x136 + x137*x3;

 sum[34]=x143;

 sum[5]=x99;

 sum[11]=x124;

 sum[17]=x132;

 sum[23]=x142;

 sum[29]=x143;

 sum[35]=x136 + x138;
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  auto x0 = std::pow(p[0], 2)*std::pow(p[1], 8)/(std::pow(M_PI, 2)*std::pow(p[3], 4));
  auto x1 = std::pow(f[0], 3.0/2.0);
  auto x2 = (1.0/2.0)*f[3] - 1.0/2.0*f[4] + (1.0/2.0)*f[5];
  auto x3 = std::cos(x2);
  auto x4 = p[7]*x3;
  auto x5 = std::sin(x2);
  auto x6 = p[8]*x5;
  auto x7 = x4 - x6;
  auto x8 = std::sin(f[1]);
  auto x9 = std::sqrt(x8);
  auto x10 = std::sin(f[2]);
  auto x11 = std::sqrt(x10);
  auto x12 = x11*x9;
  auto x13 = x12*x7;
  auto x14 = std::sqrt(3);
  auto x15 = std::tan((1.0/2.0)*f[1]);
  auto x16 = std::sqrt(x15);
  auto x17 = 1.0/x16;
  auto x18 = std::tan((1.0/2.0)*f[2]);
  auto x19 = std::sqrt(x18);
  auto x20 = std::pow(p[3], 2)*(2*p[0]*std::pow(p[1], 4) + p[4]);
  auto x21 = p[2]/std::pow(M_PI, 3.0/2.0);
  auto x22 = x20*x21;
  auto x23 = x14*x17*x19*x22;
  auto x24 = (9.0/16.0)*x23;
  auto x25 = x13*x24;
  auto x26 = 1.0/f[0];
  auto x27 = p[6]*x8;
  auto x28 = std::pow(f[0], 2);
  auto x29 = x22*x28;
  auto x30 = (9.0/4.0)*x29;
  auto x31 = x27*x30;
  auto x32 = (3.0/2.0)*x1;
  auto x33 = x14*x22;
  auto x34 = x13*x33;
  auto x35 = x32*x34;
  auto x36 = std::pow(x15, -3.0/2.0);
  auto x37 = std::pow(x15, 2);
  auto x38 = -1.0/4.0*x37 - 1.0/4.0;
  auto x39 = x19*x38;
  auto x40 = x36*x39;
  auto x41 = x35*x40;
  auto x42 = (3.0/4.0)*x1;
  auto x43 = 1.0/x9;
  auto x44 = x11*x43;
  auto x45 = x23*x44;
  auto x46 = std::cos(f[1]);
  auto x47 = x46*x7;
  auto x48 = x45*x47;
  auto x49 = -x31 + x41 + x42*x48;
  auto x50 = f[0]*x20;
  auto x51 = x21*x50;
  auto x52 = (9.0/2.0)*x51;
  auto x53 = std::sqrt(f[0]);
  auto x54 = (9.0/4.0)*x53;
  auto x55 = x34*x40;
  auto x56 = (9.0/8.0)*x53;
  auto x57 = -x27*x52 + x48*x56 + x54*x55;
  auto x58 = -x26*x49 + x57;
  auto x59 = x26*x58;
  auto x60 = (9.0/2.0)*x22;
  auto x61 = p[6]*x60;
  auto x62 = 1.0/x53;
  auto x63 = (9.0/8.0)*x62;
  auto x64 = x24*x62;
  auto x65 = x44*x47;
  auto x66 = x55*x63 - x61*x8 + x64*x65;
  auto x67 = p[5]*x10;
  auto x68 = x30*x67;
  auto x69 = 1.0/x19;
  auto x70 = std::pow(x18, 2);
  auto x71 = (1.0/4.0)*x70;
  auto x72 = x71 + 1.0/4.0;
  auto x73 = x69*x72;
  auto x74 = x17*x73;
  auto x75 = x23*x7;
  auto x76 = 1.0/x11;
  auto x77 = x76*x9;
  auto x78 = std::cos(f[2]);
  auto x79 = x42*x78;
  auto x80 = x77*x79;
  auto x81 = x35*x74 - x68 + x75*x80;
  auto x82 = x34*x74;
  auto x83 = x77*x78;
  auto x84 = x56*x83;
  auto x85 = -x52*x67 + x54*x82 + x75*x84;
  auto x86 = -x26*x81 + x85;
  auto x87 = x26*x86;
  auto x88 = p[5]*x60;
  auto x89 = x7*x83;
  auto x90 = -x10*x88 + x63*x82 + x64*x89;
  auto x91 = 6*x26;
  auto x92 = std::pow(x8, 2);
  auto x93 = 1.0/x92;
  auto x94 = (1.0/6.0)*x92;
  auto x95 = std::pow(x46, 2);
  auto x96 = (1.0/9.0)*x95;
  auto x97 = x94 + x96;
  auto x98 = x93*x97;
  auto x99 = x91*x98;
  auto x100 = (2.0/3.0)*x93;
  auto x101 = x100*x95;
  auto x102 = x101*x26;
  auto x103 = -x102;
  auto x104 = x103 + x99;
  auto x105 = p[7]*x5;
  auto x106 = (1.0/2.0)*x105;
  auto x107 = p[8]*x3;
  auto x108 = (1.0/2.0)*x107;
  auto x109 = -x106 - x108;
  auto x110 = x12*x23;
  auto x111 = x109*x110;
  auto x112 = x111*x54;
  auto x113 = x111*x32;
  auto x114 = std::pow(x10, 2);
  auto x115 = 1.0/x114;
  auto x116 = 2*x115 + 2*x93 - 1;
  auto x117 = (1.0/3.0)*x116;
  auto x118 = x117*x26;
  auto x119 = std::pow(x78, 2);
  auto x120 = (2.0/3.0)*x115;
  auto x121 = x119*x120;
  auto x122 = x121*x26;
  auto x123 = x118*x46 - x122*x46 - x46*x99;
  auto x124 = -x104*x113 + x112 - x113*x123;
  auto x125 = x104*x124;
  auto x126 = -x122;
  auto x127 = x103 + x118 + x126;
  auto x128 = x112 - x113*x127;
  auto x129 = x123*x128;
  auto x130 = x111*x63;
  auto x131 = (1.0/6.0)*x114;
  auto x132 = (1.0/9.0)*x119;
  auto x133 = x131 + x132;
  auto x134 = x115*x133;
  auto x135 = x134*x91;
  auto x136 = x126 + x135;
  auto x137 = x106 + x108;
  auto x138 = x110*x137;
  auto x139 = x138*x54;
  auto x140 = x138*x32;
  auto x141 = -x102*x78 + x118*x78 - x135*x78;
  auto x142 = -x113*x141 - x136*x140 + x139;
  auto x143 = x136*x142;
  auto x144 = x128*x141;
  auto x145 = x138*x63;
  auto x146 = x127*x128;
  auto x147 = 1.0/x28;
  auto x148 = x147*x49 - x26*x57 - x59 + x66;
  auto x149 = p[6]*x46;
  auto x150 = x13*x23;
  auto x151 = -x150*x42;
  auto x152 = (3.0/8.0)*x1;
  auto x153 = x152*x75;
  auto x154 = std::pow(x8, 3.0/2.0);
  auto x155 = 1.0/x154;
  auto x156 = x11*x155;
  auto x157 = x156*x95;
  auto x158 = x37 + 1;
  auto x159 = x150*x152;
  auto x160 = std::pow(x15, -5.0/2.0);
  auto x161 = -3.0/4.0*x37 - 3.0/4.0;
  auto x162 = x160*x161*x39;
  auto x163 = x33*x40;
  auto x164 = x163*x65;
  auto x165 = -x149*x30 + x151 - x153*x157 - x158*x159 + x162*x35 + x164*x32;
  auto x166 = (2.0/3.0)*x20;
  auto x167 = x0/std::pow(f[0], 5);
  auto x168 = p[5]*x78;
  auto x169 = x168*x52;
  auto x170 = x149*x52;
  auto x171 = f[0]*x166 + x150*x54 + (27.0/4.0)*x167 + x169 + x170;
  auto x172 = (1.0/6.0)*f[0];
  auto x173 = x171*x172;
  auto x174 = x165 + x173;
  auto x175 = x174*x26;
  auto x176 = x150*x63 + x166 + x46*x61 + x78*x88 - 135.0/4.0*x0/std::pow(f[0], 6);
  auto x177 = x172*x176;
  auto x178 = x150*x53;
  auto x179 = x177 - 9.0/8.0*x178;
  auto x180 = x24*x7;
  auto x181 = x180*x53;
  auto x182 = x25*x53;
  auto x183 = x34*x54;
  auto x184 = -x157*x181 - x158*x182 + x162*x183 + x164*x54;
  auto x185 = -x165*x26 - x170 - x175 + x179 + x184;
  auto x186 = x36*x38*x73;
  auto x187 = x163*x7;
  auto x188 = x33*x74;
  auto x189 = x188*x65;
  auto x190 = x23*x47;
  auto x191 = x43*x76;
  auto x192 = x191*x78;
  auto x193 = x190*x192;
  auto x194 = x152*x193 + x186*x35 + x187*x80 + x189*x42;
  auto x195 = 2*x26;
  auto x196 = x183*x186 + x187*x84 + x189*x56 + x192*x24*x47*x53 - x194*x195;
  auto x197 = x12*x33;
  auto x198 = x197*x40;
  auto x199 = x198*x54;
  auto x200 = x109*x199;
  auto x201 = x1*x109;
  auto x202 = (3.0/2.0)*x201;
  auto x203 = x198*x202;
  auto x204 = x109*x45;
  auto x205 = x204*x46;
  auto x206 = x205*x42;
  auto x207 = x203 + x206;
  auto x208 = x205*x56;
  auto x209 = -x207*x26 + x208;
  auto x210 = 1.0/x8;
  auto x211 = (1.0/3.0)*x210;
  auto x212 = (1.0/6.0)*x116;
  auto x213 = (1.0/3.0)*x115;
  auto x214 = x119*x213;
  auto x215 = -x211*x95 - x212*x8 + x214*x8;
  auto x216 = (1.0/4.0)*x23;
  auto x217 = x201*x44;
  auto x218 = x217*x46;
  auto x219 = x216*x218;
  auto x220 = -x113*x215 + x203 - x219;
  auto x221 = (1.0/2.0)*x23;
  auto x222 = x217*x221;
  auto x223 = x203 + x219 + x222;
  auto x224 = -x104*x220 - x123*x223;
  auto x225 = (2.0/3.0)*x46;
  auto x226 = x210*x225;
  auto x227 = -x124*x226 - x128*x215;
  auto x228 = x200 + x209 + x224 + x227;
  auto x229 = x137*x32;
  auto x230 = x198*x229;
  auto x231 = x45*x46;
  auto x232 = x137*x231;
  auto x233 = x232*x42;
  auto x234 = x230 + x233;
  auto x235 = x222*x78;
  auto x236 = x234 - x235*x46 + x235;
  auto x237 = x137*x199;
  auto x238 = x232*x56;
  auto x239 = -x136*x236 - x141*x223 + x237 + x238;
  auto x240 = x124*x211;
  auto x241 = x211*x46;
  auto x242 = x128*x241;
  auto x243 = x240*x78 - x242*x78;
  auto x244 = -x234*x26 + x239 + x243;
  auto x245 = -x127*x223 + x200;
  auto x246 = x240 - x242;
  auto x247 = x209 + x245 + x246;
  auto x248 = x147*x81 - x26*x85 - x87 + x90;
  auto x249 = std::pow(x10, 3.0/2.0);
  auto x250 = 1.0/x249;
  auto x251 = x250*x9;
  auto x252 = x119*x251;
  auto x253 = x70 + 1;
  auto x254 = std::pow(x18, 3.0/2.0);
  auto x255 = 1.0/x254;
  auto x256 = -x71 - 1.0/4.0;
  auto x257 = x17*x255*x256*x72;
  auto x258 = x188*x89;
  auto x259 = x151 - x153*x252 + x159*x253 - x168*x30 + x257*x35 + x258*x32;
  auto x260 = x173 + x259;
  auto x261 = x26*x260;
  auto x262 = -x181*x252 + x182*x253 + x183*x257 + x258*x54;
  auto x263 = -x169 + x179 - x259*x26 - x261 + x262;
  auto x264 = x197*x74;
  auto x265 = x202*x264;
  auto x266 = x109*x23;
  auto x267 = x266*x80;
  auto x268 = x265 + x267;
  auto x269 = x266*x84;
  auto x270 = -x26*x268 + x269;
  auto x271 = x264*x54;
  auto x272 = x109*x271;
  auto x273 = x1*x221;
  auto x274 = x137*x77;
  auto x275 = x273*x274;
  auto x276 = x201*x221;
  auto x277 = x46*x83;
  auto x278 = x268 + x275*x46 - x276*x277;
  auto x279 = x201*x216;
  auto x280 = x265 + x275 + x279*x83;
  auto x281 = -x104*x278 - x123*x280 + x272;
  auto x282 = 1.0/x10;
  auto x283 = (1.0/3.0)*x282;
  auto x284 = x142*x283;
  auto x285 = x283*x78;
  auto x286 = x128*x285;
  auto x287 = x284*x46 - x286*x46;
  auto x288 = x270 + x281 + x287;
  auto x289 = x137*x271;
  auto x290 = x229*x264;
  auto x291 = x23*x274;
  auto x292 = x291*x79;
  auto x293 = x290 + x292;
  auto x294 = (1.0/3.0)*x93;
  auto x295 = x294*x95;
  auto x296 = -x10*x212 + x10*x295 - x119*x283;
  auto x297 = (1.0/4.0)*x1;
  auto x298 = x291*x78;
  auto x299 = -x113*x296 + x290 - x297*x298;
  auto x300 = -x136*x299 - x141*x280;
  auto x301 = (2.0/3.0)*x78;
  auto x302 = x282*x301;
  auto x303 = -x128*x296 - x142*x302 + x298*x56;
  auto x304 = -x26*x293 + x289 + x300 + x303;
  auto x305 = -x127*x280 + x272;
  auto x306 = x284 - x286;
  auto x307 = x270 + x305 + x306;
  auto x308 = 6*x147;
  auto x309 = x308*x98;
  auto x310 = x101*x147;
  auto x311 = x117*x147;
  auto x312 = x121*x147;
  auto x313 = -x104*x112 - x112*x123 - x113*(-x309 + x310) - x113*(x309*x46 - x311*x46 + x312*x46) - x125 - x129 + x130;
  auto x314 = x210*x46;
  auto x315 = (4.0/3.0)*x26;
  auto x316 = std::pow(x46, 3);
  auto x317 = std::pow(x8, -3);
  auto x318 = x316*x317;
  auto x319 = x315*x318;
  auto x320 = 12*x26;
  auto x321 = x317*x320*x97;
  auto x322 = x315*x317;
  auto x323 = (2.0/3.0)*x26;
  auto x324 = x200 - x220*x26;
  auto x325 = -x104*x203 - x104*x206 - x113*(x195*x314 + x319 - x321*x46) - x113*(-x118*x8 + x122*x8 - x210*x323*x95 + x210*x91*x97 + x321*x95 - x322*x95) - x123*x203 - x123*x206 + x208 + x227 + x324;
  auto x326 = std::pow(x10, -3);
  auto x327 = x315*x326;
  auto x328 = x327*x78;
  auto x329 = std::pow(x78, 3);
  auto x330 = x327*x329;
  auto x331 = x282*x78;
  auto x332 = x315*x331;
  auto x333 = -x26*x278 + x269;
  auto x334 = -x104*x265 - x104*x267 - x113*(-x328*x46 + x330*x46 + x332*x46) - x123*x265 - x123*x267 + x272 + x287 + x333;
  auto x335 = f[0]*x97;
  auto x336 = x176*x335;
  auto x337 = (1.0/3.0)*x8;
  auto x338 = x337*x49;
  auto x339 = (1.0/4.0)*x4;
  auto x340 = (1.0/4.0)*x6;
  auto x341 = -x339 + x340;
  auto x342 = x110*x341;
  auto x343 = x32*x342;
  auto x344 = x171*x335 + x338*x46 + x343;
  auto x345 = x104*x344;
  auto x346 = (1.0/9.0)*f[0];
  auto x347 = x171*x346;
  auto x348 = x347*x46;
  auto x349 = -x338 + x343 + x348;
  auto x350 = x123*x349;
  auto x351 = x337*x58;
  auto x352 = x342*x54;
  auto x353 = -x104*x343 - x123*x343 + x352;
  auto x354 = x336 - x345 - x350 + x351*x46 + x353;
  auto x355 = x339 - x340;
  auto x356 = x110*x355;
  auto x357 = x32*x356;
  auto x358 = x356*x54;
  auto x359 = (1.0/3.0)*x10;
  auto x360 = x359*x86;
  auto x361 = x176*x346;
  auto x362 = x361*x46;
  auto x363 = x362*x78;
  auto x364 = -x351*x78 + x358 - x360*x46 + x363;
  auto x365 = x359*x81;
  auto x366 = -x338*x78 + x348*x78 + x357 - x365*x46;
  auto x367 = -x136*x366 - x141*x349;
  auto x368 = -x104*x357 - x123*x357 + x364 + x367;
  auto x369 = -x127*x349;
  auto x370 = -x351 + x362;
  auto x371 = x353 + x369 + x370;
  auto x372 = x134*x308;
  auto x373 = -x112*x141 - x113*(x310*x78 - x311*x78 + x372*x78) - x136*x139 - x140*(x312 - x372) - x143 - x144 + x145;
  auto x374 = -x236*x26;
  auto x375 = x322*x46;
  auto x376 = x314*x315;
  auto x377 = -x113*(x319*x78 - x375*x78 + x376*x78) - x136*x230 - x136*x233 - x141*x203 - x141*x206 + x237 + x238 + x243 + x374;
  auto x378 = x133*x320*x326;
  auto x379 = -x26*x299 + x289;
  auto x380 = -x113*(x10*x102 - x10*x118 - x119*x282*x323 - x119*x327 + x119*x378 + x133*x282*x91) - x136*x290 - x136*x292 - x140*(x195*x331 + x330 - x378*x78) - x141*x265 - x141*x267 + x303 + x379;
  auto x381 = -x136*x357 - x141*x343;
  auto x382 = x347*x78 + x357 - x365;
  auto x383 = -x104*x366 - x123*x382;
  auto x384 = x364 + x381 + x383;
  auto x385 = f[0]*x133;
  auto x386 = x171*x385 + x343 + x365*x78;
  auto x387 = x136*x386;
  auto x388 = x141*x382;
  auto x389 = x176*x385 + x352;
  auto x390 = -x136*x343 - x141*x357 + x360*x78 - x387 - x388 + x389;
  auto x391 = -x127*x382;
  auto x392 = x361*x78;
  auto x393 = x358 - x360 + x392;
  auto x394 = x381 + x391 + x393;
  auto x395 = -x112*x127 - x113*(-x147*(x100 + x120 - 1.0/3.0) + x310 + x312) + x130 - x146;
  auto x396 = -x223*x26;
  auto x397 = -x113*(x319 - x375 + x376) - x127*x203 - x127*x206 + x200 + x208 + x246 + x396;
  auto x398 = -x26*x280;
  auto x399 = -x113*(-x328 + x330 + x332) - x127*x265 - x127*x267 + x269 + x272 + x306 + x398;
  auto x400 = -x127*x343 + x352;
  auto x401 = x343 + x347;
  auto x402 = -x104*x349 - x123*x401;
  auto x403 = x370 + x400 + x402;
  auto x404 = -x136*x382 - x141*x401;
  auto x405 = -x127*x357 + x393 + x404;
  auto x406 = x127*x401;
  auto x407 = x361 + x400 - x406;
  auto x408 = (15.0/4.0)*x51;
  auto x409 = (3.0/4.0)*x51;
  auto x410 = (9.0/8.0)*x167 + x177 - 3.0/4.0*x178 + x20*x346;
  auto x411 = (1.0/3.0)*f[0];
  auto x412 = x172*x57;
  auto x413 = (9.0/4.0)*x1;
  auto x414 = x152*x34;
  auto x415 = x158*x19*x414;
  auto x416 = x152*x48;
  auto x417 = x1*x180;
  auto x418 = (9.0/8.0)*x1;
  auto x419 = x1*x158;
  auto x420 = x161*x35;
  auto x421 = x33*x65;
  auto x422 = x172*x85;
  auto x423 = x153*x83;
  auto x424 = x188*x7;
  auto x425 = (3.0/16.0)*x1;
  auto x426 = x75*x83;
  auto x427 = x186*x32;
  auto x428 = x191*x47*x79;
  auto x429 = -x152*x157*x424 - x152*x158*x82 - x155*x425*x75*x76*x78*x95 + x160*x38*x420*x73 + x162*x33*x7*x80 + x163*x428 - 3.0/16.0*x419*x426 - x42*x82 + x421*x427 - x423;
  auto x430 = x215*x223;
  auto x431 = (4.0/3.0)*x314;
  auto x432 = x111*x152;
  auto x433 = -x432;
  auto x434 = x157*x23;
  auto x435 = x201*x434;
  auto x436 = x163*x218;
  auto x437 = x162*x197;
  auto x438 = -x158*x432 + x202*x437;
  auto x439 = x433 - 3.0/8.0*x435 + (3.0/2.0)*x436 + x438;
  auto x440 = x220*x301;
  auto x441 = x223*x226;
  auto x442 = x138*x152;
  auto x443 = -x442;
  auto x444 = x137*x152;
  auto x445 = x44*x46;
  auto x446 = -x158*x442 + x163*x229*x445 + x229*x437 - x434*x444;
  auto x447 = (2.0/3.0)*x220;
  auto x448 = x172*x86 + x429;
  auto x449 = x256*x35*x72;
  auto x450 = x33*x89;
  auto x451 = -x119*x190*x250*x425*x43 - x152*x187*x252 + x152*x253*x55 + x188*x428 + x253*x425*x48 + x255*x36*x38*x449 + x257*x42*x421 - x416 - x42*x55 + x427*x450;
  auto x452 = x172*x58 + x451;
  auto x453 = x236*x283;
  auto x454 = x223*x285;
  auto x455 = x453*x46 - x454*x46;
  auto x456 = x188*x42*x445;
  auto x457 = x109*x456;
  auto x458 = x201*x43;
  auto x459 = x458*x76;
  auto x460 = x23*x46;
  auto x461 = x459*x460*x78;
  auto x462 = x186*x197;
  auto x463 = x202*x462;
  auto x464 = x109*x163*x80 + x463;
  auto x465 = x457 + (3.0/8.0)*x461 + x464;
  auto x466 = -x215*x280 - x226*x278 + x465;
  auto x467 = x455 + x466;
  auto x468 = x211*x278;
  auto x469 = x241*x280;
  auto x470 = x229*x462;
  auto x471 = x137*x456;
  auto x472 = x468*x78 - x469*x78 + x470 + x471;
  auto x473 = x163*x274;
  auto x474 = x192*x460;
  auto x475 = -x223*x296 - x236*x302 + x444*x474 + x473*x79;
  auto x476 = x472 + x475;
  auto x477 = x453 - x454;
  auto x478 = x468 - x469;
  auto x479 = x465 + x477 + x478;
  auto x480 = (3.0/8.0)*x53;
  auto x481 = x205*x480;
  auto x482 = -x112*x215 + x224 + x324 - x481;
  auto x483 = x111*x297;
  auto x484 = x124*x172;
  auto x485 = -x213 - x294 + 1.0/6.0;
  auto x486 = (1.0/8.0)*x435;
  auto x487 = (1.0/2.0)*x436;
  auto x488 = -x113*(x100*x46 + x214*x46 + x225 + x294*x316 + x46*x485) - x203*x215 - x206*x215 - x220*x226 - x430 + x438 + x483 + x484 + x486 + x487;
  auto x489 = (2.0/3.0)*x8;
  auto x490 = x301*x8;
  auto x491 = x188*x218;
  auto x492 = (1.0/4.0)*x491;
  auto x493 = (1.0/8.0)*x461;
  auto x494 = -x113*(-x302*x8 - x326*x329*x489 + x326*x490) - x215*x265 - x215*x267 + x455 + x464 - x492 - x493;
  auto x495 = x215*x349;
  auto x496 = x174*x337;
  auto x497 = x198*x32;
  auto x498 = x341*x497;
  auto x499 = x231*x297;
  auto x500 = x341*x499;
  auto x501 = -x215*x343 + x498 - x500;
  auto x502 = -x226*x344 + x335*x58 + x46*x496 - x495 + x501;
  auto x503 = x355*x499;
  auto x504 = x355*x497;
  auto x505 = x211*x344;
  auto x506 = x211*x349;
  auto x507 = x46*x506;
  auto x508 = x504 + x505*x78 - x507*x78;
  auto x509 = x194*x359;
  auto x510 = x346*x58;
  auto x511 = x46*x510;
  auto x512 = -x46*x509 - x496*x78 + x511*x78;
  auto x513 = -x215*x357 - x503 + x508 + x512;
  auto x514 = -x496 + x511;
  auto x515 = x501 + x505 - x507 + x514;
  auto x516 = (3.0/4.0)*x53;
  auto x517 = x204*x516;
  auto x518 = x517*x78;
  auto x519 = x239 + x374 - x46*x518 + x518;
  auto x520 = x142*x172;
  auto x521 = x211*x220;
  auto x522 = x223*x241;
  auto x523 = (1.0/2.0)*x1;
  auto x524 = x111*x523;
  auto x525 = (1.0/2.0)*x217;
  auto x526 = x163*x525;
  auto x527 = x279*x46;
  auto x528 = x156*x527;
  auto x529 = -x138*x42 + x157*x279*x78 + x446 - x487*x78 + x520 + x521*x78 - x522*x78 + x524*x78 + x526*x78 - x528*x78;
  auto x530 = x221*x249*x458;
  auto x531 = x216*x459;
  auto x532 = x119*x531;
  auto x533 = x188*x525;
  auto x534 = x46*x530 - x46*x532 + x470 + x471 + x475 - 1.0/2.0*x491*x78 - x530 + x532 + x533*x78;
  auto x535 = x273*x44;
  auto x536 = x341*x535;
  auto x537 = x536*x78;
  auto x538 = -x46*x537 + x537;
  auto x539 = x231*x42;
  auto x540 = x355*x539;
  auto x541 = -x215*x382 - x226*x366 + x540;
  auto x542 = x504 + x512 + x538 + x541;
  auto x543 = x211*x366;
  auto x544 = x241*x382;
  auto x545 = x355*x535;
  auto x546 = x545*x78;
  auto x547 = x341*x539 + x498;
  auto x548 = x509*x78 + x547;
  auto x549 = x385*x58 - x46*x546 + x543*x78 - x544*x78 + x546 + x548;
  auto x550 = -x509;
  auto x551 = x504 + x510*x78 + x550;
  auto x552 = x540 + x543 - x544;
  auto x553 = x538 + x551 + x552;
  auto x554 = x245 + x396 + x481 + x517;
  auto x555 = x128*x172 - x483;
  auto x556 = x436 + x438 - x486 + x521 - x522 + x526 - x528 + x555;
  auto x557 = x531*x78;
  auto x558 = x464 + x477 + x492 + x493 + x533 + x557;
  auto x559 = -x215*x401;
  auto x560 = x210*x349;
  auto x561 = x498 + x500 + x536;
  auto x562 = -x225*x560 + x514 + x559 + x561;
  auto x563 = x241*x401;
  auto x564 = x506*x78 - x563*x78;
  auto x565 = x503 + x545 + x551 + x564;
  auto x566 = x506 + x510 + x561 - x563;
  auto x567 = x17*x253*x414;
  auto x568 = x225*x282;
  auto x569 = x225*x331;
  auto x570 = x23*x252;
  auto x571 = x201*x570;
  auto x572 = x188*x83;
  auto x573 = x197*x257;
  auto x574 = x202*x573 + x253*x432;
  auto x575 = x202*x572 - 3.0/8.0*x571 + x574;
  auto x576 = x433 + x575;
  auto x577 = x280*x296;
  auto x578 = (4.0/3.0)*x331;
  auto x579 = x229*x573 + x253*x442;
  auto x580 = (2.0/3.0)*x299;
  auto x581 = x291*x516;
  auto x582 = x266*x83;
  auto x583 = x281 + x333 - x46*x516*x582 + x46*x581;
  auto x584 = x154*x76;
  auto x585 = x23*x297;
  auto x586 = x137*x585;
  auto x587 = x191*x586;
  auto x588 = x274*x523;
  auto x589 = x163*x588;
  auto x590 = x163*x201;
  auto x591 = -x137*x273*x584 + x276*x584*x78 - 1.0/2.0*x277*x590 + x46*x589 + x466 - x557*x95 + x587*x95;
  auto x592 = x283*x299;
  auto x593 = x280*x285;
  auto x594 = x188*x588;
  auto x595 = x251*x586*x78;
  auto x596 = x201*x572;
  auto x597 = -x111*x42 + x252*x527 + x46*x524 + x46*x592 - x46*x593 + x46*x594 - x46*x595 - 1.0/2.0*x46*x596 + x484 + x575;
  auto x598 = x283*x366;
  auto x599 = x285*x349;
  auto x600 = x194*x337;
  auto x601 = x264*x32;
  auto x602 = x341*x601;
  auto x603 = x23*x80;
  auto x604 = x341*x603 + x602;
  auto x605 = x46*x600 + x604;
  auto x606 = x273*x77;
  auto x607 = x355*x606;
  auto x608 = x46*x607;
  auto x609 = x341*x606;
  auto x610 = x46*x609;
  auto x611 = x608 - x610*x78;
  auto x612 = x335*x86 + x46*x598 - x46*x599 + x605 + x611;
  auto x613 = x355*x601;
  auto x614 = x260*x359;
  auto x615 = x346*x86;
  auto x616 = x46*x615;
  auto x617 = -x46*x614 - x600*x78 + x616*x78;
  auto x618 = x355*x603;
  auto x619 = -x296*x349 - x302*x366 + x618;
  auto x620 = -x608*x78 + x610 + x613 + x617 + x619;
  auto x621 = -x600;
  auto x622 = x616 + x621;
  auto x623 = x598 - x599 + x604;
  auto x624 = x611 + x622 + x623;
  auto x625 = -x112*x296 - x298*x480 + x300 + x379;
  auto x626 = (2.0/3.0)*x10;
  auto x627 = x10*x225;
  auto x628 = (1.0/8.0)*x1;
  auto x629 = x137*x628;
  auto x630 = -x113*(-x10*x226 + x317*x627 - x318*x626) - x203*x296 - x206*x296 - x297*x473*x78 + x472 - x474*x629;
  auto x631 = x138*x297;
  auto x632 = -x113*(x120*x78 + x213*x329 + x295*x78 + x301 + x485*x78) - x265*x296 - x267*x296 - x299*x302 + x520 + x570*x629 - x577 + x579 + x594*x78 + x631;
  auto x633 = x585*x83;
  auto x634 = x355*x633;
  auto x635 = -x296*x343 - x634;
  auto x636 = x283*x386;
  auto x637 = x283*x382;
  auto x638 = x637*x78;
  auto x639 = x46*x636 - x46*x638 + x613;
  auto x640 = x617 + x635 + x639;
  auto x641 = x296*x382;
  auto x642 = x341*x633;
  auto x643 = -x296*x357 - x302*x386 + x385*x86 + x602 + x614*x78 - x641 - x642;
  auto x644 = x613 + x636;
  auto x645 = -x614 + x615*x78;
  auto x646 = x635 - x638 + x644 + x645;
  auto x647 = x305 + x398 + x480*x582 + x581;
  auto x648 = x457 + x46*x587 + x463 + x478 + x493 + x589 + (1.0/4.0)*x590*x83;
  auto x649 = x555 - 1.0/8.0*x571 + x574 + x592 - x593 + x594 - x595 + x596;
  auto x650 = x602 + x607 + x642;
  auto x651 = x285*x401;
  auto x652 = x46*x637 - x46*x651;
  auto x653 = x622 + x650 + x652;
  auto x654 = -x296*x401;
  auto x655 = x282*x382;
  auto x656 = -x301*x655 + x609 + x613 + x634 + x645 + x654;
  auto x657 = x615 + x637 + x650 - x651;
  auto x658 = x337*x57;
  auto x659 = (3.0/4.0)*x29;
  auto x660 = x27*x659;
  auto x661 = x523*x55;
  auto x662 = x297*x48;
  auto x663 = x660 - x661 - x662;
  auto x664 = x149*x659;
  auto x665 = x150*x297;
  auto x666 = x628*x75;
  auto x667 = x157*x666;
  auto x668 = x150*x628;
  auto x669 = x158*x668;
  auto x670 = x34*x523;
  auto x671 = x162*x670;
  auto x672 = x164*x523;
  auto x673 = x46*x8;
  auto x674 = (1.0/8.0)*x105;
  auto x675 = (1.0/8.0)*x107;
  auto x676 = x110*x32;
  auto x677 = x676*(x674 + x675);
  auto x678 = x207*x337;
  auto x679 = x413*x97;
  auto x680 = x111*x679 + x46*x678 + x677;
  auto x681 = x234*x337;
  auto x682 = (2.0/9.0)*f[0];
  auto x683 = x46*x682;
  auto x684 = x124*x683;
  auto x685 = x676*(-x674 - x675);
  auto x686 = (2.0/27.0)*x50;
  auto x687 = (3.0/4.0)*x167;
  auto x688 = (1.0/2.0)*x51;
  auto x689 = x168*x688;
  auto x690 = x149*x688;
  auto x691 = (1.0/4.0)*x178;
  auto x692 = -x686 - x687 - x689 - x690 - x691;
  auto x693 = x46*x692;
  auto x694 = x359*x85;
  auto x695 = x358 + x363 + x367 + x383 - x46*x694 - x658*x78 - x693*x78;
  auto x696 = x659*x67;
  auto x697 = x523*x82;
  auto x698 = x297*x426;
  auto x699 = -x696 + x697 + x698;
  auto x700 = x10*x8;
  auto x701 = -x660 + x661 + x662;
  auto x702 = x46*x701;
  auto x703 = x163*x297*x89 + x186*x670 + x189*x297 + x193*x628;
  auto x704 = x10*x46;
  auto x705 = -x665;
  auto x706 = x8*(-x664 - x667 - x669 + x671 + x672 + x705);
  auto x707 = x78*x8;
  auto x708 = x346*x57;
  auto x709 = x46*x708;
  auto x710 = -x347*x707 + x508 + x541 + x699*x700 - x702*x78 - x703*x704 - x706*x78 + x709*x78;
  auto x711 = x699*x78;
  auto x712 = x168*x659;
  auto x713 = x252*x666;
  auto x714 = x253*x668;
  auto x715 = x257*x670;
  auto x716 = x258*x523;
  auto x717 = x10*(x705 - x712 - x713 + x714 + x715 + x716);
  auto x718 = x346*x85;
  auto x719 = x46*x718;
  auto x720 = -x347*x704 - x46*x711 - x46*x717 + x619 + x639 + x700*x701 - x703*x707 + x719*x78;
  auto x721 = x220*x337;
  auto x722 = x278*x359;
  auto x723 = x236*x337;
  auto x724 = x124*x346;
  auto x725 = x724*x78;
  auto x726 = x268*x359;
  auto x727 = x46*x483;
  auto x728 = -x46*x726 - x678*x78 + x727*x78;
  auto x729 = x142*x335 - x46*x722 + x46*x723 + x46*x725 + x685 - x721*x78 + x728;
  auto x730 = x293*x359;
  auto x731 = x299*x359;
  auto x732 = x142*x346;
  auto x733 = x46*x732;
  auto x734 = x46*x631;
  auto x735 = x124*x385 - x46*x730 - x46*x731 + x677 - x681*x78 + x722*x78 - x723*x78 + x733*x78 + x734*x78;
  auto x736 = -x723 + x733;
  auto x737 = x685 - x722 + x725;
  auto x738 = x728 + x736 + x737;
  auto x739 = x352 + x362 + x369 + x402 - x658 - x693;
  auto x740 = -x347*x8 - x46*x560 + x505 + x547 + x559 - x702 - x706 + x709;
  auto x741 = x621 + x623 + x652 + x719;
  auto x742 = x223*x337;
  auto x743 = x677 - x678 + x727;
  auto x744 = x128*x335 + x46*x724 + x46*x742 - x721 + x743;
  auto x745 = x280*x359;
  auto x746 = x128*x346;
  auto x747 = x46*x746;
  auto x748 = -x46*x745 - x742*x78 + x747*x78;
  auto x749 = -x681 + x734 + x737 + x748;
  auto x750 = x724 - x742 + x743 + x747;
  auto x751 = x696 - x697 - x698;
  auto x752 = x10*x78;
  auto x753 = x142*x78;
  auto x754 = x133*x413;
  auto x755 = x111*x754 + x677 + x726*x78;
  auto x756 = x358 + x391 + x392 + x404 - x692*x78 - x694;
  auto x757 = x504 + x550 + x552 + x564 + x708*x78;
  auto x758 = -x10*x347 + x618 + x644 + x654 - x655*x78 - x711 - x717 + x718*x78;
  auto x759 = x483*x78 + x685 - x726;
  auto x760 = x736 + x748 + x759;
  auto x761 = x128*x385 + x631*x78 + x677 - x730 - x731 + x732*x78 + x745*x78;
  auto x762 = x732 - x745 + x746*x78 + x759;
  auto x763 = x128*x682;
  auto x764 = x483 + x677;

 sum[0]=-x25/x1 + (405.0/2.0)*x0/std::pow(f[0], 7);

 sum[36]=x148;

 sum[72]=x248;

 sum[108]=x313;

 sum[144]=x373;

 sum[180]=x395;

 sum[6]=x148;

 sum[42]=-x149*x408 + x168*x409 - 2*x175 + x184 + x410;

 sum[78]=x196;

 sum[114]=x482;

 sum[150]=x519;

 sum[186]=x554;

 sum[12]=x248;

 sum[48]=x196;

 sum[84]=x149*x409 - x168*x408 - 2*x261 + x262 + x410;

 sum[120]=x583;

 sum[156]=x625;

 sum[192]=x647;

 sum[18]=x313;

 sum[54]=x482;

 sum[90]=x583;

 sum[126]=-x171*(-x94 - x96) + x336 - 2*x345 - 2*x350 + x352 + x46*x658;

 sum[162]=x695;

 sum[198]=x739;

 sum[24]=x373;

 sum[60]=x519;

 sum[96]=x625;

 sum[132]=x695;

 sum[168]=-x171*(-x131 - x132) - 2*x387 - 2*x388 + x389 + x694*x78;

 sum[204]=x756;

 sum[30]=x395;

 sum[66]=x554;

 sum[102]=x647;

 sum[138]=x739;

 sum[174]=x756;

 sum[210]=x352 + x361 - 2*x406 + x686 + x687 + x689 + x690 + x691;

 sum[1]=-2*x59 + x66;

 sum[37]=x185;

 sum[73]=x196;

 sum[109]=x325;

 sum[145]=x377;

 sum[181]=x397;

 sum[7]=x185;

 sum[43]=x11*x316*x417/std::pow(x8, 5.0/2.0) - x157*x187*x418 - x158*x41 - x16*x415 - x161*x36*x415 + x162*x413*x421 - x24*x419*x65 + x31 + x411*x58 + x412 - x413*x55 + x416 + x39*x420*(-5.0/4.0*x37 - 5.0/4.0)/std::pow(x15, 7.0/2.0);

 sum[79]=x448;

 sum[115]=x488;

 sum[151]=x529;

 sum[187]=x556;

 sum[13]=x196;

 sum[49]=x448;

 sum[85]=x412 + x451;

 sum[121]=x591;

 sum[157]=x630;

 sum[193]=x648;

 sum[19]=x325;

 sum[55]=x488;

 sum[91]=x591;

 sum[127]=x335*x57 - x344*x431 + x347*x673 - 2*x495 + x547 + x663*x92 - x663*x95 - x673*(x664 + x665 + x667 + x669 - x671 - x672);

 sum[163]=x710;

 sum[199]=x740;

 sum[25]=x377;

 sum[61]=x529;

 sum[97]=x630;

 sum[133]=x710;

 sum[169]=x210*x301*x366 - x226*x382*x78 + x385*x57 + x548;

 sum[205]=x757;

 sum[31]=x397;

 sum[67]=x556;

 sum[103]=x648;

 sum[139]=x740;

 sum[175]=x757;

 sum[211]=-x226*x401 + x547 + (2.0/3.0)*x560 + x708;

 sum[2]=-2*x87 + x90;

 sum[38]=x196;

 sum[74]=x263;

 sum[110]=x334;

 sum[146]=x380;

 sum[182]=x399;

 sum[8]=x196;

 sum[44]=x422 + x429;

 sum[80]=x452;

 sum[116]=x494;

 sum[152]=x534;

 sum[188]=x558;

 sum[14]=x263;

 sum[50]=x452;

 sum[86]=x1*x24*x253*x89 + x17*x449*(-3.0/4.0*x70 - 3.0/4.0)/std::pow(x18, 5.0/2.0) - x252*x418*x424 + x254*x567 + x256*x567*x69 + x257*x413*x450 + x411*x86 - x413*x82 + x422 + x423 + x68 + x329*x417*x9/std::pow(x10, 5.0/2.0);

 sum[122]=x597;

 sum[158]=x632;

 sum[194]=x649;

 sum[20]=x334;

 sum[56]=x494;

 sum[92]=x597;

 sum[128]=x335*x85 - x349*x569 + x366*x568 + x605;

 sum[164]=x720;

 sum[200]=x741;

 sum[26]=x380;

 sum[62]=x534;

 sum[98]=x632;

 sum[134]=x720;

 sum[170]=x114*x751 - x119*x751 + x347*x752 + x385*x85 - x386*x578 + x604 - 2*x641 - x752*(x665 + x712 + x713 - x714 - x715 - x716);

 sum[206]=x758;

 sum[32]=x399;

 sum[68]=x558;

 sum[104]=x649;

 sum[140]=x741;

 sum[176]=x758;

 sum[212]=-x302*x401 + x604 + (2.0/3.0)*x655 + x718;

 sum[3]=-2*x125 - 2*x129 + x130;

 sum[39]=x228;

 sum[75]=x288;

 sum[111]=x354;

 sum[147]=x384;

 sum[183]=x403;

 sum[9]=x228;

 sum[45]=-x220*x431 - 2*x430 + x439;

 sum[81]=x467;

 sum[117]=x502;

 sum[153]=x542;

 sum[189]=x562;

 sum[15]=x288;

 sum[51]=x467;

 sum[87]=-x280*x569 + x299*x568 + x576;

 sum[123]=x612;

 sum[159]=x640;

 sum[195]=x653;

 sum[21]=x354;

 sum[57]=x502;

 sum[93]=x612;

 sum[129]=2*x124*x335 + x220*x225*x8 + x680;

 sum[165]=x729;

 sum[201]=x744;

 sum[27]=x384;

 sum[63]=x542;

 sum[99]=x640;

 sum[135]=x729;

 sum[171]=-x236*x490 - x299*x627 + x683*x753 + x755;

 sum[207]=x760;

 sum[33]=x403;

 sum[69]=x562;

 sum[105]=x653;

 sum[141]=x744;

 sum[177]=x760;

 sum[213]=-x223*x489 + x46*x763 + x764;

 sum[4]=-2*x143 - 2*x144 + x145;

 sum[40]=x244;

 sum[76]=x304;

 sum[112]=x368;

 sum[148]=x390;

 sum[184]=x405;

 sum[10]=x244;

 sum[46]=x210*x440 - x441*x78 + x443 + x446;

 sum[82]=x476;

 sum[118]=x513;

 sum[154]=x549;

 sum[190]=x565;

 sum[16]=x304;

 sum[52]=x476;

 sum[88]=x188*x274*x32*x78 - x299*x578 + x443 - x444*x570 - 2*x577 + x579;

 sum[124]=x620;

 sum[160]=x643;

 sum[196]=x656;

 sum[22]=x368;

 sum[58]=x513;

 sum[94]=x620;

 sum[130]=x138*x679 - x278*x627 - x440*x8 + x46*x681 + x684*x78 + x685;

 sum[166]=x735;

 sum[202]=x749;

 sum[28]=x390;

 sum[64]=x549;

 sum[100]=x643;

 sum[136]=x735;

 sum[172]=x10*x299*x301 + x138*x754 + 2*x142*x385 + x685 + x730*x78;

 sum[208]=x761;

 sum[34]=x405;

 sum[70]=x565;

 sum[106]=x656;

 sum[142]=x749;

 sum[178]=x761;

 sum[214]=-x280*x626 + x631 + x685 + x763*x78;

 sum[5]=x130 - 2*x146;

 sum[41]=x247;

 sum[77]=x307;

 sum[113]=x371;

 sum[149]=x394;

 sum[185]=x407;

 sum[11]=x247;

 sum[47]=x210*x447 + x439 - x441;

 sum[83]=x479;

 sum[119]=x515;

 sum[155]=x553;

 sum[191]=x566;

 sum[17]=x307;

 sum[53]=x479;

 sum[89]=-x280*x302 + x282*x580 + x576;

 sum[125]=x624;

 sum[161]=x646;

 sum[197]=x657;

 sum[23]=x371;

 sum[59]=x515;

 sum[95]=x624;

 sum[131]=-x447*x8 + x680 + x684;

 sum[167]=x738;

 sum[203]=x750;

 sum[29]=x394;

 sum[65]=x553;

 sum[101]=x646;

 sum[137]=x738;

 sum[173]=-x10*x580 + x682*x753 + x755;

 sum[209]=x762;

 sum[35]=x407;

 sum[71]=x566;

 sum[107]=x657;

 sum[143]=x750;

 sum[179]=x762;

 sum[215]=x763 + x764;
       
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