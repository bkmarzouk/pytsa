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
// Potential file rewriten at Sat May 21 20:33:28 2022

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
nF=2;
nP=3;

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
  sum=(1.0/2.0)*std::pow(f[0], 2)*std::pow(p[0], 2) + (1.0/2.0)*std::pow(f[1], 2)*std::pow(p[1], 2);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot

 sum[0]=f[0]*std::pow(p[0], 2);

 sum[1]=f[1]*std::pow(p[1], 2);
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot
  auto x0 = std::pow(p[0], 2);
  auto x1 = std::pow(p[1], 2);

 sum[0]=x0;

 sum[2]=-f[1]*x1/std::tan(f[0]);

 sum[1]=0;

 sum[3]=(1.0/2.0)*f[0]*x0*std::sin(2*f[0]) + x1;
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  auto x0 = std::tan(f[0]);
  auto x1 = 2/std::pow(x0, 2);
  auto x2 = std::pow(p[1], 2);
  auto x3 = f[1]*x2;
  auto x4 = x3*(x1 + 1);
  auto x5 = 1.0/x0;
  auto x6 = std::pow(p[0], 2);
  auto x7 = 2*f[0];
  auto x8 = std::sin(x7);
  auto x9 = x6*x8;
  auto x10 = -f[0]*x9;
  auto x11 = (1.0/2.0)*x5*(x0*x9 + x10 - 4*x2);

 sum[0]=0;

 sum[4]=x1*x3;

 sum[2]=x4;

 sum[6]=x11;

 sum[1]=x4;

 sum[5]=x11;

 sum[3]=x5*((1.0/2.0)*x0*x6*(x7*std::cos(x7) + x8) + x10 - 2*x2);

 sum[7]=0;
       
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