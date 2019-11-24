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
// Potential file rewriten at Sun Nov 24 13:03:49 2019

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
nF=3;
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
  sum=(1.0/2.0)*std::pow(f[0], 2)*std::pow(p[0], 2) + (1.0/2.0)*std::pow(f[1], 2)*std::pow(p[1], 2) + (1.0/2.0)*std::pow(f[2], 2)*std::pow(p[2], 2);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot

 sum[0]=f[0]*std::pow(p[0], 2);

 sum[1]=f[1]*std::pow(p[1], 2);

 sum[2]=f[2]*std::pow(p[2], 2);
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot

 sum[0]=std::pow(p[0], 2);

 sum[3]=0;

 sum[6]=0;

 sum[1]=0;

 sum[4]=std::pow(p[1], 2);

 sum[7]=0;

 sum[2]=0;

 sum[5]=0;

 sum[8]=std::pow(p[2], 2);
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot

 sum[0]=0;

 sum[9]=0;

 sum[18]=0;

 sum[3]=0;

 sum[12]=0;

 sum[21]=0;

 sum[6]=0;

 sum[15]=0;

 sum[24]=0;

 sum[1]=0;

 sum[10]=0;

 sum[19]=0;

 sum[4]=0;

 sum[13]=0;

 sum[22]=0;

 sum[7]=0;

 sum[16]=0;

 sum[25]=0;

 sum[2]=0;

 sum[11]=0;

 sum[20]=0;

 sum[5]=0;

 sum[14]=0;

 sum[23]=0;

 sum[8]=0;

 sum[17]=0;

 sum[26]=0;
       
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