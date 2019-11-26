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
// Potential file rewriten at Tue Nov 26 17:08:17 2019

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
nF=5;
nP=5;

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
  sum=(1.0/2.0)*std::pow(f[0], 2)*std::pow(p[0], 2) + (1.0/2.0)*std::pow(f[1], 2)*std::pow(p[1], 2) + (1.0/2.0)*std::pow(f[2], 2)*std::pow(p[2], 2) + (1.0/2.0)*std::pow(f[3], 2)*std::pow(p[3], 2) + (1.0/2.0)*std::pow(f[4], 2)*std::pow(p[4], 2);
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

 sum[3]=f[3]*std::pow(p[3], 2);

 sum[4]=f[4]*std::pow(p[4], 2);
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot

 sum[0]=std::pow(p[0], 2);

 sum[5]=0;

 sum[10]=0;

 sum[15]=0;

 sum[20]=0;

 sum[1]=0;

 sum[6]=std::pow(p[1], 2);

 sum[11]=0;

 sum[16]=0;

 sum[21]=0;

 sum[2]=0;

 sum[7]=0;

 sum[12]=std::pow(p[2], 2);

 sum[17]=0;

 sum[22]=0;

 sum[3]=0;

 sum[8]=0;

 sum[13]=0;

 sum[18]=std::pow(p[3], 2);

 sum[23]=0;

 sum[4]=0;

 sum[9]=0;

 sum[14]=0;

 sum[19]=0;

 sum[24]=std::pow(p[4], 2);
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot

 sum[0]=0;

 sum[25]=0;

 sum[50]=0;

 sum[75]=0;

 sum[100]=0;

 sum[5]=0;

 sum[30]=0;

 sum[55]=0;

 sum[80]=0;

 sum[105]=0;

 sum[10]=0;

 sum[35]=0;

 sum[60]=0;

 sum[85]=0;

 sum[110]=0;

 sum[15]=0;

 sum[40]=0;

 sum[65]=0;

 sum[90]=0;

 sum[115]=0;

 sum[20]=0;

 sum[45]=0;

 sum[70]=0;

 sum[95]=0;

 sum[120]=0;

 sum[1]=0;

 sum[26]=0;

 sum[51]=0;

 sum[76]=0;

 sum[101]=0;

 sum[6]=0;

 sum[31]=0;

 sum[56]=0;

 sum[81]=0;

 sum[106]=0;

 sum[11]=0;

 sum[36]=0;

 sum[61]=0;

 sum[86]=0;

 sum[111]=0;

 sum[16]=0;

 sum[41]=0;

 sum[66]=0;

 sum[91]=0;

 sum[116]=0;

 sum[21]=0;

 sum[46]=0;

 sum[71]=0;

 sum[96]=0;

 sum[121]=0;

 sum[2]=0;

 sum[27]=0;

 sum[52]=0;

 sum[77]=0;

 sum[102]=0;

 sum[7]=0;

 sum[32]=0;

 sum[57]=0;

 sum[82]=0;

 sum[107]=0;

 sum[12]=0;

 sum[37]=0;

 sum[62]=0;

 sum[87]=0;

 sum[112]=0;

 sum[17]=0;

 sum[42]=0;

 sum[67]=0;

 sum[92]=0;

 sum[117]=0;

 sum[22]=0;

 sum[47]=0;

 sum[72]=0;

 sum[97]=0;

 sum[122]=0;

 sum[3]=0;

 sum[28]=0;

 sum[53]=0;

 sum[78]=0;

 sum[103]=0;

 sum[8]=0;

 sum[33]=0;

 sum[58]=0;

 sum[83]=0;

 sum[108]=0;

 sum[13]=0;

 sum[38]=0;

 sum[63]=0;

 sum[88]=0;

 sum[113]=0;

 sum[18]=0;

 sum[43]=0;

 sum[68]=0;

 sum[93]=0;

 sum[118]=0;

 sum[23]=0;

 sum[48]=0;

 sum[73]=0;

 sum[98]=0;

 sum[123]=0;

 sum[4]=0;

 sum[29]=0;

 sum[54]=0;

 sum[79]=0;

 sum[104]=0;

 sum[9]=0;

 sum[34]=0;

 sum[59]=0;

 sum[84]=0;

 sum[109]=0;

 sum[14]=0;

 sum[39]=0;

 sum[64]=0;

 sum[89]=0;

 sum[114]=0;

 sum[19]=0;

 sum[44]=0;

 sum[69]=0;

 sum[94]=0;

 sum[119]=0;

 sum[24]=0;

 sum[49]=0;

 sum[74]=0;

 sum[99]=0;

 sum[124]=0;
       
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