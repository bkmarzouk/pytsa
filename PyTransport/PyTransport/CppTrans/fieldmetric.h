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

#ifndef FIELDMETRIC_H  // Prevents the class being re-defined
#define FIELDMETRIC_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;


class fieldmetric
{
private:
	int nF; // field number
    int nP; // params number which definFs potential

public:
    fieldmetric()
   {
// #FP
nF=3;
nP=3;
	
   }
	
	
	//calculates fieldmetic()
	vector<double> fmetric(vector<double> f, vector<double> p)
	{
		vector<double> FM((2*nF)*(2*nF),0.0) ;
        
// metric

 FM[0]=1;

 FM[1]=0;

 FM[2]=0;

 FM[3]=1;

 FM[4]=0;

 FM[5]=0;

 FM[6]=0;

 FM[7]=1;

 FM[8]=0;

 FM[9]=0;

 FM[10]=1;

 FM[11]=0;

 FM[12]=0;

 FM[13]=0;

 FM[14]=1;

 FM[15]=0;

 FM[16]=0;

 FM[17]=1;

 FM[18]=1;

 FM[19]=0;

 FM[20]=0;

 FM[21]=1;

 FM[22]=0;

 FM[23]=0;

 FM[24]=0;

 FM[25]=1;

 FM[26]=0;

 FM[27]=0;

 FM[28]=1;

 FM[29]=0;

 FM[30]=0;

 FM[31]=0;

 FM[32]=1;

 FM[33]=0;

 FM[34]=0;

 FM[35]=1;

         return FM;
	}
	
	
	
	//calculates ChristoffelSymbole()
	vector<double> Chroff(vector<double> f, vector<double> p)
	{
		vector<double> CS((2*nF)*(2*nF)*(2*nF),0.0);
	
// Christoffel

 CS[0]=0;

 CS[1]=0;

 CS[2]=0;

 CS[3]=0;

 CS[4]=0;

 CS[5]=0;

 CS[6]=0;

 CS[7]=0;

 CS[8]=0;

 CS[9]=0;

 CS[10]=0;

 CS[11]=0;

 CS[12]=0;

 CS[13]=0;

 CS[14]=0;

 CS[15]=0;

 CS[16]=0;

 CS[17]=0;

 CS[18]=0;

 CS[19]=0;

 CS[20]=0;

 CS[21]=0;

 CS[22]=0;

 CS[23]=0;

 CS[24]=0;

 CS[25]=0;

 CS[26]=0;

 CS[27]=0;

 CS[28]=0;

 CS[29]=0;

 CS[30]=0;

 CS[31]=0;

 CS[32]=0;

 CS[33]=0;

 CS[34]=0;

 CS[35]=0;

 CS[36]=0;

 CS[37]=0;

 CS[38]=0;

 CS[39]=0;

 CS[40]=0;

 CS[41]=0;

 CS[42]=0;

 CS[43]=0;

 CS[44]=0;

 CS[45]=0;

 CS[46]=0;

 CS[47]=0;

 CS[48]=0;

 CS[49]=0;

 CS[50]=0;

 CS[51]=0;

 CS[52]=0;

 CS[53]=0;

 CS[54]=0;

 CS[55]=0;

 CS[56]=0;

 CS[57]=0;

 CS[58]=0;

 CS[59]=0;

 CS[60]=0;

 CS[61]=0;

 CS[62]=0;

 CS[63]=0;

 CS[64]=0;

 CS[65]=0;

 CS[66]=0;

 CS[67]=0;

 CS[68]=0;

 CS[69]=0;

 CS[70]=0;

 CS[71]=0;

 CS[72]=0;

 CS[73]=0;

 CS[74]=0;

 CS[75]=0;

 CS[76]=0;

 CS[77]=0;

 CS[78]=0;

 CS[79]=0;

 CS[80]=0;

 CS[81]=0;

 CS[82]=0;

 CS[83]=0;

 CS[84]=0;

 CS[85]=0;

 CS[86]=0;

 CS[87]=0;

 CS[88]=0;

 CS[89]=0;

 CS[90]=0;

 CS[91]=0;

 CS[92]=0;

 CS[93]=0;

 CS[94]=0;

 CS[95]=0;

 CS[96]=0;

 CS[97]=0;

 CS[98]=0;

 CS[99]=0;

 CS[100]=0;

 CS[101]=0;

 CS[102]=0;

 CS[103]=0;

 CS[104]=0;

 CS[105]=0;

 CS[106]=0;

 CS[107]=0;

 CS[108]=0;

 CS[109]=0;

 CS[110]=0;

 CS[111]=0;

 CS[112]=0;

 CS[113]=0;

 CS[114]=0;

 CS[115]=0;

 CS[116]=0;

 CS[117]=0;

 CS[118]=0;

 CS[119]=0;

 CS[120]=0;

 CS[121]=0;

 CS[122]=0;

 CS[123]=0;

 CS[124]=0;

 CS[125]=0;

 CS[126]=0;

 CS[127]=0;

 CS[128]=0;

 CS[129]=0;

 CS[130]=0;

 CS[131]=0;

 CS[132]=0;

 CS[133]=0;

 CS[134]=0;

 CS[135]=0;

 CS[136]=0;

 CS[137]=0;

 CS[138]=0;

 CS[139]=0;

 CS[140]=0;

 CS[141]=0;

 CS[142]=0;

 CS[143]=0;

 CS[144]=0;

 CS[145]=0;

 CS[146]=0;

 CS[147]=0;

 CS[148]=0;

 CS[149]=0;

 CS[150]=0;

 CS[151]=0;

 CS[152]=0;

 CS[153]=0;

 CS[154]=0;

 CS[155]=0;

 CS[156]=0;

 CS[157]=0;

 CS[158]=0;

 CS[159]=0;

 CS[160]=0;

 CS[161]=0;

 CS[162]=0;

 CS[163]=0;

 CS[164]=0;

 CS[165]=0;

 CS[166]=0;

 CS[167]=0;

 CS[168]=0;

 CS[169]=0;

 CS[170]=0;

 CS[171]=0;

 CS[172]=0;

 CS[173]=0;

 CS[174]=0;

 CS[175]=0;

 CS[176]=0;

 CS[177]=0;

 CS[178]=0;

 CS[179]=0;

 CS[180]=0;

 CS[181]=0;

 CS[182]=0;

 CS[183]=0;

 CS[184]=0;

 CS[185]=0;

 CS[186]=0;

 CS[187]=0;

 CS[188]=0;

 CS[189]=0;

 CS[190]=0;

 CS[191]=0;

 CS[192]=0;

 CS[193]=0;

 CS[194]=0;

 CS[195]=0;

 CS[196]=0;

 CS[197]=0;

 CS[198]=0;

 CS[199]=0;

 CS[200]=0;

 CS[201]=0;

 CS[202]=0;

 CS[203]=0;

 CS[204]=0;

 CS[205]=0;

 CS[206]=0;

 CS[207]=0;

 CS[208]=0;

 CS[209]=0;

 CS[210]=0;

 CS[211]=0;

 CS[212]=0;

 CS[213]=0;

 CS[214]=0;

 CS[215]=0;
        
		return CS;
	}
    

	
	// calculates RiemannTensor()
	vector<double> Riemn(vector<double> f, vector<double> p)
	{
		vector<double> RM((nF)*(nF)*(nF)*(nF),0.0);
		
// Riemann

 RM[0]=0;

 RM[1]=0;

 RM[2]=0;

 RM[3]=0;

 RM[4]=0;

 RM[5]=0;

 RM[6]=0;

 RM[7]=0;

 RM[8]=0;

 RM[9]=0;

 RM[10]=0;

 RM[11]=0;

 RM[12]=0;

 RM[13]=0;

 RM[14]=0;

 RM[15]=0;

 RM[16]=0;

 RM[17]=0;

 RM[18]=0;

 RM[19]=0;

 RM[20]=0;

 RM[21]=0;

 RM[22]=0;

 RM[23]=0;

 RM[24]=0;

 RM[25]=0;

 RM[26]=0;

 RM[27]=0;

 RM[28]=0;

 RM[29]=0;

 RM[30]=0;

 RM[31]=0;

 RM[32]=0;

 RM[33]=0;

 RM[34]=0;

 RM[35]=0;

 RM[36]=0;

 RM[37]=0;

 RM[38]=0;

 RM[39]=0;

 RM[40]=0;

 RM[41]=0;

 RM[42]=0;

 RM[43]=0;

 RM[44]=0;

 RM[45]=0;

 RM[46]=0;

 RM[47]=0;

 RM[48]=0;

 RM[49]=0;

 RM[50]=0;

 RM[51]=0;

 RM[52]=0;

 RM[53]=0;

 RM[54]=0;

 RM[55]=0;

 RM[56]=0;

 RM[57]=0;

 RM[58]=0;

 RM[59]=0;

 RM[60]=0;

 RM[61]=0;

 RM[62]=0;

 RM[63]=0;

 RM[64]=0;

 RM[65]=0;

 RM[66]=0;

 RM[67]=0;

 RM[68]=0;

 RM[69]=0;

 RM[70]=0;

 RM[71]=0;

 RM[72]=0;

 RM[73]=0;

 RM[74]=0;

 RM[75]=0;

 RM[76]=0;

 RM[77]=0;

 RM[78]=0;

 RM[79]=0;

 RM[80]=0;
     
        return RM;
	}

	// calculates RiemannTensor() covariant derivatives
	vector<double> Riemncd(vector<double> f, vector<double> p)
	{
		vector<double> RMcd((nF)*(nF)*(nF)*(nF)*(nF),0.0);
		
// Riemanncd

 RMcd[0]=0;

 RMcd[1]=0;

 RMcd[2]=0;

 RMcd[3]=0;

 RMcd[4]=0;

 RMcd[5]=0;

 RMcd[6]=0;

 RMcd[7]=0;

 RMcd[8]=0;

 RMcd[9]=0;

 RMcd[10]=0;

 RMcd[11]=0;

 RMcd[12]=0;

 RMcd[13]=0;

 RMcd[14]=0;

 RMcd[15]=0;

 RMcd[16]=0;

 RMcd[17]=0;

 RMcd[18]=0;

 RMcd[19]=0;

 RMcd[20]=0;

 RMcd[21]=0;

 RMcd[22]=0;

 RMcd[23]=0;

 RMcd[24]=0;

 RMcd[25]=0;

 RMcd[26]=0;

 RMcd[27]=0;

 RMcd[28]=0;

 RMcd[29]=0;

 RMcd[30]=0;

 RMcd[31]=0;

 RMcd[32]=0;

 RMcd[33]=0;

 RMcd[34]=0;

 RMcd[35]=0;

 RMcd[36]=0;

 RMcd[37]=0;

 RMcd[38]=0;

 RMcd[39]=0;

 RMcd[40]=0;

 RMcd[41]=0;

 RMcd[42]=0;

 RMcd[43]=0;

 RMcd[44]=0;

 RMcd[45]=0;

 RMcd[46]=0;

 RMcd[47]=0;

 RMcd[48]=0;

 RMcd[49]=0;

 RMcd[50]=0;

 RMcd[51]=0;

 RMcd[52]=0;

 RMcd[53]=0;

 RMcd[54]=0;

 RMcd[55]=0;

 RMcd[56]=0;

 RMcd[57]=0;

 RMcd[58]=0;

 RMcd[59]=0;

 RMcd[60]=0;

 RMcd[61]=0;

 RMcd[62]=0;

 RMcd[63]=0;

 RMcd[64]=0;

 RMcd[65]=0;

 RMcd[66]=0;

 RMcd[67]=0;

 RMcd[68]=0;

 RMcd[69]=0;

 RMcd[70]=0;

 RMcd[71]=0;

 RMcd[72]=0;

 RMcd[73]=0;

 RMcd[74]=0;

 RMcd[75]=0;

 RMcd[76]=0;

 RMcd[77]=0;

 RMcd[78]=0;

 RMcd[79]=0;

 RMcd[80]=0;

 RMcd[81]=0;

 RMcd[82]=0;

 RMcd[83]=0;

 RMcd[84]=0;

 RMcd[85]=0;

 RMcd[86]=0;

 RMcd[87]=0;

 RMcd[88]=0;

 RMcd[89]=0;

 RMcd[90]=0;

 RMcd[91]=0;

 RMcd[92]=0;

 RMcd[93]=0;

 RMcd[94]=0;

 RMcd[95]=0;

 RMcd[96]=0;

 RMcd[97]=0;

 RMcd[98]=0;

 RMcd[99]=0;

 RMcd[100]=0;

 RMcd[101]=0;

 RMcd[102]=0;

 RMcd[103]=0;

 RMcd[104]=0;

 RMcd[105]=0;

 RMcd[106]=0;

 RMcd[107]=0;

 RMcd[108]=0;

 RMcd[109]=0;

 RMcd[110]=0;

 RMcd[111]=0;

 RMcd[112]=0;

 RMcd[113]=0;

 RMcd[114]=0;

 RMcd[115]=0;

 RMcd[116]=0;

 RMcd[117]=0;

 RMcd[118]=0;

 RMcd[119]=0;

 RMcd[120]=0;

 RMcd[121]=0;

 RMcd[122]=0;

 RMcd[123]=0;

 RMcd[124]=0;

 RMcd[125]=0;

 RMcd[126]=0;

 RMcd[127]=0;

 RMcd[128]=0;

 RMcd[129]=0;

 RMcd[130]=0;

 RMcd[131]=0;

 RMcd[132]=0;

 RMcd[133]=0;

 RMcd[134]=0;

 RMcd[135]=0;

 RMcd[136]=0;

 RMcd[137]=0;

 RMcd[138]=0;

 RMcd[139]=0;

 RMcd[140]=0;

 RMcd[141]=0;

 RMcd[142]=0;

 RMcd[143]=0;

 RMcd[144]=0;

 RMcd[145]=0;

 RMcd[146]=0;

 RMcd[147]=0;

 RMcd[148]=0;

 RMcd[149]=0;

 RMcd[150]=0;

 RMcd[151]=0;

 RMcd[152]=0;

 RMcd[153]=0;

 RMcd[154]=0;

 RMcd[155]=0;

 RMcd[156]=0;

 RMcd[157]=0;

 RMcd[158]=0;

 RMcd[159]=0;

 RMcd[160]=0;

 RMcd[161]=0;

 RMcd[162]=0;

 RMcd[163]=0;

 RMcd[164]=0;

 RMcd[165]=0;

 RMcd[166]=0;

 RMcd[167]=0;

 RMcd[168]=0;

 RMcd[169]=0;

 RMcd[170]=0;

 RMcd[171]=0;

 RMcd[172]=0;

 RMcd[173]=0;

 RMcd[174]=0;

 RMcd[175]=0;

 RMcd[176]=0;

 RMcd[177]=0;

 RMcd[178]=0;

 RMcd[179]=0;

 RMcd[180]=0;

 RMcd[181]=0;

 RMcd[182]=0;

 RMcd[183]=0;

 RMcd[184]=0;

 RMcd[185]=0;

 RMcd[186]=0;

 RMcd[187]=0;

 RMcd[188]=0;

 RMcd[189]=0;

 RMcd[190]=0;

 RMcd[191]=0;

 RMcd[192]=0;

 RMcd[193]=0;

 RMcd[194]=0;

 RMcd[195]=0;

 RMcd[196]=0;

 RMcd[197]=0;

 RMcd[198]=0;

 RMcd[199]=0;

 RMcd[200]=0;

 RMcd[201]=0;

 RMcd[202]=0;

 RMcd[203]=0;

 RMcd[204]=0;

 RMcd[205]=0;

 RMcd[206]=0;

 RMcd[207]=0;

 RMcd[208]=0;

 RMcd[209]=0;

 RMcd[210]=0;

 RMcd[211]=0;

 RMcd[212]=0;

 RMcd[213]=0;

 RMcd[214]=0;

 RMcd[215]=0;

 RMcd[216]=0;

 RMcd[217]=0;

 RMcd[218]=0;

 RMcd[219]=0;

 RMcd[220]=0;

 RMcd[221]=0;

 RMcd[222]=0;

 RMcd[223]=0;

 RMcd[224]=0;

 RMcd[225]=0;

 RMcd[226]=0;

 RMcd[227]=0;

 RMcd[228]=0;

 RMcd[229]=0;

 RMcd[230]=0;

 RMcd[231]=0;

 RMcd[232]=0;

 RMcd[233]=0;

 RMcd[234]=0;

 RMcd[235]=0;

 RMcd[236]=0;

 RMcd[237]=0;

 RMcd[238]=0;

 RMcd[239]=0;

 RMcd[240]=0;

 RMcd[241]=0;

 RMcd[242]=0;
     
        return RMcd;
	}
    
    int getnF()
    {
        return nF;
    }
    


};
#endif

