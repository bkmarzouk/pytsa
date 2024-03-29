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




//model class file contains the defining features of the model -- the u1, u2 flow tesors, and the A B C tesors and u3 flow tensor as well as the guage transform N tensors
#ifndef MODEL_H  // Prevents the class being re-defined
#define MODEL_H 

#include "potential.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;


class model
{
private:
	int nF;  // field number
	int nP; // params number which definFs potential
	potential pot; // potential which defines model
    

public:
	// constructor
	model()
	{
	   potential pot;
        nP=pot.getnP();
        nF=pot.getnF();
    }

    double ke(vector<double> fdf, vector<double> p)
        {

        // Unpack velocities from fields-dotfields
        vector<double> phi(nF);
        vector<double> phidot(nF);

        for(int ii = 0; ii < nF; ii++)
        {
            phi[ii] = fdf[ii];
            phidot[ii] = fdf[ii + nF];
        }

        // Double-contract field space velocities with cov. field metric
		double dotsigmasq = 0;

        // Iterate over field indices
        for (int i = 0; i < nF; i++){

            // dotphi^i
            double pd = phidot[i];

            // Compute product
            dotsigmasq += pd * pd;
        }

        // Return field kinetic energy
		return 0.5 * dotsigmasq;
    }

    // function returns Hubble rate
	double H(vector<double> f, vector<double> p   )
	{
	    // Compute potential
		double Vi;
		Vi=pot.V(f,p);

        // Compute field-space kinetic energy
        double dotsigmasq = ke(f, p);

        // Compute hubble rate directly from background EoM
		double hubble = sqrt((dotsigmasq + Vi)/3.);

		return hubble;
	}


    // function returns H dot
    double Hdot(vector<double> fdf, vector<double> p)
	{
        // Hdot is simply minus the kinetic energy of the fields
        double dotsigmasq = ke(fdf, p);

		return -1. * dotsigmasq;
	}

	// function returns H dot dot
	double Hddot(vector<double> fdf, vector<double> p)
	{
	    // We can get the second derivative from the background equations

        // Unpack fieldsdotfields
        vector<double> f(nF);
        vector<double> v(nF);

        for(int ii = 0; ii < nF; ii++){
            f[ii] = fdf[ii];
            v[ii] = fdf[ii + nF];
        }

	    // first derivative of potential
	    vector<double> dv = pot.dV(fdf, p);

	    // hubble rate
	    double hubble = H(fdf, p);

	    double out = 0;

	    for (int ii=0; ii < nF; ii++){
	        out +=  v[ii] * (3 * hubble * v[ii] + dv[ii]);
	    }

	    return out;
	}


    // function returns a double dot
    double addot(vector<double> f, vector<double> p)
    {
        double sum=0.;
        double addot;
        double Vi;
        Vi=pot.V(f,p);
        for(int i=0; i<nF; i++)
        {
            sum= sum - 1./2.*(f[nF+i]*f[nF+i]);
        }
        addot=-1./3.*(sum+Vi);
        return addot;
    }
    
    // function returns epsilon
	double Ep(vector<double> f,vector<double> p)
	{
        // Compute potential and Hubble rate
        double Vi = pot.V(f, p);
        double Hubble = H(f, p);

        // Compute epsilon via definition
        return 3. - Vi / (Hubble*Hubble);
	}

	// funtion returns eta
	double Eta(vector<double> fdf, vector<double> p)
	{

        // Unpack fieldsdotfields
        vector<double> f(nF);
        vector<double> v(nF);

        for(int ii = 0; ii < nF; ii++){
            f[ii] = fdf[ii];
            v[ii] = fdf[ii + nF];
        }


        double hub = H(f, p);
        double dothub = Hdot(fdf, p);
        double ddothub = Hddot(fdf, p);

        double t1 = ddothub / (dothub * hub);
        double t2 = 2 * dothub / (hub * hub);

        return t1 - t2;
	}

    // function returns mass-squared-matrix
    vector<double> Mij(vector<double> fdf, vector<double> p, bool hessianApprox, bool covariantExpression)
    {
        // Retains covariant output kwarg for signature consistency with NC

        // Unpack fieldsdotfields
        vector<double> f(nF);
        vector<double> v(nF);

        for(int ii = 0; ii < nF; ii++){
            f[ii] = fdf[ii];
            v[ii] = fdf[ii + nF];
        }

        // Potential derivatives
        vector<double> ddVi;
        vector<double> dVi;
        ddVi = pot.dVV(fdf,p);
		dVi =  pot.dV(fdf,p);

        // Get epsilon
        double eps = Ep(fdf, p);

        // Get hubble rate
        double Hi = H(fdf, p);

		// double time derivatives of fields (just background eom)
		vector<long double> cdt_phidot_d(nF);
		for (int ii = 0; ii < nF; ii++)
		{
            cdt_phidot_d[ii] = -3.*Hi*v[ii] - dVi[ii];
		}

        // Output
        vector<double> mijout(nF*nF);


		for (int ii=0; ii < nF; ii++)
		{

		    for (int jj=0; jj < nF; jj++)
		    {

                double hij_sum = 0;  // Hessian term
                double kij_sum = 0;  // Kinetic term

		        hij_sum += ddVi[ii*nF + jj];

		        if (hessianApprox == false)
		        {

                    kij_sum -= (v[ii] * cdt_phidot_d[jj] + v[jj] * cdt_phidot_d[ii]) / Hi;
                    kij_sum -= (3. - eps) * v[ii] * v[jj];

                }

                mijout[ii*nF + jj] = hij_sum + kij_sum;

		    }
		}

        for (int ii=0; ii < nF; ii++){
            for (int jj=0; jj < nF; jj++){
                int mij_idx = ii*nF + jj;
                mijout[mij_idx] /= Hi*Hi;
            }
        }

        return mijout;
    }

    // function returns number of fields
    int getnF()
    {
        return nF;
    }
    
    // function returns number of fields
    int getnP()
    {
        return nP;
    }
	
    
    // a rescale function for delta dot phi to improve performance
    double scale(vector<double> f, vector<double> p,  double N)
    {
        double k=1.0;
        double a =exp(N);
        double Hi=H(f,p);
    //    return a;
    //    return a + (1.0-a)/pow(1. + k/a/Hi,2.0);
        return  a/(1.+a*Hi/k)/Hi; //+ 1./pow(1. + k/a/Hi,10.0);
    }
    
    // derivative of rescaling function divided by rescaling function
    double dscale(vector<double> f, vector<double> p, double N)
    {
        double k = 1.0;
        double a = exp(N);
        double Hi = H(f,p);
        double Hdi = Hdot(f, p);
        //return a*Hi;
       // return  a*Hi - a *Hi/pow(1. + k/a/Hi,2.0) + (1.0-a)*2.0*k/a/a/Hi/Hi*(a*Hi*Hi+a*Hdi)/pow(1. + k/a/Hi,3.0);
        return  -Hdi/Hi/Hi*a/(1.+a*Hi/k) + a/(1.+a*Hi/k) -a*(a*Hi*Hi/k + a*Hdi/k)/(1.+a*Hi/k)/(1.+a*Hi/k)/Hi; //+ 10.0*k/a/a/Hi/Hi*(a*Hi*Hi+a*Hdi)/pow(1. + k/a/Hi,11.0);
    }
    
    // calculates u1
	vector<double> u(vector<double> f,vector<double> p)
	{
		vector<double> u1out(2*nF);
		vector<double> dVi;
		double Hi;
		Hi=H(f,p);
		
		for(int i=0; i<nF; i++)
		{
			u1out[i]  = f[nF+i]/Hi;
		}
		
		dVi=pot.dV(f,p);

		for(int i=0; i<nF; i++)	
		{
			u1out[nF+i]  = -3.*Hi*f[nF+i]/Hi-dVi[i]/Hi;
		}
		return u1out;
	}

	// calculates u2
	vector<double> u(vector<double> f,vector<double> p, double k1, double N)
	{
		double a = exp(N);
		double ep = Ep(f,p);
		vector<double> u2out(2*nF*2*nF);
		double Hi=H(f,p);
        double s=scale(f,p,N);
        double ds=dscale(f,p,N);
		vector<double> dVVi;
		dVVi = pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);

		for(int i = 0; i<nF; i++){for(int j = 0; j<nF; j++){
            u2out[i+ j*2*nF]=0.;
            u2out[i+(j+nF)*2*nF]=0.;
            u2out[i+nF+(j)*2*nF]=(-dVVi[i + nF*j] + (-3.+ep)*f[nF+i]*f[nF+j] + 1./Hi*(-dVi[i])*f[nF+j] + 1./Hi*f[nF+i]*(-dVi[j]) )/Hi *s; // *a;
            u2out[i+nF+(j+nF)*2*nF]=0.;
            if(i==j){
                u2out[i+nF+(j)*2*nF]=u2out[i+nF+(j)*2*nF]-1.0*(k1*k1)/(a*a)/Hi  * s ;// *a;
                u2out[i+(j+nF)*2*nF]=u2out[i+(j+nF)*2*nF] + 1./Hi  /s;// /a ;
                      u2out[i+nF+(j+nF)*2*nF]= u2out[i+nF+(j+nF)*2*nF] - 3.0*Hi/Hi  + ds/s/Hi; // - 2.0*Hi/Hi;
            }
        }}

		return u2out;
	}

    
    //calculates A (the field field field term of action)
    vector<double> Acalc(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
		double a = exp(N);
        double Vi=pot.V(f,p);
		double Hi=H(f,p);
     
        vector<double> dVVi;
		dVVi=pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		vector<double> dVVVi;
		dVVVi=pot.dVVV(f,p);
        vector<double> Xi(nF); vector<double> A(nF*nF*nF);
        
        double sum1=0;
		for(int i=0;i<nF;i++){sum1=sum1+f[nF+i]*f[nF+i];}
		for(int i=0;i<nF;i++){Xi[i] = 2.*(-dVi[i]-3.*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			A[i + j*nF +k* nF*nF] = -1./3. * dVVVi[i + j*nF +k* nF*nF]
			- 1./3.*f[nF + i]/2./Hi* dVVi[j + k*nF]
            - 1./3.*f[nF + j]/2./Hi* dVVi[i + k*nF]
            - 1./3.*f[nF + k]/2./Hi* dVVi[i + j*nF]
			+ 1./3.*f[nF + i] * f[nF + j]/8./Hi/Hi * Xi[k]
            + 1./3.*f[nF + i] * f[nF + k]/8./Hi/Hi * Xi[j]
            + 1./3.*f[nF + k] * f[nF + j]/8./Hi/Hi * Xi[i]
			+ 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xi[j] *Xi[k]
            + 1./3.*f[nF + j]/32./Hi/Hi/Hi * Xi[i] *Xi[k]
            + 1./3.*f[nF + k]/32./Hi/Hi/Hi * Xi[i] *Xi[j]
			+ 1.*f[nF + i]*f[nF + j]*f[nF + k]/8./Hi/Hi/Hi*2.*Vi
			- 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xi[j] * Xi[k] * (k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            - 1./3.*f[nF + j]/32./Hi/Hi/Hi * Xi[i] * Xi[k] * (k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.
            - 1./3.*f[nF + k]/32./Hi/Hi/Hi * Xi[i] * Xi[j] * (k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4.;
    		if(j==k){A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*f[nF+i]/2./Hi*(-k2*k2-k3*k3+k1*k1)/a/a/2.;}
			if(i==k){A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*f[nF+j]/2./Hi*(-k1*k1-k3*k3+k2*k2)/a/a/2.;}
			if(i==j){A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*f[nF+k]/2./Hi*(-k2*k2-k1*k1+k3*k3)/a/a/2.;}
            }}}

        return A;
    }

    //calculates AS (the "slow" parts of the field field field term of action -- this is used only for initial conditions)
    vector<double> AScalc(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
    {
        double Vi=pot.V(f,p);
        double Hi=H(f,p);
        
        vector<double> dVVi;
        dVVi=pot.dVV(f,p);
        vector<double> dVi;
        dVi =  pot.dV(f,p);
        vector<double> dVVVi;
        dVVVi=pot.dVVV(f,p);
        vector<double> Xi(nF); vector<double> AS(nF*nF*nF);
        
        double sum1=0;
        for(int i=0;i<nF;i++){sum1=sum1+f[nF+i]*f[nF+i];}
        for(int i=0;i<nF;i++){Xi[i] = 2.*(-dVi[i]-3.*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
        
        for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
            AS[i + j*nF +k* nF*nF] = -1./3. * dVVVi[i + j*nF +k* nF*nF]
            - 1./3.*f[nF + i]/2./Hi* dVVi[j + k*nF]
            - 1./3.*f[nF + j]/2./Hi* dVVi[i + k*nF]
            - 1./3.*f[nF + k]/2./Hi* dVVi[i + j*nF]
            + 1./3.*f[nF + i] * f[nF + j]/8./Hi/Hi * Xi[k]
            + 1./3.*f[nF + i] * f[nF + k]/8./Hi/Hi * Xi[j]
            + 1./3.*f[nF + k] * f[nF + j]/8./Hi/Hi * Xi[i]
            + 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xi[j] *Xi[k]
            + 1./3.*f[nF + j]/32./Hi/Hi/Hi * Xi[i] *Xi[k]
            + 1./3.*f[nF + k]/32./Hi/Hi/Hi * Xi[i] *Xi[j]
            + 1.*f[nF + i]*f[nF + j]*f[nF + k]/8./Hi/Hi/Hi*2.*Vi
            - 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xi[j] * Xi[k] * (k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            - 1./3.*f[nF + j]/32./Hi/Hi/Hi * Xi[i] * Xi[k] * (k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.
            - 1./3.*f[nF + k]/32./Hi/Hi/Hi * Xi[i] * Xi[j] * (k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4.;
        }}}
        
        return AS;
    }
    
    
    //Calculates B term of action
   vector<double> Bcalc(vector<double> f,vector<double> p, double k1, double k2, double k3,double N)
	{
		
        double Hi=H(f,p);
        
		vector<double> dVVi;
       // dVVi = new double[nF*nF];
		dVVi=pot.dVV(f,p);
		vector<double> dVi; //dVi = new double[nF];
		dVi =  pot.dV(f,p);
		vector<double> dVVVi; //dVVVi = new double[nF*nF*nF];
		dVVVi=pot.dVVV(f,p);
        vector<double> Xi(nF);vector<double> B(nF*nF*nF);
		
        double sum1=0;
		for(int i=0;i<nF;i++){sum1=sum1+f[nF+i]*f[nF+i];}
		for(int i=0;i<nF;i++){Xi[i] = 2.0*(-dVi[i]-3.0*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
		
        
        for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			B[i + j*nF +k* nF*nF] = 1.*f[nF + i]*f[nF+j]*f[nF+k]/4./Hi/Hi
			- 1./2.*f[nF + i] * f[nF + k]/8./Hi/Hi/Hi * Xi[j]
            - 1./2.*f[nF + j] * f[nF + k]/8./Hi/Hi/Hi * Xi[i]
			+ 1./2.*f[nF + i] * f[nF + k]/8./Hi/Hi/Hi * Xi[j]*(k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            + 1./2.*f[nF + j] * f[nF + k]/8./Hi/Hi/Hi * Xi[i]*(k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.;
			if(j==k){B[i + j*nF +k* nF*nF] = B[i + j*nF +k* nF*nF] - 1.*Xi[i]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k1/k1/2.;}
			if(i==k){B[i + j*nF +k* nF*nF] = B[i + j*nF +k* nF*nF] - 1.*Xi[j]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k2/k2/2.;}
		}}}
        return B;
    }

    //Calculates C term of action
    vector<double> Ccalc(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
		double Hi=H(f,p);
        
     	vector<double> dVVi; //dVVi = new double[nF*nF];
		dVVi=pot.dVV(f,p);
		vector<double> dVi; //dVi = new double[nF];
		dVi =  pot.dV(f,p);
		vector<double> dVVVi; //dVVVi = new double[nF*nF*nF];
		dVVVi=pot.dVVV(f,p);
        vector<double> Xi(nF); vector<double> C(nF*nF*nF);
		
        double sum1=0;
		for(int i=0;i<nF;i++){sum1=sum1+f[nF+i]*f[nF+i];}
		for(int i=0;i<nF;i++){Xi[i] = 2.*(-dVi[i]-3.*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			C[i + j*nF +k* nF*nF] = 1.*f[nF + i]*f[nF+j]*f[nF+k]/8./Hi/Hi/Hi
			- 1.*f[nF + i] * f[nF+j] *f[nF+k]/8./Hi/Hi/Hi *(k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4. ;
			if(i==j){C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] - 1.*f[nF+k]/2./Hi;}
			if(j==k){C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] + f[nF+i]/2./Hi*(-k1*k1-k3*k3+k2*k2)/k1/k1/2.;}
			if(i==k){C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] + f[nF+j]/2./Hi*(-k2*k2-k3*k3+k1*k1)/k2/k2/2.;}
		}}}
        return C;
    }
    

    
	//calculates u3
	vector<double> u(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
        vector<double>  A, B,B2, B3, C, C2,C3;
        double Hi;
		Hi=H(f,p);
        double s=scale(f,p,N);
        
        A = Acalc(f,p, k1, k2, k3 ,N);
        B= Bcalc(f,p, k2, k3, k1 ,N);
        B2=  Bcalc(f,p, k1, k2, k3 ,N);
        B3=Bcalc(f,p, k1, k3, k2 ,N);
        C=  Ccalc(f,p, k1, k2, k3 ,N);
        C2=  Ccalc(f,p, k1, k3, k2 ,N);
        C3 = Ccalc(f,p, k3, k2, k1 ,N);

        vector<double> u3out(2*nF*2*nF*2*nF);
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			u3out[i+j*2*nF+k*2*nF*2*nF]= -B[j+k*nF+i*nF*nF]/Hi;
			
            u3out[(i)+(nF+j)*2*nF+k*2*nF*2*nF]= -C[i+j*nF+k*nF*nF]/Hi  /s; // /a;
            u3out[(i)+j*2*nF+(k+nF)*2*nF*2*nF]= -C2[i+k*nF+j*nF*nF]/Hi /s;// /a;
			
			u3out[(i)+(j+nF)*2*nF+(k+nF)*2*nF*2*nF]= 0.;
			
            u3out[(nF+i) + j*2*nF + k*2*nF*2*nF]= 3.*A[i+j*nF+k*nF*nF]/Hi  *s;// *a;
		
			u3out[(nF+i)+(nF+j)*2*nF+k*2*nF*2*nF]=B3[i+k*nF+j*nF*nF]/Hi ;
			u3out[(nF+i)+(j)*2*nF+(k+nF)*2*nF*2*nF]=B2[i+j*nF+k*nF*nF]/Hi ;
			
            u3out[(nF+i)+(j+nF)*2*nF + (k+nF)*2*nF*2*nF]=C3[k+j*nF+i*nF*nF]/Hi  /s;// /a;

		}}}
        return u3out;
	}


//calculates N1
vector<double> N1(vector<double> f,vector<double> p, double N)
{
    double Hd=Hdot(f,p);
    double Hi=H(f,p);
    //double a = exp(N);
    vector<double> dVi;
    vector<double> Ni(2*nF);
    dVi=pot.dV(f,p);
 
    for(int i=0;i<nF;i++){
        Ni[i] = 1./2.*Hi/Hd * f[nF+i];
     
        Ni[nF+i] = 0. ;
        }

    return Ni;
}

vector<double> N2(vector<double> f, vector<double> p, double k1, double k2, double k3, double N)
{
    double Hd=Hdot(f, p);
    double Hin=H(f,p);
    vector<double> dVi, dVVi;
    vector<double> Nii(2*nF*2*nF);
    double s = scale(f,p,N);
    dVi=pot.dV(f,p);
    dVVi=pot.dVV(f,p);
    
    double sum3 = 0.0;
    for(int i=0;i<nF;i++){sum3=sum3+dVi[i]*f[nF+i]/Hin/Hin/Hin;}
     
     
    double ep = -Hd/Hin/Hin;
    for(int i=0;i<nF;i++){for(int j=0; j<nF; j++){
    Nii[i + (j) * 2*nF]= 2./ep/Hin/Hin/6. * (f[nF+i]*f[nF+j] *(-3./2. + 9./2./ep + 3./4.*sum3/ep/ep));
    Nii[i + (j+nF) * 2*nF]=2./ep/Hin/Hin/6.*3./2.*f[i+nF]*f[j+nF]/Hin/ep  /s;// /a;
    Nii[i+nF + (j) * 2*nF]=2./ep/Hin/Hin/6.*3./2.*f[i+nF]*f[j+nF]/Hin/ep  /s; // /a;
    Nii[i+nF + (j+nF) * 2*nF]=0.;
        if(i==j){Nii[i+nF+(j)*2*nF] = Nii[i+nF + (j) * 2*nF] - 2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k3*k3)  /s;// /a ;
                 Nii[i+(j+nF)*2*nF] = Nii[i + (j+nF) * 2*nF] - 2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k2*k2)  /s;}// /a;}
    }}

    return Nii;

}





};
#endif
