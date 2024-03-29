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


// C++ file which defines the functions make available to Python through the MTeasy module.
#include <Python.h>
#include <iostream>
#include "numpy/arrayobject.h"

//don't adjust the labels at the end of the 4 lines below (they are used to fix directory structure)
#include"CPPT_DIR/evolve.h"//evolve
#include"CPPT_DIR/moments.h"//moments
#include"CPPT_DIR/model.h"//model
#include"STEP_DIR/rkf45.hpp"//stepper
//*************************************************************************************************

#include <math.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <time.h>
# include <iomanip>
# include <cmath>

using namespace std;

// The line below is updated evey time the moduleSetup file is run.
// Package recompile attempted at: Wed Nov  3 15:56:05 2021


// Look for environment variable that indicates we should be running in sampler mode, i.e. flag returns
bool _SamplerMode(){
    char const* tmp = getenv( "SAMPLER_MODE" );

    bool sampler_mode;

    if ( tmp == NULL ) {
        sampler_mode = false;
    }
    else {
        std::string s( tmp );
        if (s == "TRUE") {
            sampler_mode = true;
        }
        else {
            sampler_mode = false;
        }
    }
    return sampler_mode;
}

// Changes python array into C array (or rather points to pyarray data)
//    Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carray(PyArrayObject *arrayin)
{
    return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

int size_pyvector(PyArrayObject *arrayin)
{
      return arrayin->dimensions[0];  /* pointer to arrayin data as double */
}

// function to retun amplitude of potential
static PyObject* MT_V(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *params;
    double *Cfields,*Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);
    potential pp;
    int nF = pp.getnF(); if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}

    vector<double> vectIn;
    vectIn = vector<double>(Cfields, Cfields + nF);

    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    return Py_BuildValue("d", pp.V(vectIn,Vparams));
}

// function to calculate derivatives of potential
static PyObject* MT_dV(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *dVI, *params;
    double *Cfields, *dVC, *Cparams ;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);

    potential pp;
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}
    vector<double> vectIn;
    npy_intp dims[1];
    dims[0]=nF;

    dVI = (PyArrayObject*) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    dVC = (double*) dVI->data;

    vectIn = vector<double>(Cfields, Cfields +  nF);

    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    vector<double> dVect = pp.dV(vectIn,Vparams);
    for(int i=0; i<nF;i++){dVC[i] = dVect[i];}

    return PyArray_Return(dVI);
}

// function to calculate derivatives of potential
static PyObject* MT_ddV(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *ddVI, *params;
    double *Cfields, *ddVC, *Cparams ;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);

    potential pp;
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}

    vector<double> vectIn;
    npy_intp dims[2];
    dims[0]=nF; dims[1]=nF;

    ddVI = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    ddVC = (double*) ddVI->data;

    vectIn = vector<double>(Cfields, Cfields +  nF);

    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    vector<double> ddVect = pp.dVV(vectIn,Vparams);
    for(int i=0; i<nF;i++){for(int j=0; j<nF;j++){ddVC[i+j*nF] = ddVect[i+j*nF];}}

    return PyArray_Return(ddVI);
}

// function to calculate Hubble rate
static PyObject* MT_H(PyObject* self,  PyObject *args)
{
    PyArrayObject *fields_dfieldsIn, *params;
    double *Cfields_dfields, *Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fields_dfieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields_dfields = pyvector_to_Carray(fields_dfieldsIn);
    model mm;
    int nF = mm.getnF(); if (2*nF!=size_pyvector(fields_dfieldsIn)){cout<< "\n \n \n field space array not of correct length\n \n \n ";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);
    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    return Py_BuildValue("d", mm.H(vectIn, Vparams));
}

// function to return number of fields (useful for cross checks)
static PyObject* MT_fieldNumber(PyObject* self)
{
    model mm;
    return Py_BuildValue("i",mm.getnF());
}

// function to return number of parameters (useful for cross checks)
static PyObject* MT_paramNumber(PyObject* self)
{
    model mm;
    return Py_BuildValue("i",mm.getnP());
}


// Compute slow roll epsilon param
static PyObject* MT_Ep(PyObject* self, PyObject *args)
{
    PyArrayObject *fields_dfieldsIn, *params;

    double *Cfields_dfields, *Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fields_dfieldsIn,&PyArray_Type,&params)) {
        return NULL;}

    Cfields_dfields = pyvector_to_Carray(fields_dfieldsIn);
    model mm;

    int nF = mm.getnF(); if (2*nF!=size_pyvector(fields_dfieldsIn)){cout<< "\n \n \n field space array not of correct length\n \n \n ";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);

    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    return Py_BuildValue("d", mm.Ep(vectIn, Vparams));

}

// Compute slow roll epsilon param
static PyObject* MT_Eta(PyObject* self, PyObject *args)
{
    PyArrayObject *fields_dfieldsIn, *params;

    double *Cfields_dfields, *Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fields_dfieldsIn,&PyArray_Type,&params)) {
        return NULL;}

    Cfields_dfields = pyvector_to_Carray(fields_dfieldsIn);
    model mm;

    int nF = mm.getnF(); if (2*nF!=size_pyvector(fields_dfieldsIn)){cout<< "\n \n \n field space array not of correct length\n \n \n ";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);

    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    return Py_BuildValue("d", mm.Eta(vectIn, Vparams));

}


// evaluates the full mass-matrix of the model subject to params + fieldsdotfields
static PyObject* MT_massMatrix(PyObject* self, PyObject *args)
{
    // Initialize python objects
    PyArrayObject *fieldsdotfields;
    PyArrayObject *params;
    PyArrayObject *mij_out;

    // Initialize corresponding cpp objects
    double *Cfieldsdotfields, *Cparams;

    bool hessapprox = false;
    bool covexp = false;

    // Demand python array objects for fields, dotfields and params
    if(!PyArg_ParseTuple(args, "O!O!|bb", &PyArray_Type, &fieldsdotfields, &PyArray_Type, &params,
                         &hessapprox, &covexp))
    {
        return NULL;
    }

    // Check input arrays are the correct length
    model mm;
    int nF = mm.getnF();

    if(size_pyvector(fieldsdotfields) != 2*nF)
      {
        cout << "\n \n \n fieldsdotfields array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // convert parameter list to a C array
    Cparams = pyvector_to_Carray(params);
    int nP = mm.getnP();

    if(size_pyvector(params) != nP)
      {
        cout << "\n \n \n parameters array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // convert fieldsdotfields to Cvectors

    // First convert python object to c array
    Cfieldsdotfields = pyvector_to_Carray(fieldsdotfields);

    // Create vectors from c arrays
    vector<double> Cvec_fieldsdotfields;
    vector<double> Cvec_params;
    Cvec_fieldsdotfields = vector<double>(Cfieldsdotfields, Cfieldsdotfields + 2*nF);
    Cvec_params = vector<double>(Cparams, Cparams + nP);

    // Computing mass matrix C vector
    vector<double> Cmij_flat = mm.Mij(Cvec_fieldsdotfields, Cvec_params, hessapprox, covexp);

    // Build 2 d array
    double Cmij[nF][nF];

    for (int ii = 0; ii < nF; ii++)
    {
        for (int jj = 0; jj < nF; jj++)
        {
            Cmij[ii][jj] = Cmij_flat[ii*nF + jj];
        }
    }

    // Define dimensions of output array (nF*nF)
    npy_intp dims[2];
    dims[0] = nF;
    dims[1] = nF;
    mij_out = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);

    // Create cpp pointer to pyarray obj
    double * Cmij_out;
    Cmij_out = (double *) PyArray_DATA(mij_out);


    for (int ii = 0; ii < nF; ii++)
    {
        for (int jj = 0; jj < nF; jj++)
        {
            Cmij_out[ii*nF + jj] = Cmij[ii][jj]; /// (Hub*Hub);
        }
    }



    return PyArray_Return(mij_out);
}


// detect end of inflation within a suitable search window
static PyObject* MT_findEndOfInflation(PyObject* self, PyObject* args)
{
    // initialize python array objects: initial conditions, parameters & tols
    PyArrayObject* initialCs;
    PyArrayObject* params;
    PyArrayObject* tols;

    // Set Nstart, max number of efoldings for search window and max integration times
    double Ninit     = 0.0;     // value of N corresponding to initial conditions
    double DeltaN    = 3000.0;  // default number of e-folds to search through
    double tmax_feoi = -1;      // max integration time, tmax_feoi < 0 => no limit

    // Setup integration flags
    int timeoutFlag     = -33;        // Integration timeout flag
    int integrationFlag = -32;        // Integration failure flag
    int eternalFlag     = -35;        // Unable to find end of inflation flag

    // Indicate whether to return flag tuple, or eFold at failure
    bool flagReturn = _SamplerMode();

    // parse arguments; args after | are optional (https://docs.python.org/3/c-api/arg.html)
    if(!PyArg_ParseTuple(args, "O!O!O!|dddb", &PyArray_Type, &initialCs, &PyArray_Type, &params, &PyArray_Type, &tols,
                         &Ninit, &DeltaN, &tmax_feoi, &flagReturn))
                         {
                            return NULL;
                         }

    // convert requested tolerances to a C array and extract absolute & relative error targets
    double* tolsC = pyvector_to_Carray(tols);
    double abserr = 1E-8;
    double relerr = 1E-8;

    // check input tols are formatted correctly, revert to default if required
    if(size_pyvector(tols) != 2) {
        cout << "\n \n \n incorrect tolerances input, using defaults \n \n \n";
      }
    else {
        abserr = tolsC[0];
        relerr = tolsC[1];
      }

    // convert initial conditions to a C array
    double* CinitialCs = pyvector_to_Carray(initialCs);

    // check whether the expected number of initial conditions were supplied
    model mm;
    int nF = mm.getnF();

    if(size_pyvector(initialCs) != 2*nF)
      {
        cout << "\n \n \n field space array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // convert parameter list to a C array
    double* Cparams = pyvector_to_Carray(params);
    int nP = mm.getnP();

    if(size_pyvector(params) != nP)
      {
        cout << "\n \n \n parameters array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // allocate working space for the stepper
    double* y  = new double[2*nF];      // current values
    double* dy = new double[2*nF];      // derivatives

    // initialize y using supplied initial conditions
    for(int i = 0; i < 2*nF; ++i)
      {
        y[i] = CinitialCs[i];
      }

    // populate dy for initial step
    double N = Ninit;
    const double Nstop = Ninit + DeltaN;
    evolveB(N, y, dy, Cparams);

    // Initialize time difference (will measure elapsed time with every step of integration)
    double deltat;

    // Build time instances
    time_t t_ini, t_fin;

    // Record current time before proc.
    time(&t_ini);

    // integrate background until we encounter the end-of-inflation, or the end of the search window
    int flag = -1;      // '-1' puts the integrator into 'single-step' mode, so it returns after taking one stride
    while(N < Nstop)
      {

        // clock time at start of loop
        time(&t_fin);

        // measure time elapsed since first time instance
        deltat = difftime(t_fin,t_ini);

        // if deltat exceeds allowed time and tmax_feoi > 0

        if((deltat > tmax_feoi) && (tmax_feoi>0)){

                delete[] y;
                delete[] dy;

                if(flagReturn == true){
                    return Py_BuildValue("i", timeoutFlag);

                }

                else{
                    return Py_BuildValue("d", N);
                }

            }

        flag = r8_rkf45(evolveB, 2*nF, y, dy, &N, Nstop, &relerr, abserr, flag, Cparams);

        // detect some error conditions
        if(flag == 50){

            delete[] y;
            delete[] dy;

            if(flagReturn == true){
                return Py_BuildValue("i", integrationFlag);
            }

            else{
                return Py_BuildValue("d", N);
            }

          }

        // compute value of epsilon
        vector<double> vecy(y, y+2*nF);
        vector<double> vecParams(Cparams, Cparams+nP);

        double eps = mm.Ep(vecy, vecParams);

        // If we have an unacceptable result, or if we are now past the end of inflation
        if(eps < 0 || eps > 1){

            delete[] y;
            delete[] dy;
            return Py_BuildValue("d", N); // Report integration lim. under epsilon condition

          }

        // '-2' means continue as normal in 'single-step' mode
        flag = -2;
      }

    // deallocate workspace
    delete[] y;
    delete[] dy;

    // No end to inflation
    if(flagReturn == true){
        return Py_BuildValue("i", eternalFlag);
    }
    else{
        return Py_BuildValue("d", N);
    }
}

// function to calculate background evolution
static PyObject* MT_backEvolve(PyObject* self,  PyObject *args)
{

    // Initialize python (array) objects; these should be parsed from the python side
    PyArrayObject *initialCs, *t, *backOut, *backOutT, *params, *tols;

    // Initialize corresponding ctype arrays
    double *CinitialCs, *tc, *Cparams, *tolsC ;

    // exit: if true -> exit when epsilon > 1
    bool exit;

    // setup flags
    int timeoutFlag     = -33;  // timeout flag
    int integrationFlag = -32;  // integration flag

    // Maximum integration time: Default -1 implies no max
    double tmax_back = -1;

    // If true, returns tuple with error flag, else N
    bool flagReturn = _SamplerMode();

    // Parse python arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!b|db",&PyArray_Type, &t, &PyArray_Type, &initialCs, &PyArray_Type, &params,
        &PyArray_Type, &tols, &exit, &tmax_back, &flagReturn)) {
        return NULL;
    }

    // Convert tols to c array
    tolsC = pyvector_to_Carray(tols);

    // If python array args are garbage, revert to default abs. and relative integrator tols.
    double abserr, relerr;
    if (2!=size_pyvector(tols)) {
        cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        abserr = pow(10,-8.); relerr = pow(10,-8.);
    }
    else {
        abserr =tolsC[0];relerr = tolsC[1];
    }

    // Convert initial conditions & parameters to c array
    CinitialCs = pyvector_to_Carray(initialCs);
    Cparams    = pyvector_to_Carray(params);

    // Cross check number of ICs and params against model
    model mm;
    int nF=mm.getnF();
    int nP = mm.getnP();

    if (2*nF!=size_pyvector(initialCs)) {
        cout<< "\n \n \n field space array not of correct length \n \n \n";
        Py_RETURN_NONE;
    }
    if (nP!=size_pyvector(params)) {
        cout<< "\n \n \n parameters array not of correct length \n \n \n";
        Py_RETURN_NONE;
    }

    // Convert time series (for background steps) to c array
    tc = pyvector_to_Carray(t);

    // Get Nstart for integrator
    double N=tc[0];

    // Build vecotr inputs for background integrator
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);
    back b(nF, vectIn);

    // Initialize integrator steps (vals and derivs.) as dynamic vars.
    int flag=-1;
    double *y;
    y = new double[2*nF];
    double *yp;
    yp= new double[2*nF];

    // Populate initial step with ICs
    for (int i=0;i<2*nF;i++) {
        y[i] = CinitialCs[i];
    }

    // Initialize time objects for measuring relative integrator times
    double deltat;
    time_t t_ini, t_fin;

    // Record start time
    time(&t_ini);

    // Get number of integration steps
    int nt = t->dimensions[0];

    // Create 2 dimensional array to populate with background calculation
    double backC[nt][1+2*nF];

    // Start counter for integration step
    int ii =0;

    // Evolve background
    evolveB(N, y, yp, Cparams);


    // If exit routine is false integrate target number of background steps as far as possible
    if (exit == 0){

//        cout << "Entering extended back routine" << endl;

        // While N < Nmax
        while (ii < nt){

            // Record current time
            time(&t_fin);

            // Measure time difference since start
            deltat = difftime(t_fin,t_ini);

            // If time exceeds max integration time
            if((deltat > tmax_back) && (tmax_back>0)){

//                    // Deallocate integrator work space and return timeout flag
//                    delete[] y;
//                    delete[] yp;
//
//                    // Change to return background!
//                    if(flagReturn == true){
//                        return Py_BuildValue("id", timeoutFlag, N);
//                    }
//
//                    else{
//                        return Py_BuildValue("d", N);
//                    }
                    break;

                }

            // Otherwise compute background step and get (flag for) integrator status
            flag = r8_rkf45(evolveB , 2*nF, y, yp, &N, tc[ii], &relerr, abserr, flag, Cparams);

            // If integrator error is raised
            if (flag== 50){
//                // Deallocate integrator work space and return integration failure flag for background
//                delete[] y;
//                delete[] yp;
//
//                if(flagReturn == true){
//                    return Py_BuildValue("id", integrationFlag, N);
//                }
//
//                else{
//                    return Py_BuildValue("d", N);
//                }
                break;

            }

            // Otherwise continue in single step mode
            flag=-2;

            // Populate background C array for current step (row) value
            backC[ii][0]=N;
            for (int i = 0; i < 2*nF; i++){
                backC[ii][1+i] = y[i];
            }

            // Increment counter
            ii += 1;

        }
    }

    // If exit is true, integrate up until \epsilon = 1 or N = Nend, whichever first
    if (exit == 1){

//        cout << "Entering epsilon back routine" << endl;

        // Define vector values for step values and model params (used for computing epsilon)
        vector<double> vecy;
        vector<double> Vparams;

        // Initialize \epsilon value
        double eps=0.0;

        // While \epsilon < 1 and N < Nmax
        while (eps<1 && ii < nt){

            // Record current time
            time(&t_fin);

            // Measure time difference since start of integration
            deltat = difftime(t_fin,t_ini);

            // If time difference is greater than maximum integration time
            if((deltat > tmax_back) && (tmax_back>0)){

                // Deallocate integrator workspace and return timeout flag
                delete[] y;
                delete[] yp;

                if(flagReturn == true){
                    return Py_BuildValue("i", timeoutFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }

            }

            // Compute integration step and get (flag for) integrator status
            flag = r8_rkf45(evolveB , 2*nF, y, yp, &N, tc[ii], &relerr, abserr, flag, Cparams );

            // Integration error
            if (flag== 50){

                // deallocate integrator workspace and return background integration failure flag
                delete[] y;
                delete[] yp;

                if(flagReturn == true){
                    return Py_BuildValue("i", integrationFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }
            }

            // integration continues in single step-mode
            flag=-2;

            // Populate background C array for current step (row) value
            backC[ii][0]=N;
            for (int i = 0; i < 2*nF; i++){
                backC[ii][1+i] = y[i];
            }

            // Compute \epsilon and increment counter
            vecy = vector<double>(y, y + 2*nF);
            Vparams = vector<double>(Cparams, Cparams +  nP);
            eps = mm.Ep(vecy,Vparams);
            ii = ii+1;

            }
    }


    // Now create python array object with background computation

    // Build appropriate dimensions for numpy array
    npy_intp dims[2];
    dims[1]=1+2*nF;
    dims[0]=ii; // Last iteration count

    // Create python array object
    backOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);

    // Create c++ pointer to python array
    double * backOutC;
    backOutC = (double *) PyArray_DATA(backOut);

    // Populate python memory values with background data
    for (int jj = 0; jj < ii; jj++){
        backOutC[jj*(2*nF+1)]=backC[jj][0];
        for (int kk = 0; kk < 2*nF; kk++){
             backOutC[jj*(2*nF+1)+kk+1]=backC[jj][1+kk];
        }
    }

    // Deallocate integrator workspace and return python array for background evolution
    delete[] y; delete[] yp;
    return PyArray_Return(backOut);
}


// function to calculate 2pt evolution
static PyObject* MT_sigEvolve(PyObject* self,  PyObject *args)
{

    PyArrayObject *initialCs, *t, *sigOut, *params, *tols;
    double *CinitialCs, *tc, k, *Cparams, *tolsC;
    bool full;

    // Setup flags
    int integrationFlag = -52;
    int timeoutFlag = -53;

    // Report tuple flag or N at fail
    bool flagReturn = _SamplerMode();

    // Set max integration time
    double tmax_2pf = -1;

    if (!PyArg_ParseTuple(args, "O!dO!O!O!bdb", &PyArray_Type, &t, &k, &PyArray_Type, &initialCs,&PyArray_Type, &params,
        &PyArray_Type, &tols, &full, &tmax_2pf, &flagReturn)) {
        return NULL;
    }

    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);

    tolsC = pyvector_to_Carray(tols);
    double rtol, atol;
    if (2!=size_pyvector(tols)){cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        atol = pow(10,-8.); rtol = pow(10,-8.);}
    else {
        atol =tolsC[0];rtol = tolsC[1];}


    model mm;
    potential pott;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "\n \n \n field space array not of correct length, not proceeding further \n \n \n";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);

    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length, not proceeding further \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

    // we use a scaling below that we rescale back at the end (so the final answer is as if the scaling was never there -- this helps standarise the rtol and atol needed for the same model run with differnet initial conditions
    double kn = 1.0;
    double kscale = k;
    double Nstart=tc[0] - log(kscale);


    sigma sig(nF, kn, Nstart, vectIn, Vparams) ; // instance of sigma object which fixs ics

    double* y; // set up array for ics
    y = new double[2*nF + 2*nF*2*nF];

    for(int i=0; i<2*nF;i++){y[i] = CinitialCs[i];} // fix values of input array
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF+ i+2*nF*j] = sig.getS(i,j);}}

    double* paramsIn; // array of parameters to pass to LHS of ODE routine

    paramsIn = new double[1+nP];
    for(int i=0; i<nP;i++) paramsIn[i]=Vparams[i];
    paramsIn[nP]=kn;

    // evolve a 2pt run **************************

    double N=Nstart;
    double* yp ; yp = new double [2*nF + 2*nF*2*nF];
    vector<double> Ni;
    double zz=0;

    int flag=-1;
    evolveSig(N, y, yp, paramsIn);
    vector<double> fieldIn(2*nF);
    fieldIn = vector<double>(y,y+2*nF);
    Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
    zz=0;
    for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
        zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}
    }

    int nt = t->dimensions[0];

    int size;
    if (full ==true){size = 1+2*nF + 1+ 2*nF*2*nF;}
    if (full ==false){size = 1 + 1;}

    npy_intp dims[2];
    dims[1]=size; dims[0]=nt;
    double * sigOutC;
    sigOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    sigOutC = (double *) PyArray_DATA(sigOut);

    // Initialize time objects for measuring relative integrator times
    double deltat;
    time_t t_ini, t_fin;

    // Record start time
    time(&t_ini);

    // Start integration
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]-log(kscale)){

            // Check whether maximum time has exceeded
            time(&t_fin);
            deltat = difftime(t_fin,t_ini);

            // If time difference is greater than maximum integration time
            if((deltat > tmax_2pf) && (tmax_2pf>0)){

                // Deallocate integrator workspace and return timeout flag
                delete [] y;
                delete [] yp;
                delete [] paramsIn;

                if(flagReturn == true){
                    return Py_BuildValue("i", timeoutFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }

            }

            // Compute step
            flag = r8_rkf45(evolveSig , 2*nF+2*nF*2*nF, y, yp, &N, tc[ii]-log(kscale), &rtol, atol, flag, paramsIn );

            if (flag== 50){

                delete [] y;
                delete [] yp;
                delete [] paramsIn;

                if(flagReturn == true){
                    return Py_BuildValue("i", integrationFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }

            }
            flag = -2;

        }

        fieldIn = vector<double>(y,y+2*nF);

        sigOutC[ii*size] = N+log(kscale);

        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        zz=0;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}}

        sigOutC[ii*size+1] = zz/kscale/kscale/kscale;


        if(full==true){
            for(int i=0;i<2*nF;i++)
            {
                sigOutC[ii*(size)+i+2]=y[i];
            }
            for(int i=2*nF;i<2*nF+ 2*nF*2*nF;i++)
            {
                sigOutC[ii*(size)+i+2]=y[i]/kscale/kscale/kscale;
            }
        }

    }

    delete [] y;
    delete [] yp;
    delete [] paramsIn;

    return PyArray_Return(sigOut);
}


// function to calculate 3pt evolution
static PyObject* MT_alphaEvolve(PyObject* self,  PyObject *args)
{

    PyArrayObject *initialCs, *t, *alpOut,*params, *tols;

    double k1, k2, k3, Nstart, *CinitialCs, *tc,*Cparams, *tolsC;

    bool full;

    // Setup flags
    int integrationFlag = -62;
    int timeoutFlag = -63;
    bool flagReturn = _SamplerMode();

    // Maximum integration time, if < 0 then no max
    double tmax_3pf = -1;

    if (!PyArg_ParseTuple(args, "O!dddO!O!O!bdb", &PyArray_Type, &t, &k1, &k2, &k3, &PyArray_Type, &initialCs,
        &PyArray_Type,&params,&PyArray_Type,&tols, &full, &tmax_3pf, &flagReturn)) {
        return NULL;
    }

    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    int nt = t->dimensions[0];

    tolsC = pyvector_to_Carray(tols);
    double rtol, atol;
    if (2!=size_pyvector(tols)){cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        atol = pow(10,-8.); rtol = pow(10,-8.);}
    else {
        atol =tolsC[0];rtol = tolsC[1];}

 //   if (full != 0 && full !=1 ){ full=1;cout << "\n \n \n Number out of range, defaulted to full outout mode \n \n \n";}
    model mm;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "\n \n \n field space array not of correct length, not proceeding further \n \n \n";    Py_RETURN_NONE;}

    int nP = mm.getnP();if (nP!=size_pyvector(params)){cout<< "\n \n \n  parameters array not of correct length, not proceeding further \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);

// we use a scaling below that we rescale back at the end (so the final answer is as if the scaling was never there -- this helps standarises the rtol and atol needed for the same model run with differnet initial conditions

    double kscale = (k1+k2+k3)/3.;
    double k1n = k1/kscale; double k2n = k2/kscale; double k3n = k3/kscale;
    Nstart=tc[0] -log(kscale);
    double N=Nstart; // reset N

    // do not alter the comment at the end of the next line -- used by preprocessor
    // ****************************************************************************

    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs+2*nF);

    sigma sig1(nF, k1n, Nstart, vectIn,Vparams)  ; // 3 instances of sigmmas
    sigma sig2(nF, k2n, Nstart, vectIn,Vparams)  ;
    sigma sig3(nF, k3n, Nstart, vectIn,Vparams)  ;
    sigmaI sig1I(nF, k1n, Nstart, vectIn,Vparams)  ; // 3 instances of sigma imaginary
    sigmaI sig2I(nF, k2n, Nstart, vectIn,Vparams)  ;
    sigmaI sig3I(nF, k3n, Nstart, vectIn,Vparams)  ;
    alpha alp(nF, k1n, k2n, k3n, Nstart, vectIn, Vparams); // instance of alpha

    double* y; // array for initial conditions
    y = new double[2*nF + 6*2*nF*2*nF + 2*nF*2*nF*2*nF];

    for(int i=0; i<2*nF;i++){y[i] = CinitialCs[i];}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF+ i+2*nF*j] = sig1.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 1*(2*nF*2*nF)+ i+2*nF*j] = sig2.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 2*(2*nF*2*nF)+ i+2*nF*j] = sig3.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 3*(2*nF*2*nF)+ i+2*nF*j] = sig1I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 4*(2*nF*2*nF)+ i+2*nF*j] = sig2I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 5*(2*nF*2*nF)+ i+2*nF*j] = sig3I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){for(int k=0; k<2*nF;k++){y[2*nF + 6*(2*nF*2*nF)+ i+2*nF*j + 2*nF*2*nF*k] = alp.getA(i,j,k);}}}


    double* paramsIn2; // array for parameters of RHS of ODE routine
    paramsIn2 = new double[3+nP];
    for(int i=0; i<nP;i++) paramsIn2[i]=Vparams[i];
    paramsIn2[nP]=k1n;
    paramsIn2[nP+1]=k2n;
    paramsIn2[nP+2]=k3n;


    double ZZZ=0., ZZ1=0., ZZ2=0., ZZ3=0.; //  for zeta zeta calcs
    vector<double> Ni, Nii1, Nii2, Nii3 ; // for N transforms to get to zeta
    double *yp; yp=new double[2*nF +6*2*nF*2*nF+  2*nF*2*nF*2*nF];

    npy_intp dims[2];
    int size;
    if (full==false){size =   5;}
    if (full==true){size =  5+  2*nF + 6*2*nF*2*nF+2*nF*2*nF*2*nF;}
    dims[1]=size; dims[0]=nt;
    double * alpOutC;
    alpOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    alpOutC = (double *) PyArray_DATA(alpOut);

    evolveAlp(N, y, yp, paramsIn2);
    int flag=-1;

    // Initialize time objects for measuring relative integrator times
    double deltat;
    time_t t_ini, t_fin;

    // Record start time
    time(&t_ini);

    // run alpha *******************************************
    vector<double> fieldIn(2*nF);

    for(int ii=0; ii<nt; ii++ ){
        while (N<(tc[ii]-log(kscale))){

            // Check whether maximum time has exceeded
            time(&t_fin);
            deltat = difftime(t_fin,t_ini);

            // If time difference is greater than maximum integration time
            if((deltat > tmax_3pf) && (tmax_3pf>0)){

                // Deallocate integrator workspace and return timeout flag
                delete [] y;
                delete [] yp;
                delete [] paramsIn2;

                if(flagReturn == true){
                    return Py_BuildValue("i", timeoutFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }

            }

            // Compute integration step
            flag = r8_rkf45(evolveAlp, 2*nF + 6*(2*nF*2*nF) + 2*nF*2*nF*2*nF, y, yp, &N,
                            tc[ii]-log(kscale), &rtol, atol, flag, paramsIn2);

            // If integration failure
            if (flag== 50){

                // Deallocate integrator workspace and return integrator flag
                delete [] y;
                delete [] yp;
                delete [] paramsIn2;

                if(flagReturn == true){
                    return Py_BuildValue("i", integrationFlag);
                }

                else{
                    return Py_BuildValue("d", N);
                }

            }

            flag=-2;
        }

        fieldIn = vector<double>(y,y+2*nF);
        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        Nii1=mm.N2(fieldIn,Vparams,k1n,k2n,k3n,N); // claculate N,ij array for first arrangement of ks
        Nii2=mm.N2(fieldIn,Vparams,k2n,k1n,k3n,N); // for second
        Nii3=mm.N2(fieldIn,Vparams,k3n,k1n,k2n,N); // etc

        ZZ1=0.;
        ZZ2=0.;
        ZZ3=0.;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            ZZ1=ZZ1+Ni[i]*Ni[j]*(y[2*nF + i + j*2*nF] );
            ZZ2=ZZ2+Ni[i]*Ni[j]*y[2*nF + (2*nF*2*nF) + i + j*2*nF];
            ZZ3=ZZ3+Ni[i]*Ni[j]*y[2*nF + 2*(2*nF*2*nF) + i + j*2*nF];
        }}


        ZZZ=0.;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF;j++){for(int k=0; k<2*nF;k++){
            ZZZ=ZZZ + Ni[i]*Ni[j]*Ni[k]*y[2*nF + 6*(2*nF*2*nF) + i + j*2*nF+ k*2*nF*2*nF];
            for(int l=0; l<2*nF;l++){ZZZ=ZZZ+(Nii1[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 1*(2*nF*2*nF) + i+k*2*nF]*y[2*nF+2*(2*nF*2*nF)+j+l*2*nF]
                                              +Nii2[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 0*(2*nF*2*nF) + i+k*2*nF]*y[2*nF + 2*(2*nF*2*nF) + j+l*2*nF]
                                              +Nii3[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 0*(2*nF*2*nF) + i+k*2*nF]*y[2*nF + 1*(2*nF*2*nF) + j+l*2*nF]);
            }}}}

        alpOutC[ii*size] =  N+log(kscale);
        alpOutC[ii*size+1] = ZZ1/kscale/kscale/kscale;
        alpOutC[ii*size+2] = ZZ2/kscale/kscale/kscale;
        alpOutC[ii*size+3] = ZZ3/kscale/kscale/kscale;
        alpOutC[ii*size+4] = ZZZ/kscale/kscale/kscale/kscale/kscale/kscale;

        if(full==true){
            for(int i=0;i<2*nF ;i++){
                alpOutC[ii*size+5+i] =  y[i] ;   }

            for(int i=2*nF;i<2*nF + 6*(2*nF*2*nF);i++){
                alpOutC[ii*size+5+i] =  y[i]/kscale/kscale/kscale ;   }

            for(int i=2*nF + 6*(2*nF*2*nF);i<2*nF + 6*(2*nF*2*nF)+ 2*nF*2*nF*2*nF;i++){
                alpOutC[ii*size+5+i] =  y[i]/kscale/kscale/kscale/kscale/kscale/kscale ;   }
        }

    }

    delete [] y;
    delete [] yp;
    delete [] paramsIn2;

    return PyArray_Return(alpOut);
}

PyDoc_STRVAR(doc_H, "H(fields_dotfields: np.ndarray, params: np.ndarray) -> float\n\nComputes the Hubble rate of expansion.\n\nParameters\n----------\nfields_dotfields : np.ndarray\nnumpy array of shape (1, 2*num_fields) of dtype float. the field values followed by field velocities.\nparams : np.ndarray,\nnumpy array of shape (1, num_params) of dtype float containing the model parameters.\n\nReturns\n-------\nHubbleRate : float\nHubble rate at input parameters.");
PyDoc_STRVAR(doc_nF, "nF() -> int\n\nGet number of fields in model.\n\nReturns\n-------\nnF: int\nNumber of fields associated with model.");
PyDoc_STRVAR(doc_nP, "nP() -> int\n\nGet number of parameters in model.\n\nReturns\n-------\nnP: int\nNumber of parameters associated with model.");
PyDoc_STRVAR(doc_Epsilon, "Eps(fields_dotfields: np.ndarray, params: np.ndarray) -> float\n\nCompute epsilon slow-roll parameter in terms of the Hubble rate.\n\nParameters\n----------\nfields_dotfields : numpy array of shape (1, 2*num_fields) of dtype float\nthe field values followed by field velocities.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n-------\nEpsilon : float\nEpsilon slow-roll parameter.");
PyDoc_STRVAR(doc_Eta, "Eta(fields_dotfields: np.ndarray, params: np.ndarray) -> float\n\nCompute eta slow-roll parameter in terms of the Hubble rate.\n\nParameters\n----------\nfields_dotfields : numpy array of shape (1, 2*num_fields) of dtype float\nthe field values followed by field velocities.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n-------\nEta : float\nEta slow-roll parameter.");
PyDoc_STRVAR(doc_V, "V(fields: np.ndarray, params: np.ndarray) -> float\n\nCompute the potential at given field values\n\nParameters\n----------\nfields : np.ndarray of shape(1, num_fields) of dtype float\nthe field values.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n------\npotential : float\nValue of potential with input field and parameter combination.");
PyDoc_STRVAR(doc_dV, "dV(fields: np.ndarray, params: np.ndarray) -> np.ndarray\n\nCompute 1st (covariant) derivative of potential at given field values\n\nParameters\n----------\nfields : np.ndarray of shape(1, num_fields) of dtype float\nthe field values.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n------\n1st potential derivative : numpy array of shape (1, num_fields) of type float\nValue of potential derivatives with input field and parameter combination.");
PyDoc_STRVAR(doc_ddV, "ddV(fields: np.ndarray, params: np.ndarray) -> np.ndarray\n\nCompute 2nd (covariant) derivative of potential at given field values\n\nParameters\n----------\nfields : np.ndarray of shape(1, num_fields) of dtype float\nthe field values.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n------\n2nd potential derivative : numpy array of shape (num_fields, num_fields) of type float\nValue of potential derivatives with input field and parameter combination.");
PyDoc_STRVAR(doc_findEndOfInflation, "findEndOfInflation(fields_dotfields: np.ndarray, params: np.ndarray, tols: np.ndarray, Ninit=0.0, Nmax=3000.0, tmax=-1, flag_return=False) -> float\n\nFinds end of inflation; i.e. the first efold when epsilon > 1.\n\nParameters\n----------\nfields_dotfields : numpy array of shape (1, 2*num_fields) of dtype float\nthe initial conditions for evolution.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\n\nReturns\n-------\nNend : float\nEnd of inflation.");
PyDoc_STRVAR(doc_massMatrix, "massMatrix(fields_dotfields: np.ndarray, params: np.ndarray, hessian=False, covariant=False) -> np.ndarray\n\nCompute mass (squared) matrix from background\n\nParameters\n----------\nfields_dotfields : numpy array of shape (1, 2*num_fields) of dtype float\nthe field values followed by field velocities.\nparams : numpy array of shape (1, num_params) of dtype float\nthe model parameters.\nhessian : if True, returns Hessian approximation to matrix, otherwise retains curvature and kinetic terms.\ncovariant : if True, evaluates index down-down expression; else up-down.\n\nReturns\n-------\nMassMatrix : numpy array of shape (num_fields, num_fields) of dtype float\nThe mass matrix evaluated at input configuration.");
PyDoc_STRVAR(doc_backEvolve, "backEvolve(Narray: np.ndarray, ics: np.ndarray, params: np.ndarray, tols: np.ndarray, full: bool, tmax=-1, flag_return=False) -> np.ndarray\n\nCompute background evolution for a range of efoldings\n\nParameters\n----------\nNarray : np.ndarray,\n1D numpy array corresponding to target efolds for integration.\nics : numpy array of dtype float of size (1, 2*num_fields),\ncorresponding to the initial conditions of evolution.\nparams : numpy array of dtype float of size (1, num_params),\ncorresponding to the model parameters.\ntols : numpy array of dtype float of size 2,\ncorresponding to the abs. and rel. integration tolerances.\nexit : boolean,\nif True, integrate until slow roll violation, else attempt maximal integration.\ntmax : int,\nmaximum time (wall-clock) to attempt integration before bad exit. Defaults to -1, corresponding to unconstrained.\nflag_return : bool,\nif True, returns flag containing error code for bad integration procedure. Otherwise returns Null.\n\nReturns\n-------\nback : np.ndarray\nA two dimensional numpy array: contains the fields, and field velocities at the times contained in Narray. The format is that the array has 1 + 2nF columns, with the zeroth column the times (Narray), the next columns are the field values and field velocity values at those times.");
PyDoc_STRVAR(doc_sigEvolve, "sigEvolve(Narray: np.ndarray, k: float, ics: np.ndarray, tols: np.ndarray, full: bool, tmax=-1, flag_return=False) -> np.ndarray\n\nEvolve 2pt correlation functions\n\nParameters\n----------\nNarray : np.ndarray\n1D numpy array of type float, containing the efold times for wich outputs should be evaluated and returned. First and last elements should correspond to initial and final efolding of evolution.\nk : float,\nFourier mode.\nics : numpy array of dtype float of size (1, 2*num_fields),\ncorresponding to the initial conditions of evolution.\ntols : numpy array of dtype float of size 2,\ncorresponding to the abs. and rel. integration tolerances.\nfull : boolean,\nThe function returns a two dimensional numpy array the zeroth column of which contains the times (Narray). If full=False the next column contains the power spectrum of ζ at these times, and this is the final column (there are therefore only 2 columns in total). If full = True then there is in addition 2 nF + 2 nF * 2 nF columns. See docs @ https://arxiv.org/pdf/1609.00381.pdf.\nflag_return : bool,\nif True, returns flag containing error code for bad integration procedure. Otherwise returns Null.\n\nReturns\n-------\nsigma : np.ndarray\nTwo dimensional array representing evolution of 2-point correlators. (see description of 'full' parameter.)");
PyDoc_STRVAR(doc_alphaEvolve, "alphaEvolve(Narray: np.ndarray, k1: float, k2: float, k3: float, ics: np.ndarray, tols: np.ndarray, full: bool, tmax=-1, flag_return=False) -> np.ndarray\n\nEvolve 3pt correlation functions\n\nParameters\n----------\nNarray : np.ndarray\n1D numpy array of type float, containing the efold times for wich outputs should be evaluated and returned. First and last elements should correspond to initial and final efolding of evolution.\nk1 : float,\nfirst Fourier mode.\nk2 : float,\nsecondFourier mode.\nk3 : float,\nthird Fourier mode.\nics : numpy array of dtype float of size (1, 2*num_fields),\ncorresponding to the initial conditions of evolution.\ntols : numpy array of dtype float of size 2,\ncorresponding to the abs. and rel. integration tolerances.\nfull : boolean,\nThe function returns a two dimensional numpy array, the first column of which contains the times (Narray). If full=False the next four columns contains the power spectrum of zeta for each of the three k values input, and the value of the bispectrum of zeta for a triangle with sides of length of these k values. A total of [1+ 4] columns. If full=True, there are an additional 2*nF + 6 *2 nF* 2 nF + 2 nF * 2 nF * 2 nF columns. The first 2nF of these columns contain the fields, and field velocities at the time steps requested (the background cosmology). See docs in appendix of @ https://arxiv.org/pdf/1609.00381.pdf.\ntmax : int,\nmaximum time (wall-clock) to attempt integration before bad exit. Defaults to -1, corresponding to unconstrained.\nflag_return : bool,\nif True, returns flag containing error code for bad integration procedure. Otherwise returns Null.\n\nReturns\n-------\nalpha : np.ndarray\nTwo dimensional array representing evolution of 3-point correlators. (see description of 'full' parameter.)");

static char PyTrans_docs[] =
"This is a custom PyTransport package for solving the moment transport equations of inflationary cosmology\n";

// **************************************************************************************
static PyMethodDef dquad_euclidean_funcs[] = {{"H", (PyCFunction)MT_H,   METH_VARARGS, doc_H},{"nF", (PyCFunction)MT_fieldNumber,   METH_NOARGS, doc_nF},{"nP", (PyCFunction)MT_paramNumber,   METH_NOARGS, doc_nP},{"Epsilon", (PyCFunction)MT_Ep,   METH_VARARGS, doc_Epsilon},{"Eta", (PyCFunction)MT_Eta,   METH_VARARGS, doc_Eta},{"V", (PyCFunction)MT_V,   METH_VARARGS, doc_V},{"dV", (PyCFunction)MT_dV,   METH_VARARGS, doc_dV},{"ddV", (PyCFunction)MT_ddV,   METH_VARARGS, doc_ddV},{"massMatrix", (PyCFunction)MT_massMatrix,   METH_VARARGS, doc_massMatrix},{"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,   METH_VARARGS, doc_findEndOfInflation},{"backEvolve", (PyCFunction)MT_backEvolve,   METH_VARARGS, doc_backEvolve},{"sigEvolve", (PyCFunction)MT_sigEvolve,   METH_VARARGS, doc_sigEvolve},{"alphaEvolve", (PyCFunction)MT_alphaEvolve,   METH_VARARGS, doc_alphaEvolve},   {NULL}};//FuncDef
// do not alter the comment at the end of preceeding line -- it is used by preprocessor

#ifdef __cplusplus
extern "C" {
#endif

// **************************************************************************************
//modDef
// do not alter the comment at the end of preceeding line -- it is used by preprocessor

// **************************************************************************************
//initFunc
// do not alter the comment at the end of preceeding line -- it is used by preprocessor

#ifdef __cplusplus
}
#endif