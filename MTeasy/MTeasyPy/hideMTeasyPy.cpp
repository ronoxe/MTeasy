// C++ file which defines the functions make available to Python through the MTeasy module.
#include <Python.h>
#include <iostream>
#include "numpy/arrayobject.h"
//don't adjust the labels at the end of the 4 lines below (they are used to fix directory structure)
#include"/home/jwr/Code/MTeasy/MTeasyC/evolve.h"//evolve
#include"/home/jwr/Code/MTeasy/MTeasyC/moments.h"//moments
#include"/home/jwr/Code/MTeasy/MTeasyC/model.h"//model
#include"/home/jwr/Code/MTeasy/MTeasyC/stepper/rkf45.hpp"//stepper
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
// Package recompiled: Thu May 19 12:00:45 2016


// Changes python array into C array (or rather points to pyarray data)
//    Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carray(PyArrayObject *arrayin)
{
    //  int i,n;
    //	n=arrayin->dimensions[0]; n is length of python array
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
    int nF = pp.getnF(); if (nF!=size_pyvector(fieldsIn)){cout<< "field space array not of correct length";    Py_RETURN_NONE;}
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length";  Py_RETURN_NONE;}
        
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
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "field space array not of correct length";    Py_RETURN_NONE;}
    vector<double> vectIn;
    int dims[1];
    dims[0]=nF;
    
    dVI = (PyArrayObject*) PyArray_FromDims(1,dims,NPY_DOUBLE);
    dVC = (double*) dVI->data;
    
    vectIn = vector<double>(Cfields, Cfields +  nF);
    
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length";  Py_RETURN_NONE;}
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
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "field space array not of correct length";    Py_RETURN_NONE;}
    
    vector<double> vectIn;
    int dims[2];
    dims[0]=nF; dims[1]=nF;
    
    ddVI = (PyArrayObject*) PyArray_FromDims(2,dims,NPY_DOUBLE);
    ddVC = (double*) ddVI->data;
    
    vectIn = vector<double>(Cfields, Cfields +  nF);
    
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length";  Py_RETURN_NONE;}
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
    int nF = mm.getnF(); if (2*nF!=size_pyvector(fields_dfieldsIn)){cout<< "field space array not of correct length";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);
    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    return Py_BuildValue("d", mm.H(vectIn, Vparams));
}

// function to return number of fields (useful for cross checks)
static PyObject* MT_fieldNumber(PyObject* self,  PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;}
    model mm;
    return Py_BuildValue("i",mm.getnF());
}

// function to return number of parameters (useful for cross checks)
static PyObject* MT_paramNumber(PyObject* self,  PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;}
    model mm;
    return Py_BuildValue("i",mm.getnP());
}

// function to calculate background evolution
static PyObject* MT_backEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject *initialCs, *t, *backOut, *params;
    double *CinitialCs, *tc, *Cparams ;
    //char *file;
    if (!PyArg_ParseTuple(args, "O!O!O!",&PyArray_Type, &t, &PyArray_Type, &initialCs,&PyArray_Type,&params)) {
        return NULL;}
    //ofstream outfile(file); outfile.precision(20);
    //outfile.close();
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    model mm;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "field space array not of correct length";    Py_RETURN_NONE;}

    double N=tc[0];
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);
    back b(nF, vectIn);
    
       int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    //   vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    double abserr = pow(10,-10.); double relerr = pow(10,-14.);
    
    int flag=-1;
    double *y; y = new double[2*nF];
    double *yp; yp= new double[2*nF];
    
    //   ofstream outfile("../../runData/back.dat"); outfile.precision(20);
    
    for (int i=0;i<2*nF;i++){y[i] = CinitialCs[i];}
    
    int nt = t->dimensions[0];
    
    int dims[2];
    dims[1]=1+2*nF; dims[0]=nt;
    double * backOutC;
    backOut = (PyArrayObject*) PyArray_FromDims(2,dims,NPY_DOUBLE);
    backOutC = (double *) PyArray_DATA(backOut);
    evolveB(N, y, yp, Cparams);
    // run background *********************
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]){
            flag = r8_rkf45(evolveB , 2*nF, y, yp, &N, tc[ii], &relerr, abserr, flag, Cparams );
            flag=-2;
        }
        backOutC[ii*(2*nF+1)]=N;
        //        outfile<<N<<"\t";
        for(int i=0;i< 2*nF;i++){//*outfile << y[i]  <<"\t";
            backOutC[ii*(2*nF+1)+i+1]=y[i];} // outputs to file at each step
        //        outfile <<endl;
    }
    //   outfile.close();
    delete[] y; delete[] yp;
    return PyArray_Return(backOut);
}


static PyObject* MT_sigEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject  *initialCs, *fields_dfieldsIn, *t, *sigOut, *params;
    double *Cfields_dfields, *CinitialCs, *tc, k, *Cparams;
    int full;
    if (!PyArg_ParseTuple(args, "O!dO!O!i", &PyArray_Type, &t, &k, &PyArray_Type, &initialCs, &PyArray_Type, &params,&full)) {
        return NULL;}
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    
    if (full != 0 && full !=1 ){ full=1; cout << "Number out of range, defaulted to full outout mode\n";}
    
	Cfields_dfields = pyvector_to_Carray(fields_dfieldsIn);
    model mm;
    potential pott;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "field space array not of correct length, not proceeding further\n";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);
    
    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length, not proceeding further\n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    //    ofstream outfileZZ("../../runData/zz.dat"); outfileZZ.precision(35);
    //    ofstream outfile2("../../runData/sigma.dat"); outfile2.precision(35);
        
    double Nstart = tc[0];
    sigma sig(nF, k, Nstart, vectIn, Vparams) ; // instance of sigma object which fixs ics
    
    double* y; // set up array for ics
    y = new double[2*nF + 2*nF*2*nF];
    
    for(int i=0; i<2*nF;i++){y[i] = CinitialCs[i];} // fix values of input array
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF+ i+2*nF*j] = sig.getS(i,j);}}
    
    double* paramsIn; // array of parameters to pass to LHS of ODE routine
    
    paramsIn = new double[1+nP];
    for(int i=0; i<nP;i++) paramsIn[i]=Vparams[i];
    paramsIn[nP]=k;
    
    // evolve a 2pt run **************************
	vector<double> f;
	f = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);
    double N=tc[0];
    double rtol=pow(10,-8.); double atol=pow(10,-14.) ; 
    double* yp ; yp = new double [2*nF + 2*nF*2*nF];
    vector<double> Ni;
    double zz=0;
    
    int flag=-1;
    evolveSig(N, Cfields_dfields, y, yp, paramsIn);
    vector<double> fieldIn(2*nF);
    fieldIn = vector<double>(y,y+2*nF);
    Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
    zz=0;
    for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
        zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}
    }
    
    int nt = t->dimensions[0];
    
    int size;
    if (full ==1){size = 1+2*nF + 1+ 2*nF*2*nF;}
    if (full ==0){size = 1 + 1;}
    
    int dims[2];
    dims[1]=size; dims[0]=nt;
    double * sigOutC;
    sigOut = (PyArrayObject*) PyArray_FromDims(2,dims,NPY_DOUBLE);
    sigOutC = (double *) PyArray_DATA(sigOut);
    
    
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]){
            flag = r8_rkf45(evolveSig , 2*nF+2*nF*2*nF, y, yp, &N, tc[ii], &rtol, atol, flag, paramsIn );
            flag = -2;
        }
        fieldIn = vector<double>(y,y+2*nF);
        
        sigOutC[ii*size] = N;
        
        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        zz=0;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}}
        //        outfileZZ << N <<"\t" << zz;
        //        outfileZZ << endl;
        
        sigOutC[ii*size+1] = zz;
        
        
        if(full==1){
            //        outfile2<<N<<"\t";
            for(int i=0;i<2*nF+ 2*nF*2*nF;i++)
            {
                sigOutC[ii*(size)+i+2]=y[i];
                //            outfile2<< y[i] <<"\t";
            }
            //		outfile2<<endl;
        }
        
        
        
    }
    
    //    outfileZZ.close();
    //    outfile2.close();
    
    delete [] y; delete [] yp;
    delete [] paramsIn;
    return PyArray_Return(sigOut);
}

//more functions
static PyObject* MT_alphaEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject *initialCs, *t, *alpOut,*params;
    double k1, k2, k3, Nstart, *CinitialCs, *tc,*Cparams;
    int full;
    if (!PyArg_ParseTuple(args, "O!dddO!O!i", &PyArray_Type, &t, &k1,&k2,&k3, &PyArray_Type, &initialCs,&PyArray_Type,&params, &full)) {
        return NULL; }
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    int nt = t->dimensions[0];
    if (full != 0 && full !=1 ){ full=1;cout << "Number out of range, defaulted to full outout mode\n";}
    model mm;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "field space array not of correct length, not proceeding further\n";    Py_RETURN_NONE;}
    
    int nP = mm.getnP();if (nP!=size_pyvector(params)){cout<< "parameters array not of correct length, not proceeding further\n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    
    //ofstream outfile3("../../runData/alpha.dat"); // file for raw 3 pt
    //outfile3.precision(20);
    //ofstream outfile4("../../runData/zzz.dat"); // file for zeta correlations
    //outfile3.precision(20);
    //
    Nstart=tc[0];
    double N=Nstart; // reset
    double rtol=1e-12, atol=1e-18;//tols
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs+2*nF);
    
    sigma sig1(nF, k1, Nstart, vectIn,Vparams)  ; // 3 instances of sigmmas
    sigma sig2(nF, k2, Nstart, vectIn,Vparams)  ;
    sigma sig3(nF, k3, Nstart, vectIn,Vparams)  ;
    sigmaI sig1I(nF, k1, Nstart, vectIn,Vparams)  ; // 3 instances of sigma imaginary
    sigmaI sig2I(nF, k2, Nstart, vectIn,Vparams)  ;
    sigmaI sig3I(nF, k3, Nstart, vectIn,Vparams)  ;
    alpha alp(nF, k1, k2, k3, Nstart, vectIn,Vparams); // instance of alpha
    
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
    paramsIn2[nP]=k1;
    paramsIn2[nP+1]=k2;
    paramsIn2[nP+2]=k3;
    
    
    double ZZZ=0., ZZ1=0., ZZ2=0., ZZ3=0.; //  for zeta zeta calcs
    vector<double> Ni, Nii1, Nii2, Nii3 ; // for N transforms to get to zeta
    double *yp; yp=new double[2*nF +6*2*nF*2*nF+  2*nF*2*nF*2*nF];
    
    int dims[2];
    int size;
    if (full==0){size =   5;}
    if (full==1){size =  5+  2*nF + 6*2*nF*2*nF+2*nF*2*nF*2*nF;}
    dims[1]=size; dims[0]=nt;
    double * alpOutC;
    alpOut = (PyArrayObject*) PyArray_FromDims(2,dims,NPY_DOUBLE);
    alpOutC = (double *) PyArray_DATA(alpOut);
    
    
    evolveAlp(N, y, yp, paramsIn2);
    int flag=-1;
    
    // run alpha *******************************************
    vector<double> fieldIn(2*nF);
    
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]){
            flag = r8_rkf45(evolveAlp, 2*nF + 6*(2*nF*2*nF) + 2*nF*2*nF*2*nF, y, yp, &N, tc[ii], &rtol, atol, flag, paramsIn2);
            flag=-2;
        }
        
        fieldIn = vector<double>(y,y+2*nF);
        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        Nii1=mm.N2(fieldIn,Vparams,k1,k2,k3,N); // claculate N,ij array for first arrangement of ks
        Nii2=mm.N2(fieldIn,Vparams,k2,k1,k3,N); // for second
        Nii3=mm.N2(fieldIn,Vparams,k3,k1,k2,N); // etc
        
        ZZ1=0.;
        ZZ2=0.;
        ZZ3=0.;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            ZZ1=ZZ1+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];
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
        
        alpOutC[ii*size] =  N;
        alpOutC[ii*size+1] = ZZ1;
        alpOutC[ii*size+2] = ZZ2;
        alpOutC[ii*size+3] = ZZ3;
        alpOutC[ii*size+4] = ZZZ;
        //outfile4<<N<<"\t";
        //outfile4<<ZZ1<<"\t"<<ZZ2<<"\t"<<ZZ3<<"\t"<<ZZZ;
        //outfile4<<endl;
        
        if(full==1){
            //outfile3<<N<<"\t";
            for(int i=0;i<2*nF + 6*(2*nF*2*nF)+ 2*nF*2*nF*2*nF;i++){
                //outfile3<< y[i] <<"\t";
                alpOutC[ii*size+5+i] =  y[i] ;   }
            //outfile3 <<endl;
        }
        
        //outfile4<<N<<"\t";
        //outfile4<<ZZ1<<"\t"<<ZZ2<<"\t"<<ZZ3<<"\t"<<ZZZ;
        //outfile4<<endl;
    }
    
    
    //outfile3.close();
    //outfile4.close();
    
    delete [] y;  delete [] paramsIn2; delete [] yp;
    
    return PyArray_Return(alpOut);
}


static char MTeasyPy_docs[] =
"This MTeasy, a package for solving the moment transport equations of inflationary cosmology\n";

static PyMethodDef MTeasyPyAQ_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, MTeasyPy_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, MTeasyPy_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, MTeasyPy_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, MTeasyPy_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, MTeasyPy_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, MTeasyPy_docs} , {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, MTeasyPy_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, MTeasyPy_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, MTeasyPy_docs},    {NULL}};//FuncDef


#ifdef __cplusplus
extern "C" {
#endif
    
void initMTeasyPyAQ(void)    {        Py_InitModule3("MTeasyPyAQ", MTeasyPyAQ_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc

#ifdef __cplusplus
}
#endif