//model class file contains the defining features of the model -- the u1, u2 flow tesors, and the A B C tesors and u3 flow tensor as well as the guage transform N tensors  
#ifndef MODEL_H  // Prevents the class being re-defined
#define MODEL_H 

#include "fieldmetric.h"
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
    fieldmetric fmet;

public:
	// constructor
	model()
	{
	   potential pot;
	   fieldmetric fmet;
        nP=pot.getnP();
        nF=pot.getnF();
    }
    
    // function returns Hubble rate
	double H(vector<double> f, vector<double> p   )
	{
		double Hi2;
		double Vi;
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		
		Vi=pot.V(f,p);
		Hi2=0.;
        
		for(int i=0; i<nF; i++)
		{	for(int j=0; j<nF; j++)
			{
				Hi2=Hi2+1./3.*(FMi[(2*nF)*(i)+(j)]*f[nF+j]*f[nF+i]/2.);
				//if(i==j){Hi2=Hi2+1./3.*f[nF+i]*f[nF+j]/2.;}
			}
		}
		Hi2=Hi2 + 1./3.*Vi;
        return sqrt(Hi2);
	}

    // function returns H dot
    double Hdot(vector<double> f)
	{
        double sum=0.;
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		
		for(int i=0; i<nF; i++)
		{for(int j=0; j<nF; j++)
			{
				sum= sum - 1./2.*(FMi[(2*nF)*(i)+(j)]*f[nF+i]*f[nF+j]);
				//if(i==j){sum= sum - 1./2.*(f[nF+i]*f[nF+j]);}
			}
		}
		return sum;
	}
    
    // function returns epsilon
	double Ep(vector<double> f,vector<double> p)
	{
		double Hi = H(f,p);
		double mdotH=0.;
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		
		for(int i=0; i<nF; i++)
		{for(int j=0; j<nF; j++)
			{
			mdotH= mdotH + 1./2.*(FMi[(2*nF)*(i)+(j)]*f[nF+i]*f[nF+j]);
			//if(i==j){mdotH= mdotH + 1./2.*(f[nF+i]*f[nF+j]);}
			}
		}
  		return mdotH/(Hi*Hi);
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
	// calculates u1
	vector<double> u(vector<double> f,vector<double> p)
	{
		vector<double> u1out(2*nF);
		vector<double> dVi;
		double Hi;
		Hi=H(f,p);
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		
		for(int i=0; i<nF; i++)
		{
			u1out[i]  = f[nF+i]/Hi;
		}
		
		dVi=pot.dV(f,p);

		for(int i=0; i<nF; i++)			
		{	
			double sum =0.0;
			for(int j=0; j<nF; j++)
			{
				sum = sum-dVi[j]/Hi;
				//if(i==j){u1out[nF+i]  = -3.*Hi*f[nF+i]/Hi-dVi[j]/Hi;}
			}
			u1out[nF+i]  = -3.*Hi*f[nF+i]/Hi +sum;
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
    
		vector<double> dVVi;
		dVVi = pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		double sum1 = 0.0;
		
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		vector<double> RMi;
		RMi = fmet.Riemn(f);


		for(int i = 0; i<nF; i++){for(int j = 0; j<nF; j++){
			for(int l=0; l<nF; l++)
			{for(int m=0; m<nF; m++)
			{for(int n=0; n<nF; n++)
			{for(int p=0; p<nF; p++)
			{
				sum1=sum1+RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(j+nF)+(m+nF)]*Fmi[(2*nF)*(m)+(n)]*f[nF+n]*Fmi[(2*nF)*(l)+(p)]*f[nF+p];
			}
			}
			}
			}
            u2out[i+ j*2*nF]=0.;
            u2out[i+(j+nF)*2*nF]=+ Fmi[(2*nF)*(i+nF)+(j+nF)]*1./Hi /a;
            u2out[i+nF+(j)*2*nF]=(-dVVi[i + nF*j] + (-3.+ep)*f[nF+i]*f[nF+j] + 1./Hi*(-dVi[i])*f[nF+j] + 1./Hi*f[nF+i]*(-dVi[j]) -sum1)/Hi *a - Fmi[(2*nF)*(i+nF)+(j+nF)]*1.0*(k1*k1)/(a*a)/Hi *a;;
            u2out[i+nF+(j+nF)*2*nF]=- Fmi[(2*nF)*(i+nF)+(j+nF)]*2.0*Hi/Hi;
			sum1=0.0;

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
        
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		vector<double> RMi;
		RMi = fmet.Riemn(f); 
		vector<double> RMCi;
		RMCi = fmet.Riemn(f);
		
		double sum1=0;
		double sum2=0;
		for(int i=0;i<nF;i++){
		for(int j=0;j<nF;j++){
		sum1=sum1+FMi[(2*nF)*(i)+(j)]f[nF+j]*f[nF+i];}}
		for(int i=0;i<nF;i++){Xi[i] = 2.*(-dVi[i]-3.*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			for(int l=0; l<nF; l++)
			{for(int m=0; m<nF; m++)
			{for(int n=0; n<nF; n++)
			{for(int p=0; p<nF; p++)
				
			sum2 = sum2 - 1./3.*f[nF + k](1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(j+nF)+(m+nF)] 
			+ 1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(j+nF)+(2*nF)*(i+nF)+(m+nF)] *Fmi[(2*nF)*(m)+(n)]*f[nF+n]*Fmi[(2*nF)*(l)+(p)]*f[nF+p])
			- 1./3.*f[nF + i](1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(k+nF)+(2*nF)*(j+nF)+(m+nF)] 
			+ 1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(j+nF)+(2*nF)*(k+nF)+(m+nF)] *Fmi[(2*nF)*(m)+(n)]*f[nF+n]*Fmi[(2*nF)*(l)+(p)]*f[nF+p])
			- 1./3.*f[nF + j](1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(k+nF)+(2*nF)*(i+nF)+(m+nF)] 
			+ 1./2. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(k+nF)+(m+nF)] *Fmi[(2*nF)*(m)+(n)]*f[nF+n]*Fmi[(2*nF)*(l)+(p)]*f[nF+p]);
			
			sum3 = sum3 + 1./18. *( RMCi[(nF)*(nF)*(nF)*(nF)*i+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+k] + RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+i]
			+ RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*i+k] + RMCi[(nF)*(nF)*(nF)*(nF)*i+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+j]
			+ RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*i+j] + RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+i]) * Fmi[(2*nF)*(m)+(n)]*f[nF+n] * Fmi[(2*nF)*(l)+(p)]*f[nF+p];
			
			}}}
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
            - 1./3.*f[nF + k]/32./Hi/Hi/Hi * Xi[i] * Xi[j] * (k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4.
    		+ 1./3.*FMi[(2*nF)*(j+nF)+(k+nF)]*f[nF+i]/2./Hi*(-k2*k2-k3*k3+k1*k1)/a/a/2.
			+ 1./3.*FMi[(2*nF)*(i+nF)+(k+nF)]*f[nF+j]/2./Hi*(-k1*k1-k3*k3+k2*k2)/a/a/2.
			+ 1./3.*FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF+k]/2./Hi*(-k2*k2-k1*k1+k3*k3)/a/a/2.
			+ sum3
			+ sum2;
			sum2=0;
			sum3=0;
            }}}

        return A;
    }

		
    //Calculates B term of action
   vector<double> Bcalc(vector<double> f,vector<double> p, double k1, double k2, double k3,double N)
	{
		vector<double> FMi;
		FMi = fmet.fmetric(f);
		vector<double> RMi;
		RMi = fmet.Riemn(f); 

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
		double sum2=0;
		
		for(int i=0;i<nF;i++){sum1=sum1+f[nF+i]*f[nF+i];}
		for(int i=0;i<nF;i++){Xi[i] = 2.0*(-dVi[i]-3.0*Hi*f[nF+i])+f[nF+i]/Hi*sum1;}
		
        
        for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			for(int l=0; l<nF; l++)
			{for(int n=0; n<nF; n++)
			{for(int p=0; p<nF; p++)
			{
			sum2 = sum 2 + 2./9 *( RMi[(2*nF)*(2*nF)*(2*nF)*(i+nF)+(2*nF)*(2*nF)*(j+nF)+(2*nF)*(k+nF)+(l+nF)] +RMi[(2*nF)*(2*nF)*(2*nF)*(i+nF)+(2*nF)*(2*nF)*(k+nF)+(2*nF)*(j+nF)+(l+nF)] 
			+ RMi[(2*nF)*(2*nF)*(2*nF)*(j+nF)+(2*nF)*(2*nF)*(k+nF)+(2*nF)*(i+nF)+(l+nF)] + RMi[(2*nF)*(2*nF)*(2*nF)*(k+nF)+(2*nF)*(2*nF)*(j+nF)+(2*nF)*(i+nF)+(l+nF)] 
			+ RMi[(2*nF)*(2*nF)*(2*nF)*(k+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(j+nF)+(l+nF)] + RMi[(2*nF)*(2*nF)*(2*nF)*(j+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(k+nF)+(l+nF)] )*Fmi[(2*nF)*(l)+(p)]*f[nF+p])
			}
			}
			}
		
			B[i + j*nF +k* nF*nF] = 1.*f[nF + i]*f[nF+j]*f[nF+k]/4./Hi/Hi
			- 1./2.*f[nF + i] * f[nF + k]/8./Hi/Hi/Hi * Xi[j]
            - 1./2.*f[nF + j] * f[nF + k]/8./Hi/Hi/Hi * Xi[i]
			+ 1./2.*f[nF + i] * f[nF + k]/8./Hi/Hi/Hi * Xi[j]*(k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            + 1./2.*f[nF + j] * f[nF + k]/8./Hi/Hi/Hi * Xi[i]*(k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.
			- FMi[(2*nF)*(j+nF)+(k+nF)]*1.*Xi[i]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k1/k1/2.
			- FMi[(2*nF)*(i+nF)+(k+nF)]*1.*Xi[j]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k2/k2/2.
			+ sum2
			;}
		}}}
        return B;
    }

    //Calculates C term of action
    vector<double> Ccalc(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
		double Hi=H(f,p);
        vector<double> FMi;
		FMi = fmet.fmetric(f);
		
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
			- 1.*f[nF + i] * f[nF+j] *f[nF+k]/8./Hi/Hi/Hi *(k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4. 
			- FMi[(2*nF)*(i+nF)+(j+nF)]*1.*f[nF+k]/2./Hi
			+ FMi[(2*nF)*(j+nF)+(k+nF)]*f[nF+i]/2./Hi*(-k1*k1-k3*k3+k2*k2)/k1/k1/2.
			+ FMi[(2*nF)*(i+nF)+(k+nF)]*f[nF+j]/2./Hi*(-k2*k2-k3*k3+k1*k1)/k2/k2/2.;}
		}}}
        return C;
    }
    

    
	//calculates u3
	vector<double> u(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
        vector<double>  A, B,B2, B3, C, C2,C3;
        double Hi;
		Hi=H(f,p);
        double a =exp(N);
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
			
			u3out[(i)+(nF+j)*2*nF+k*2*nF*2*nF]= -C[i+j*nF+k*nF*nF]/Hi /a;
			u3out[(i)+j*2*nF+(k+nF)*2*nF*2*nF]= -C2[i+k*nF+j*nF*nF]/Hi /a;
			
			u3out[(i)+(j+nF)*2*nF+(k+nF)*2*nF*2*nF]= 0.;
			
			u3out[(nF+i) + j*2*nF + k*2*nF*2*nF]= 3.*A[i+j*nF+k*nF*nF]/Hi *a;
		
			u3out[(nF+i)+(nF+j)*2*nF+k*2*nF*2*nF]=B3[i+k*nF+j*nF*nF]/Hi ;
			u3out[(nF+i)+(j)*2*nF+(k+nF)*2*nF*2*nF]=B2[i+j*nF+k*nF*nF]/Hi ;
			
			u3out[(nF+i)+(j+nF)*2*nF + (k+nF)*2*nF*2*nF]=C3[k+j*nF+i*nF*nF]/Hi /a;

		}}}
        return u3out;
	}


//calculates N1
vector<double> N1(vector<double> f,vector<double> p, double N)
{
    double Hd=Hdot(f);
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
    double Hd=Hdot(f);
    double Hin=H(f,p);
    vector<double> dVi, dVVi;
    vector<double> Nii(2*nF*2*nF);
    double a = exp(N);
	
   	vector<double> FMi;
	FMi = fmet.fmetric(f);
	vector<double> RMi;
	RMi = fmet.Riemn(f); 
	vector<double> CHi;
	CHi = fmet.Chroff(f);
	
    dVi=pot.dV(f,p);
    dVVi=pot.dVV(f,p);
    
    double sum3 = 0.0;
	double sum4 = 0.0;
	double sum5 = 0.0;
    for(int i=0;i<nF;i++){
	for(int j=0;j<nF;j++){
	sum3=sum3+dVi[i]*Fmi[(2*nF)*(i)+(j)]*f[nF+i]/Hin/Hin/Hin;}}
     
     
    double ep = -Hd/Hin/Hin;
    for(int i=0; i<nF; i++){for(int j=0; j<nF; j++){
	for(int l=0; l<nF; l++){
		sum5 = sum5 - 1./4. * dVi[l] * CHi[2*nF)*(l)+(nF)*(i+nf)+(j+nF)];
	for(int m=0; m<nF; m++){
	for(int n=0; n<nF; n++){
		sum6 = sum6 + 1./2. * Fmi[(2*nF)*(i+nF)+(m+nF)]*Fmi[(2*nF)*(l)+(n)]*f[nf+n]*CHi[2*nF)*(m)+(nF)*(l+nf)+(j+nF)];
	for(int p=0; p<nF; p++){
		sum4 = sum4 - 1./4. * RMi[(2*nF)*(2*nF)*(2*nF)*(l+nF)+(2*nF)*(2*nF)*(i+nF)+(2*nF)*(j+nF)+(m+nF)] * Fmi[(2*nF)*(m)+(n)]*f[nF+n]*Fmi[(2*nF)*(l)+(p)]*f[nF+p];
	}
	}
	}
	}
    Nii[i + (j) * 2*nF]= 2./ep/Hin/Hin/6. * (f[nF+i]*f[nF+j] *(-3./2. + 9./2./ep + 3./4.*sum3/ep/ep) +sum4 +sum5 );
    Nii[i + (j+nF) * 2*nF]=2./ep/Hin/Hin/6.*(3./2.*f[i+nF]*f[j+nF]/Hin/ep +1./2. * sum6 ) /a - Fmi[(2*nF)*(i)+(j)]*2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k2*k2) /a;
    Nii[i+nF + (j) * 2*nF]=2./ep/Hin/Hin/6.*(3./2.*f[i+nF]*f[j+nF]/Hin/ep +1./2. * sum6 ) /a - Fmi[(2*nF)*(i)+(j)]*2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k3*k3) /a;
    Nii[i+nF + (j+nF) * 2*nF]=0.;
    }}

    return Nii;

}





};
#endif
