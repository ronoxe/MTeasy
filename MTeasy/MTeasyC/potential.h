#ifndef POTENTIAL_H  // Prevents the class being re-defined
#define POTENTIAL_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;


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

 sum=0.5*pow(f[0], 2)*pow(p[0], 2) + 0.5*pow(f[1], 2)*pow(p[1], 2) + 0.5*pow(f[2], 2)*pow(p[2], 2);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF);
	
// dPot

 sum[0]=1.0*f[0]*pow(p[0], 2);

 sum[1]=1.0*f[1]*pow(p[1], 2);

 sum[2]=1.0*f[2]*pow(p[2], 2);
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF);
		
// ddPot

 sum[0+nF*0]=1.0*pow(p[0], 2);

 sum[0+nF*1]=0;

 sum[0+nF*2]=0;

 sum[1+nF*0]=0;

 sum[1+nF*1]=1.0*pow(p[1], 2);

 sum[1+nF*2]=0;

 sum[2+nF*0]=0;

 sum[2+nF*1]=0;

 sum[2+nF*2]=1.0*pow(p[2], 2);
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF);
// dddPot

 sum[0+nF*0+nF*nF*0]=0;

 sum[0+nF*0+nF*nF*1]=0;

 sum[0+nF*0+nF*nF*2]=0;

 sum[0+nF*1+nF*nF*0]=0;

 sum[0+nF*1+nF*nF*1]=0;

 sum[0+nF*1+nF*nF*2]=0;

 sum[0+nF*2+nF*nF*0]=0;

 sum[0+nF*2+nF*nF*1]=0;

 sum[0+nF*2+nF*nF*2]=0;

 sum[1+nF*0+nF*nF*0]=0;

 sum[1+nF*0+nF*nF*1]=0;

 sum[1+nF*0+nF*nF*2]=0;

 sum[1+nF*1+nF*nF*0]=0;

 sum[1+nF*1+nF*nF*1]=0;

 sum[1+nF*1+nF*nF*2]=0;

 sum[1+nF*2+nF*nF*0]=0;

 sum[1+nF*2+nF*nF*1]=0;

 sum[1+nF*2+nF*nF*2]=0;

 sum[2+nF*0+nF*nF*0]=0;

 sum[2+nF*0+nF*nF*1]=0;

 sum[2+nF*0+nF*nF*2]=0;

 sum[2+nF*1+nF*nF*0]=0;

 sum[2+nF*1+nF*nF*1]=0;

 sum[2+nF*1+nF*nF*2]=0;

 sum[2+nF*2+nF*nF*0]=0;

 sum[2+nF*2+nF*nF*1]=0;

 sum[2+nF*2+nF*nF*2]=0;
       
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