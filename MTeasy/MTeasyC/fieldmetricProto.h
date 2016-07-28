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
    
public:
    fieldmetric()
   {
// #FP
	
   }
	
	
	//calculates fieldmetic()
	vector<double> fmetric(vector<double> f)
	{
		vector<double> FM((2*nF)*(2*nF)) ;
        
// metric

         return FM;
	}
	
	
	
	//calculates ChristoffelSymbole()
	vector<double> Chroff(vector<double> f)
	{
		vector<double> CS((2*nF)*(nF)*(nF));
	
// Christoffel
        
		return CS;
	}
    
	
	
	// calculates RiemannTensor()
	vector<double> Riemn(vector<double> f)
	{
		vector<double> RM((2*nF)*(2*nF)*(2*nF)*(2*nF));
		
// Riemann
     
        return RM;
	}

	// calculates RiemannTensor() covariant derivatives
	vector<double> Riemncd(vector<double> f)
	{
		vector<double> RMcd((nF)*(nF)*(nF)*(nF)*(nF));
		
// Riemanncd
     
        return RMcd;
	}
    
    int getnF()
    {
        return nF;
    }
    


};
#endif

