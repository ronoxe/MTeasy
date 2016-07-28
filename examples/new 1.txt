void evolveSig( double N,  double yin[], double yp[], double paramsIn[])
{
    cout << "Quack" << endl;
	//fieldmetric fmet;
    model m;
    double k;
    int nP = m.getnP();
    vector<double> p(paramsIn,paramsIn+nP);
    k=paramsIn[nP];
    int nF=m.getnF();
    vector<double> fields(yin, yin+2*nF);
    vector<double> u1=m.u(fields,p);
    vector<double> u2=m.u(fields,p,k,N);
    //vector<double> CHR;
	//CHR = fmet.Chroff(f);
    for(int i=0;i<2*nF;i++){yp[i] = u1[i];}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)// use this dummy index instead of o...
            {	
				double sum2=0.0;
				for(int n=0;n<2*nF;n++){
					for(int o=0;o<2*nF;o++){
						//sum2= sum2 + f[nF+o]*CHR[(2*nF)*(nF)*(n)+(nF)*(o)+i]*yin[2*nF+n+2*nF*j] + f[nF+o]*CHR[(2*nF)*(nF)*(n)+(nF)*(o)+j]*yin[2*nF+i+2*nF*n];
						//sum2= sum2 + CHR[(2*nF)*(nF)*(n)+(nF)*(o)+i]*yin[2*nF+n+2*nF*j] + CHR[(2*nF)*(nF)*(n)+(nF)*(o)+j]*yin[2*nF+i+2*nF*n];

					}
					
				}
                sum = sum + sum2 + u2[i+m*2*nF]*yin[2*nF+m+2*nF*j] + u2[j+m*2*nF]*yin[2*nF+m+2*nF*i];
            }
            yp[2*nF+i+2*nF*j]=sum;
        }}
}