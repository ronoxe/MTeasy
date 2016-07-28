## README ##

For license/copyright information see LICENSE.txt which was distributed with this file. 

If you use MTeasy for any published work you are kindly asked to cite the related papers: 


## MomentTransportEasy ##

This is a simple code in c++, together with python files to compile a python module which 
wraps the c++ code. The use of python enhances functionality and convenience, leading to an 
easy to use resource for inflationary cosmology. The code solves the transport equations for
the correlation functions (also known as moments) of canonical multi-field inflationary 
cosmology. This allows the power-spectrum and bispectrum of multi-field inflation to be 
calculated. The original papers on this system of equations refer to them as the moment transport 
equations. (although the original papers were limited to "super-horizon" effects, 
while the present system evolves the moments from sub-horizon scales).

The c++ code has a simple object orientated structure, but currently this is not exposed by 
the python module, which simply provides key routines. The c++ code folder also includes an rk45 
integrator routine written by John Burkardt and distributed under a GNU LGPL license described
within the relevant files. 

There are no dependancies external to the folders provided except for a working python 
installation -- this is deliberate to make using the code as easy as possible to use.

### How do I set up the python module? ###

Here follow concise instructions (we hope to add a illustrated pdf guide).

Although the c++ code can be used on its own perfectly well, the intention is that a user will
interact with the code using python. This allows the power of python math, science, symbolic and
plotting packages to be utilised to analyse an inflationary model within an integrated 
environment. This is the path described below, and therefore a working python installation is 
needed. Because the c++ code is compiled into a python module, there is no (or negligible) 
change in speed with respect to running the c++ routines natively. 

* Python:  For convenience the author recommends downloading a complete python distribution, for  
  example pythonxy, enthought, continuum anaconda, or similar distribution, which come with all 
  the python packages used by this code. They also come with interactive scientific development 
  environments (Spyder provided by continuum anaconda is used by the author). The packages 
  currently used by MTeasyPy or examples include matplotlib, scipy, sympy, distutils. 

* Download code: Once you have a working installation of python, with the packages mentioned
  (which the distributions sort for you), download the MTeasy folder containing MTeasyPy, 
  MTeasyC and runData subfolders. Maintain this file structure wherever you place the MTeasy 
  folder. You are ready to go!

* Compiling the python module: The only real setup required is to compile the python module for 
  the inflationary module you want to analysis. The first step is to edit the potential. This 
  is defined in potentialSetup.py in the python folder. This takes a user supplied potential, 
  and when the file is run automattically produces a potential class file in the c++ folder, 
  which will be compiled into the c++ code when the python module is compiled. Which is the next 
  step. This task is completed by the moduleSetup.py file. For convenience, both this file and 
  potentialSetup.py are called from MTeasySetup.py, which is therefore the only file which needs 
  to be called once the user is happy the potential they want to analysis is defined in 
  potentialSetup.py. moduleSetup.py uses MTeasyPy.cpp to compile a python module consisting of 
  a number of functions which wrap c++ code, those functions are described now: 

### How do use the python module? ###
The python module gives access to a number of functions:

  - MTeasyPy.nF takes no arguments and returns the number of fields relevant to the model.  
  
  - MTeasyPy.V(fields) takes an numpy array of length nF (the number of fields) containing a set 
    of field values and returns the value of the potential.

  - MTeasyPy.dV(fields,dV) takes an numpy array of length nF (the number of fields) containing a 
    set, and a second bumpy array of length nF (dV), and sets dV to contain the first
    derivatives of the potential wrt to each field.

  - MTeasyPy.H(fields_dotfields) takes an numpy array of length 2*nF (twice the number of 
    fields) containing a set of field values followed by the field's velocities in cosmic time 
    (field derivative wrt cosmic time) and returns the value of the Hubble rate.

  - MTeasyPy.backEvolve(Narray, fields_dotfields) takes a numpy array of times in efolds (N) 
    that output is required at (Narray), this must start with the initial N and finish with the 
    final N, and also takes the initial conditions at the initial N (field values followed by 
    field velocities) as a numpy array. It produces a file in the runData folder called 
    back.dat. This file contains the fields, and field velocities at the times requested by 
    Narray. The format is 1 + 2nF columns, with the first column the times (Narray), the next 
    columns the field values and field velocity values at those times. This file is written over 
    every time this function runs, so if data needs to be stored it should be uploaded to python 
    and stored elsewhere between calls (see the examples).

  - MTeasyPy.sigEvolve(Narray, k, fields_dotfields, full) takes times in efolds that output is   
    required at (must start at initial time and end at the final time), a fourier mode value 
    (k), and the initial conditions at the initial time, and an integer "full" set to 0 or 1 
    (if other defaults to 1). If full = 1 it produces a file in the runData folder called 
    sig.dat, and a second file zz.dat. If full=0 it just produces zz.dat. sig.dat contains the 
    field, and field velocities at the time steps requested, and the real parts of the sigma_ij 
    (field space power spectrum) matrix, where i and j run over fields and then field 
    velocities. There are [1 + 2nF + 2nF x 2nF] columns, with the first column the times, the 
    next 2nF columns the corresponding field values and field velocity values, and 
    the next 2nF x 2nF are the sigma matrix. The convention is sig^r_ij is the [1 + 2nF + i + 
    2nF x (j-1)]th column of the file. This file is written over every time the this function 
    runs, so if data needs to be stored it should be uploaded to python and stored elsewhere 
    between calls (see examples). The zz.dat file contains 2 columns, the first is the times, 
    and the second contains P_zeta(k) (the power spectrum of zeta).

 -  MTeasyPy.alpEvolve(Narray, k1, k2, k3, fields_dotfields, full) takes times in efolds that 
    output is required at (must start at initial time and end at the final time), three fourier  
    mode values (k1, k2 , k3), the initial conditions, and an integer "full" set to 0 or 1 (if 
    other defaults to 1). If full = 1 it produces a file in the runData folder called 
    alpha.dat, and a second file zzz.dat. If full=0 it just produces zzz.dat. alpha.dat 
    contains the field, and field velocities at the time steps requested, all the relevant sigma 
    data, and the field space bispectrum alpha. There are [1 + 2nF + 6 x 2nF x 2nF + 2nF x 
    2nF x 2nF] columns, with the first column the times, the next 2*nF columns the 
    corresponding field values and field velocity values, and the next 2nF x 2nF are the real 
    parts of sigma_ij(k1) in the same numbering convention as above. Then the real part of 
    sigma_ij(k2) and then the real parts of sigma_ij(k3), the following 3 x 2nF x 2nF columns are 
    the imaginary parts of the sigma_ij(k1), sigma_ij(k2) and sigma_ij(k3) matricies. So for 
    example if one wanted access to the sigma^i(k2)_ij that would be the 
    [1 + 2nF + 4 x 2nF x 2nF + i + 2nF x (j-1)]th column of the file. The final 
    2nF x 2nF x 2nF columns of this file correspond to the alpha_ijl(k1,k2,k3) matrix, such that     
    the corresponding columns would be the [1 + 2nF + 6 x 2nF x 2nF + i + 2nF x (j-1) = 
    2nF x 2nF x 2nF x (l-1)]th columns. This file is written over every time, this function runs, so 
    if data needs to be stored, it should uploaded to python and stored elsewhere between calls 
    (see examples). The zzz.dat file contains 5 columns, the first is the times, and the second 
    contains P_zeta(k1), the next P_zeta(k2) and the fourth P_zeta(k3). The fifth column is the 
    bispectrum B_zeta(k1,k2,k3).