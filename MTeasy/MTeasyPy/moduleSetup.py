#module stepup script

from distutils.core import setup, Extension
import numpy
import time
import os
dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'MTeasyPy.cpp')
f = open(filename,"r")
lines = f.readlines()
f.close()
f = open(filename,"w")

for line in lines:
    if not line.startswith("// Package recompiled"): 
        f.write(line)

    if line.startswith("// Package recompiled"):      
        f.write('// Package recompiled: '+ time.strftime("%c") +'\n')
f.close()        
filename2 = os.path.join(dir, '../MTeasyC/stepper/rkf45.cpp')
dirs = os.path.join(dir, '../MTeasyC')
# don't edit the comment at the end of the setup line below
setup(name="MTeasyPyPS", version="3.0", ext_modules=[Extension("MTeasyPyPS", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup