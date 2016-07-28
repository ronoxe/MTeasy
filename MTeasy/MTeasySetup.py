import sympy as sym
import numpy as np
from gravipy import *
import math
import subprocess
import sys
import os
import shutil



def directory():
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'MTeasyPy/MTeasyPy.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith("//model\n") and not line.endswith("//stepper\n"):
            f.write(line)
        if line.endswith("//evolve\n"):
            fileT = os.path.join(dir, 'MTeasyC/evolve.h')
            f.write('#include' + '"'+ fileT +'"' + '//evolve' +'\n')
        if line.endswith("//moments\n"):
            fileT = os.path.join(dir, 'MTeasyC/moments.h')
            f.write('#include' + '"'+ fileT +'"' + '//moments' +'\n')
        if line.endswith("//model\n"):
            fileT = os.path.join(dir, 'MTeasyC/model.h')
            f.write('#include' + '"'+ fileT +'"' + '//model' +'\n')
        if line.endswith("//stepper\n"):
            fileT = os.path.join(dir, 'MTeasyC/stepper/rkf45.hpp')
            f.write('#include' + '"'+ fileT +'"' + '//stepper' +'\n')

    f.close()

def pathSet():
    dir = os.path.dirname(__file__)
    path1 = os.path.join(dir, 'MTeasyPy/lib/python/')
    path2 = os.path.join(dir, 'MTeasyPyScripts/')
    sys.path.append(path1)
    sys.path.append(path2)
 
def compile():   
    directory()
    name = ""
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'MTeasyPy/')
    filename1 = os.path.join(dir, 'MTeasyPy/moduleSetup.py')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup"):
            f.write(line)
        if line.endswith("#setup"):
            f.write('setup(name="MTeasyPy'+name+'", version="3.0", ext_modules=[Extension("MTeasyPy'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup')
    f.close()    
    
    filename = os.path.join(dir, 'MTeasyPy/MTeasyPy.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            print "here2"
            f.write('static PyMethodDef MTeasyPy'+name+'_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, MTeasyPy_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, MTeasyPy_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, MTeasyPy_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, MTeasyPy_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, MTeasyPy_docs},  {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, MTeasyPy_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, MTeasyPy_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, MTeasyPy_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//initFunc\n"):
            f.write('void initMTeasyPy'+name+'(void)    {        Py_InitModule3("MTeasyPy'+name+'", MTeasyPy'+name+'_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')
    f.close()    
    
    subprocess.call(["python", "/MTeasyPy/moduleSetup.py", "install", "--home="+location],cwd=location)
    sys.path.append(location+"/lib/python/")
    sys.path.append(location+"../MTeasyPyScripts")
    shutil.rmtree(location+"/build/")

def compileName(name):    
    directory()
    dir = os.path.dirname(__file__)
    print "here"
    filename1 = os.path.join(dir, 'MTeasyPy/moduleSetup.py')
    location = os.path.join(dir, 'MTeasyPy/')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup"):
            f.write(line)
        if line.endswith("#setup"):
            f.write('setup(name="MTeasyPy'+name+'", version="3.0", ext_modules=[Extension("MTeasyPy'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup')
    f.close()    
    
    filename = os.path.join(dir, 'MTeasyPy/MTeasyPy.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write('static PyMethodDef MTeasyPy'+name+'_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, MTeasyPy_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, MTeasyPy_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, MTeasyPy_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, MTeasyPy_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, MTeasyPy_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, MTeasyPy_docs} , {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, MTeasyPy_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, MTeasyPy_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, MTeasyPy_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//initFunc\n"):
            f.write('void initMTeasyPy'+name+'(void)    {        Py_InitModule3("MTeasyPy'+name+'", MTeasyPy'+name+'_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')
    f.close()

    subprocess.call(["python", filename1, "install", "--home=" + location],cwd=location)
    sys.path.append(location+"/lib/python/")
    sys.path.append(location+"../MTeasyPyScripts")
    shutil.rmtree(location+"/build/")

def deleteModule(name):
    location = os.path.join(dir, 'MTeasyPy/')
    os.remove(location+"/lib/python/MTeasyPy"+name+".so")    
    os.remove(location+"/lib/python/MTeasyPy"+name+"-3.0-py2.7.egg-info")    
 
def tol(rtol, atol):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'MTeasyPy/MTeasyPy.cpp')
    f = open(filename,"r")  

    lines = f.readlines()
    f.close()
   
    f = open(filename,"w")  
    for line in lines:
        if not  line.endswith("//Tols\n"):
            f.write(line)
        if line.endswith("//Tols\n"):
            f.write('    double rtol='+str(rtol)+', atol='+str(atol)+';//tols\n')
    f.close()

def potential(G,V,nF,nP):
    f=sym.symarray('f',nF)
    p=sym.symarray('p',nP)
    
    vd=sym.symarray('vd',nF)
    vdd=sym.symarray('vdd',nF*nF)
    vddd=sym.symarray('vddd',nF*nF*nF)
    g, Ga, Ri, Rm =fieldmetric(G,nF)
    FMP=0

    for i in range(nF):
        vd[i] = sym.simplify(V.diff(f[i]))
    for i in range(nF):
        for j in range(nF):
            for l in range(nF):
                FMP=FMP-Ga(-(l+1),i+1,j+1) * vd[l]
                vdd[i+j*nF] = sym.simplify(V.diff(f[i]).diff(f[j]))+FMP
                FMP=0
    for i in range(nF):
        for j in range(nF):
            for k in range(nF):
                for l in range(nF):
                    FMP=FMP-Ga(-(l+1),k+1,i+1)*vdd[l+j*nF] - Ga(-(l+1),k+1,j+1)*vdd[i+l*nF]
                vddd[i+j*nF+k*nF*nF] = sym.simplify(V.diff(f[i]).diff(f[j]).diff(f[k]))
                FMP=0

    import os
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'MTeasyC/potentialProto.h')
    filename2 = os.path.join(dir, 'MTeasyC/potential.h')
    f = open(filename1, 'r')
    g = open(filename2, 'w')


    for line in f:
        g.write(line)
        if line == "// #FP\n":
            g.write('nF='+str(nF)+';\n'+'nP='+str(nP)+';\n')
        #if line == "// pdef\n":
            #for i in range(nP):
            #    g.write('p['+str(i)+']='+str(pp[i])+';\n')
        if line == "// Pot\n":
            expr=str(sym.ccode(V))
            for l in range(max(nP,nF)): 
                expr=expr.replace("_"+str(l),"["+str(l)+"]")
            g.write('\n sum='+str(expr)+';\n')
        if line == "// dPot\n":
            for i in range(nF):
                expr=str(sym.ccode(vd[i]))
                for l in range(max(nF,nP)): 
                    expr=expr.replace("_"+str(l),"["+str(l)+"]")
                g.write('\n sum['+str(i)+']='+str(expr)+';\n')    

        if line == "// ddPot\n":
              for i in range(nF):
                  for j in range(nF):
                      expr=str(sym.ccode(vdd[i+nF*j]))
                      for l in range(max(nF,nP)): 
                          expr=expr.replace("_"+str(l),"["+str(l)+"]")
                      g.write('\n sum['+str(i)+'+nF*'+str(j)+']='+str(expr)+';\n')    
        if line == "// dddPot\n":
              for i in range(nF):
                  for j in range(nF):
                      for k in range(nF):
                          expr=str(sym.ccode(vddd[i+nF*j+nF*nF*k]))
                          for l in range(max(nF,nP)): 
                              expr=expr.replace("_"+str(l),"["+str(l)+"]")
                          g.write('\n sum['+str(i)+'+nF*'+str(j)+'+nF*nF*'+str(k)+']='+str(expr)+';\n')    
                          
       
    g.close()
    f.close()

def fieldmetric(G,nF):
    f=symarray('f',nF)
    COR = Coordinates('\chi', f)
    g = MetricTensor('g',COR , G)
    Ga = Christoffel('Ga', g)
    Ri = Ricci('Ri', g)
    Rm = Riemann('Rm',g)
    import os
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'MTeasyC/fieldmetricProto.h')
    filename2 = os.path.join(dir, 'MTeasyC/fieldmetric.h')
    e = open(filename1, 'r')
    h = open(filename2, 'w')


    for line in e:
        h.write(line)
        if line == "// #FP\n":
            h.write('nF='+str(nF)+';\n')
        if line == "// metric\n":
            for i in  range(2*nF):
                for j in  range(2*nF):
                    if i<nF:
                        ii=-i-1
                    else:
                        ii=i-(nF-1)
                    if j<nF:
                        jj=-j-1
                    else:
                        jj=j-(nF-1)
                    expr=str(ccode(g(ii,jj)))
                    for m in range(nF): 
                        expr=expr.replace("_"+str(m),"["+str(m)+"]")
                    h.write('\n FM['+str((2*nF)*i+j)+']='+str(expr)+';\n')


        if line == "// Christoffel\n":
            for i in range(2*nF):
                for j in range(nF):
                    for k in range(nF):
                        if i<nF:
                            ii=-i-1
                        else:
                            ii=i-(nF-1)
                        jj=j+1
                        kk=k+1
                        expr=str(ccode(Ga(ii,jj,kk)))
                        for m in range(nF): 
                            expr=expr.replace("_"+str(m),"["+str(m)+"]")
                        h.write('\n CS['+str((2*nF)*i+(nF)*j+k)+']='+str(expr)+';\n')
        
    
        if line == "// Riemann\n":
                for i in range(2*nF):
                    for j in range(2*nF):
                        for k in range(2*nF):
                            for l in range(2*nF):
                                if i<nF:
                                    ii=-i-1
                                else:
                                    ii=i-(nF-1)
                                if j<nF:
                                    jj=-j-1
                                else:
                                    jj=j-(nF-1)
                                if k<nF:
                                    kk=-k-1
                                else:
                                    kk=k-(nF-1)
                                if l<nF:
                                    ll=-l-1
                                else:
                                    ll=l-(nF-1)
                                expr=str(ccode(Rm(ii,jj,kk,ll)))
                                for m in range(nF): 
                                    expr=expr.replace("_"+str(m),"["+str(m)+"]")
                                h.write('\n RM['+str((2*nF)*(2*nF)*(2*nF)*i+(2*nF)*(2*nF)*j+(2*nF)*k+l)+']='+str(expr)+';\n')
                          
       
        if line == "// Riemanncd\n":
                for i in range(nF):
                    for j in range(nF):
                        for k in range(nF):
                            for l in range(nF):
                                for m in range(nF):    
                                    ii=i+1
                                    jj=j+1
                                    kk=k+1
                                    ll=l+1
                                    mm=m+1
                                    expr=str(ccode(Rm.covariantD(ii,jj,kk,ll,mm)))
                                    for n in range(nF): 
                                        expr=expr.replace("_"+str(n),"["+str(n)+"]")
                                    h.write('\n RMcd['+str((nF)*(nF)*(nF)*(nF)*i+(nF)*(nF)*(nF)*j+(nF)*(nF)*k+(nF)*l+m)+']='+str(expr)+';\n')

    h.close()
    e.close()
    return g, Ga, Ri, Rm