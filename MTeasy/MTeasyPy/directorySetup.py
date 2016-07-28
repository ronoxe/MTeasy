import os
dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'MTeasyPy.cpp')
f = open(filename,"r")
lines = f.readlines()
f.close()
filename = os.path.join(dir, 'MTeasyPy.cpp')
f = open(filename,"w")
for line in lines:
    if not  line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith("//model\n") and not line.endswith("//stepper\n"):   
       f.write(line)
    if line.endswith("//evolve\n"):
       fileT = os.path.join(dir, '../MTeasyC/evolve.h')
       f.write('#include' + '"'+ fileT +'"' + '//evolve' +'\n')
    if line.endswith("//moments\n"):      
       fileT = os.path.join(dir, '../MTeasyC/moments.h')
       f.write('#include' + '"'+ fileT +'"' + '//moments' +'\n')
    if line.endswith("//model\n"):      
       fileT = os.path.join(dir, '../MTeasyC/model.h')
       f.write('#include' + '"'+ fileT +'"' + '//model' +'\n')
    if line.endswith("//stepper\n"):      
       fileT = os.path.join(dir, '../MTeasyC/stepper/rkf45.hpp')
       f.write('#include' + '"'+ fileT +'"' + '//stepper' +'\n')

f.close()

