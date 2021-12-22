#!/usr/bin/env python3

import sys
import subprocess
import numpy as np
from vtk import vtkXMLUnstructuredGridReader
from vtk.util import numpy_support as VN

if len(sys.argv) != 4 :
    print('Invalid number of arguments!')
    print('Usage: ' + sys.argv[0] + ' LADDSExe cfg reference.vtu')
    sys.exit(1)

# path to the ladds executable
exe_path=sys.argv[1]
# path to the yaml input file
yaml_path=sys.argv[2]

# run the simulation and capture the output. Abort if the simulation failed.
sim = subprocess.Popen([exe_path, yaml_path])
output, err = sim.communicate()
exit_code = sim.wait()
if exit_code != 0 :
    sys.exit('An error occured during the simulation')

output_path='output_99.vtu'
# path to the given reference vtu
reference_path=sys.argv[3]

# read all positions and velocities data from output and reference
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(output_path)
reader.Update()
output = reader.GetOutput()

v_out = VN.vtk_to_numpy(output.GetPointData().GetArray('velocity')).astype('double')
r_out = VN.vtk_to_numpy(output.GetPoints().GetData()).astype('double')

reader.SetFileName(reference_path)
reader.Update()
reference = reader.GetOutput()

v_ref = VN.vtk_to_numpy(reference.GetPointData().GetArray('velocity')).astype('double')
r_ref = VN.vtk_to_numpy(reference.GetPoints().GetData()).astype('double')

# perform similarity checks on all loaded data
# Tolerances are empirically determined as tight as possible.
allTestsOk = True
if not (len(r_out) == len(r_ref)) :
    print('Number of particles differs!')
    allTestsOk = False
if not np.allclose(v_out, v_ref, atol=1, rtol=0) :
    print('Velocities differ!')
    allTestsOk = False
if not np.allclose(r_out, r_ref, atol=4, rtol=0) :
    print('Positions differ!')
    allTestsOk = False

# communicate test result via exit status
if not allTestsOk :
    sys.exit(-1)
