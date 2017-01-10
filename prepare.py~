#!/usr/bin/env python
#
# Main script for submitting a RAISHIN job to a cluster.
#
# Usage: 
#   prepare.py <NDIM> <NCORES> <NCPN>
#
#	where
#		NDIM: number of dimensions
#		NCORES: number of CPUs you want to use (MPI) 
#       NCPN: number of cores per node in your cluster,
#			0 if you are running in a single machine
#
# Run this in the folder which will store the simulation results.
#
# Example: 
# >>> ./scripts/prepare.py 2 400 24
#
# Requirements:
# - RAISHIN with configuration files properly set up
# - nmmn python module installed (https://github.com/rsnmmn/nmmn)
#
# TODO
# - [x] adjust number of CPUs to number of nodes/cores in the cluster
# - [x] prepares files
# - [ ] update job file
# - [ ] correct automatic determination of CPUs for 3D case
# - [ ] use configuration file to specify parameters instead of 
#       pram.f90


import math, sys
import shutil, subprocess, re, datetime, os
import nmmn.misc

# Number of cores in each node in your cluster.
# For alphacrucis, 24 cores in 96 nodes.
#nppn=24

# ASCII art from http://ascii.co.uk/art/lightning
print("                    ,/")
print("                  ,'/")
print("                ,' /")
print("              ,'  /_____,")
print("            .'____    ,'    RAISHIN")
print("                 /  ,'")
print("                / ,'")
print("               /,'")
print("              /'")



here      = os.getcwd()

############# COMMAND-LINE ARGUMENTS
# check if there were command-line arguments
if len(sys.argv)==4: # there are command-line arguments that were actually typed
	ndim= int(sys.argv[1])
	ncpu = int(sys.argv[2])
	nppn=int(sys.argv[3])
else: # there is nothing
	print('Usage: \npreparecode.py <NDIM> <NCORES> <NCORES/NODE>')
	sys.exit(0)




############# CPU SPLITTING 
# determines splitting of computational mesh for MPI
# careful: imperfect splitting, can be improved
if ndim==1:	# 1D 
	ni=ncpu
	nj=1
	nk=1
elif ndim==2:	# 2D in i,k
	if nppn==0: # running in single machine
		ni=int(round(math.sqrt(ncpu)))
		nj=1
		nk=int(round(ncpu/float(ni)))
		ncpu=ni*nk  # not necessarily equal to the input ncpu
	else: # running in a cluster
		ni=int(round(ncpu/nppn))
		nj=1
		nk=nppn
		ncpu=ni*nk  # not necessarily equal to the input ncpu
elif ndim==3:	# 3D, NEED TO IMPROVE
	ni=int(round(ncpu**(1./3.)))
	nj=ni
	nk=int(round(ncpu/(float(ni)*nj)))
	ncpu=ni*nj*nk   # not necessarily equal to the input ncpu

print(str(ncpu)+" CPUs will be used")


############# READS/UPDATES FILES
# http://stackoverflow.com/questions/39086/search-and-replace-a-line-in-a-file-in-python
#
# input will be superseded eventually by a python configuration file
# instead of pram.f90

shutil.copy("pram.f90","pram.f90.bak") # backup copy
f = open("pram.f90.bak","r")
newf=open("pram.f90", "w")	# replaces parameter file

for line in f:
	if 'iprocs' in line:	# sets right number of CPUs
		newf.write("  integer, parameter :: iprocs="+str(ni)+", jprocs="+str(nj)+", kprocs="+str(nk)+" !- CPU number in i-,j-, and  k- direction\n")
	else:
		newf.write(line)

	# gets tmax
	if re.search(r' tmax=\d+\.\d+d', line):
		tmax=re.search(r'\d+\.\d+d0', line).group()
		tmax=float(re.sub(r"d", "e", tmax)) # replaces d with e in the number for use in python

	# gets nshot
	if re.search(r' nshot=\d+', line):
		nshot=re.search(r'\d+', line).group()
		
f.close()
newf.close()



shutil.copy("convert_vtk2dn1.f90","convert_vtk2dn1.f90.bak") # backup copy
f = open("convert_vtk2dn1.f90.bak","r")
newf=open("convert_vtk2dn1.f90", "w")	# replaces parameter file

for line in f:
	if 'integer, parameter :: ns' in line:	# sets right number of VTK snapshots
		newf.write("  integer, parameter :: ns=0, ne="+nshot+" ! start and end data file number\n")
	else:
		newf.write(line)
		
f.close()
newf.close()


# opens this file only if running in a cluster
if nppn!=0:
	shutil.copy("mpirun.job","mpirun.job.bak") # backup copy
	f = open("mpirun.job.bak","r")
	newf=open("mpirun.job", "w")	# replaces parameter file

	for line in f:
		if 'nodes=' in line:	
			newf.write("#PBS -l nodes="+str(ni)+":ppn="+str(nppn)+"\n")
		elif 'mpirun' in line:	
			newf.write("mpirun -n "+ str(ncpu) +" -machinefile $PBS_NODEFILE /sto/home/gustavo/raishin/runs/"+here+"/xgrmhd.exe\n")
		else:
			newf.write(line)

		
	f.close()
	newf.close()


print("Updated configuration files \n")



############# COMPILATION
shutil.copy("Makefile_xgrmhd","Makefile")

# removes previous datafiles
subprocess.Popen("rm -f restart*.outdat struct*.outdat ok*.vtk", shell=True)

print("Source files ready. Now proceed: \n    send source to cluster and ssh \n    make clean \n    make \n    qsub mpirun.job ")
