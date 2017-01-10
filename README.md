Running and producing data from RAISHIN GRMHD simulations using a cluster
========================================

Here are the basic steps required in order to run and produce data for RAISHIN.

1. Setup your model carefully (`pram.f90`, `mdgrmhd.f90` etc)

2. Run `prepare.py` to configure the appropriate number of cores and nodes for your MPI run. For example, prepare source for 2D, 512 cores, 24 cores/node:

```python
prepare.py 2 512 24
```

Preparing the same run in a desktop machine:

    prepare.py 2 8 0

3. If running in a cluster: I simply run the bash script `toac.sh` to send it from my workstation to the alphacrucis cluster.

4. `make` and `qsub mpirun.job` to submit the job. Wait. Use `status.sh` once in a while to check the status of your simulation.

5. Once finished, copy the raw data files to the workstation where the visualization and analysis will be carried out

6. Run `convert.py` to convert the raw data files to VTK format



