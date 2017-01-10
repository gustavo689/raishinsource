#!/usr/bin/env bash
#
# Script that tells you:
#
# § [ ] when the simulation started running
# § [x] how long (simulation `t`) it has ran so far
# § [x] storage pressure
# 

# how long has it been running?
qstat

# current status of simulation 
grep "time=" grmhd.out | tail -n1

# storage pressure
echo -n "Storage pressure = "
du  -hc *dat | tail -n1