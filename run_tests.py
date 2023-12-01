#!/usr/bin/env python3
import os
from subprocess import check_output
import re
from time import sleep


#
#  Feel free (a.k.a. you have to) to modify this to instrument your code
#

#THREADS = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]
THREADS = [1,2,4,6,8,10,12]
INPUTS = ["nb-100.txt"]
STEPS = 100
THETA = 0.5
TIMESTEP = 0.005

csvs = []
for inp in INPUTS:
    #for loop in LOOPS:
        csv = ["{}".format(inp)]
        for thr in THREADS:
            inputSize = inp.split(".")[0]
            outputName = "test" 
            cmd = "mpirun -np {} ./bin/nbody -i ./input/{} -o ./{} -s {} -t {} -d {}".format(
                thr, inp, outputName, STEPS, THETA, TIMESTEP)
            out = check_output(cmd, shell=True).decode("ascii")
            #m = re.search("*", out)
            #print(out)
            csv.append(out)
            #if m is not None:
            #    time = m.group(1)
            #    csv.append(time)
            #compareName= "results/result_%s_t0.txt" % inputSize
            #cmd = "diff {} {}".format(compareName,outputName)
            #out = check_output(cmd, shell=True).decode("ascii")
            #print(out)

        csvs.append(csv)
        sleep(0.5)

header = ["seconds"] + [str(x) for x in THREADS]

print("\n")
print(", ".join(header))
for csv in csvs:
    print (", ".join(csv))



