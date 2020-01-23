import sys, string, os
import numpy as np

if (sys.argc == 2):
    logfolder = sys.argv[1]
else:
    logfolder = "logs/"

outp = np.zeros([4,13])

for i in range(3,16):
    ndofs = 0
    setuptime = 0
    niterations = -1
    solvetime = 0
    with open(logfolder + "log" + str(i) + ".txt", "rt") as fin:
        donesfound = 0;
        for line in fin:
            stringtofind = "[GlobalAssembler] number of primal + dual dofs = "
            res = line.find(stringtofind)
            if (res != -1):
                ndofs = int(line[res + len(stringtofind):-1])

            stringtofind = "35mdone, in "
            res = line.find(stringtofind)
            if (res != -1):
                donesfound = donesfound + 1
                if (donesfound == 1):
                    setuptime = float(line[res + len(stringtofind):-9])
                elif (donesfound == 2):
                    solvetime = float(line[res + len(stringtofind):-9])

            stringtofind = "Iter "
            res = line.find(stringtofind)
            if (res != -1):
                niterations = niterations + 1

    print(i-3)
    outp[0,i-3] = ndofs
    outp[1,i-3] = setuptime
    outp[2,i-3] = niterations
    outp[3,i-3] = solvetime

np.savetxt(logfolder + "output.csv", outp, delimiter=",")
