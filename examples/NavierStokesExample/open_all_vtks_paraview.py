#### import the simple module from the paraview
import os
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

stringtosearch = "velocity." # CHANGE THIS VARIABLE TO LOAD OTHER VARIABLE

print(os.curdir)
dir = os.getcwd() + "/solutions/"
renderView1 = GetActiveViewOrCreate('RenderView')
count = 0

def getStamp(elem):
    firstdot = elem.find(".")+1
    print(elem[firstdot:firstdot+5])
    return int(elem[firstdot:firstdot+5])

listdir = os.listdir(dir)
lists = []
keep = 1
while (keep):
    thereisone = 0
    newlist = []
    for file in listdir:
        if stringtosearch in file and "block" + str(count) + "_" in file and ".pvtu" in file:
            thereisone = 1
            newlist.append(dir + file)
    lists.append(newlist)
    if thereisone is 0:
        keep = 0
        del lists[-1]
    else:
        count = count + 1

animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
renderView1 = GetActiveViewOrCreate('RenderView')

for i in range(0,count):
    lists[i].sort(key=getStamp)
    block = XMLPartitionedUnstructuredGridReader(FileName=lists[i])
    block_disp = Show(block, renderView1)
