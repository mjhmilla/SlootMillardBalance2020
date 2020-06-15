import glob
import os
import rbdl
import c3d
import argparse
import sys
from numpy import genfromtxt 
#import scipy as sci
#import csv
#
#import matplotlib.pyplot as plt

def getGroupNameList(reader, groupName, categoryName):
    nameList=[]
    groups = ((k, v) for k, v in reader.groups.items() if isinstance(k, str))
    for key, g in sorted(groups):
        if isinstance(key, str):
            if key == groupName:
                for key,p in sorted(g.params.items()):
                    if key == categoryName and len(p.dimensions) == 2:
                        C, R = p.dimensions
                        for r in range(R):
                            name=p.string_array[r]
                            nameList.append(name.strip())                                                       
    return nameList


dirPythonCode  = os.getcwd()
dirCode        = os.path.join(dirPythonCode,os.pardir)
dirProjectRoot = os.path.join(dirCode,os.pardir)
dirOutputData  = os.path.join(dirProjectRoot, 'outputData')

tag2D = "2D"
isModel2D = False

tagSequence = "motionSequence"

subjectsToProcess = ["E01"]

os.chdir(dirOutputData)



for subject in subjectsToProcess:
    os.chdir(subject)
    print(os.getcwd())
    for luaFile in glob.glob('*.lua'):
        print(luaFile)
        isModel2D=True if luaFile.find(tag2D)>=0 else False
        model = rbdl.loadModel(luaFile, kwargs={"floating_base":True,"verbose":True})

        
        qinit    = np.ndarray(shape(model.q_size), dtype=float, mode="c")
        qres     = np.ndarray(shape(model.q_size), dtype=float, mode="c")

        #body_ids = np.ndarray(shape(, dtype=float, mode="c")
        #body_point_position = np.ndarray[double, ndim=2, mode="c"] 
        #target_pos_position =np.ndarray[double, ndim=2, mode="c"] 
        #np.ndarray[double, ndim=1, mode="c"] qres,

        print("DoF: ", model.q_size)   
        for root,dirs,files in os.walk(os.getcwd()):
            for name in dirs:
                print('\t'+name)
                os.chdir(name)
                c3dFileName = ""
                for file in glob.glob('*.c3d'):
                    if file.find(tag2D) >= 0 and isModel2D == True:
                        c3dFileName=file
                    if file.find(tag2D) == -1 and isModel2D == False:
                        c3dFileName=file
                print('\t'+'\t'+c3dFileName)
                reader=c3d.Reader(open(c3dFileName,'rb'))   

                #Get the names of the c3d marker entries
                pointNames = getGroupNameList(reader, 'POINT', 'LABELS')
                #analogNames= getGroupNameList(reader,'ANALOG','LABELS')

                

                sequenceFileName = ""
                for file in glob.glob('*.csv'):
                    if file.find(tagSequence)==0:
                        sequenceFileName=file
                sequence=genfromtxt(sequenceFileName,delimiter=",")
                print('\t'+'\t'+sequenceFileName)    
                os.chdir(os.pardir)
        os.chdir(os.pardir)

#plotInputData     = True
#plotOutput        = True

#Preprocessing flags
#resampleFrequency = 100
# All of the input data is resampled to the above rate. This step is included
# so that data coming from different equipment, at different sample rates, 
# is properly handled.

#filterFreq        = 7.5;
# The inverse-kinematics data is filtered using a 2nd order Butterworth filter
# in the forwards and backwards direction (so there is no phase introduced).
# This is the 3db frequency of


#Read in the model
# The model name convention follows the one in OpenSim:
# gait: model intended for walking simulations
#    9: DoF
#   12: Number of muscles. In this case torque muscles
#model = rbdl.loadModel("gait912.lua", kwargs={"floating_base":True,"verbose":True})
#print("DoF: ", model.q_size)
#q_size    = model.q_size
#qdot_size = model.qdot_size

