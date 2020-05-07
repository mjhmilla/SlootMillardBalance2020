import glob
import os
import rbdl
import c3d
from numpy import genfromtxt 
#import scipy as sci
#import csv
#import c3d
#
#import matplotlib.pyplot as plt

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
                c3dReader=c3d.Reader(open(c3dFileName,'rb'))
                print('\t'+'\t'+c3dFileName)
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

