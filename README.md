# Description

author: M.Millard
date  : 20/1/2020 - 11/5/2022
email : matthew.millard@inspo.uni-stuttgart.de


This repository contains the collection of scripts (without data) used to process the sit-to-stand dynamic balance of younger and older adults as described in Sloot LH, Millard M, Werner C, Mombaur K. Slow but Steady: Similar Sit-to-Stand balance at Seat-Off in older vs. younger adults. Frontiers in sports and active living. 2020:144.

# Overview

This is intended as a code example to help other researchers apply technical measures to analyze the dynamic balance of people. Here the balance measures considered are: the foot-placement estimator, the capture-point, and the functional base of support as described in Sloot et al. In this repostory you will find the following pipeline to analyze dynamic balance data:

1. Kinematic data: c3d files and Visual3D data. The Visual3D data includes: the center-of-mass (COM) position, velocity, linear momentum, angular momentum, and the whole body moment of inertia about the COM. This data (not included in the repository) is stored in /inputData and loaded into a generic structure by making use of the configuration (e.g. configE01.m) files and the configuration scripts (procesConfig.m and processConfigCommon.m)

2. The collection of scripts in code/matlab/ performs 4 functions:

- main\_ProcessSubjectFunctionalBOSAnalysis.m : evalutes the generic scalable base-of-support model.
- main\_ProcessBalanceData.m : uses the kinematic data to evaluate the FPE, CAP, and functional base of support
- main\_StatisticalAnalysisPreprocessing.m : extracts measures at specific times of interest and packages these measures in a structure for later statistical analysis
- main\_PlotFrontiersSpecialIssueFiguresR2.m : generates the figures that appear in Sloot et al.

Any file in code/matlab/ that begins with main\_ is designed to be run manually. Not all main\_ files are fully independent since the data in Sloot et al. is processed in the following stages: main\_ProcessSubjectFunctionalBOSAnalysis.m, main\_ProcessBalanceData.m, main\_StatisticalAnalysisPreprocessing.m, and main\_PlotFrontiersSpecialIssueFiguresR2.m.

3. To get started in applying these technical analysis to your own work, I recommend starting with the main function that does much of the number crunching used in Sloot et al.:

- main\_ProcessBalanceData.m
  - Preprocessing
    - getC3DTrialData (line 282) 
    - getAnthropometryData (line 298) 
    - getWholeBodyTrialData (line 318)
  - Technical analysis 
    - constructFootFrames (line 411)
    - process3DFootPlacementEstimator (line 446)
    - processCapturePoint (line 470)
  - Metric evaluation
    - processDistanceToConvexHull (line 542)
    - processDistanceToConvexHullUsingNormBosModel (line 548)
  - Movement segmentation
    - segmentDataKMeans (line 664)
    - extractSitToStandSequence (line 702)
  - Generation of RBDL model (requires a special version of Manish Sreenivasa's ModelFactory)
    - writeModelFactoryModelFile (line 878)
    - writeModelFactoryEnvironmentFile (line 897)
    - createModel (line 908)

Each of these functions appear in code/matlab/ should be read with the Sloot et al. paper at hand. Some functions, particularly the technical analysis functions, are quite detailed and will take patience and time to follow.

