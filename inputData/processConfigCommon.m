%Standard Trials collected
index_Static      = 1;
index_Chest       = 2;
index_Conv        = 3; %What does this stand for?
index_Leg         = 4;
index_Side        = 5;
index_Rob         = 6; %Rob?

numberOfDataTrials = 5;
offsetTrialIndex = 1;

index_FeetForcePlate = 1;
index_ChairForcePlate = 2;

trialTypeNames = {'Static','Chest','Conv','Leg','Side','Rob'};
numberOfTrialTypes = length(trialTypeNames);

%Keywords that distinguish C3D data from each trial
%  and the folder that contains the corresponding data from Visual 3D
c3dFileKeyWords = {  'static', 'Che','Con', 'Leg','Sid','Rob'};
folderKeyWords  = {        '', 'Che','Conv','Leg','Sid','Rob'};