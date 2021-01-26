subjectId = 'P06';
inputFolder = [subjectId];
outputFolder=[subjectId];

markerSet = {'FCC','TAM','FM1','FM2','FM5','FAL'};
markerLeftFootPrefix = 'L_';
markerRightFootPrefix = 'R_';

bodyWeight = 74*9.81;




inputC3DFiles = { 'Bar_1F_med.c3d',...  
                  'Bar_2F_far.c3d',...  
                  'Hik_1F_med.c3d',...  
                  'Run_1F_med.c3d',...
                  'Bar_2F_clo.c3d',...  
                  'Bar_2F_med.c3d',...  
                  'Hik_2F_med.c3d',...  
                  'Run_2F_med.c3d'};


withShoes = zeros(length(inputC3DFiles),1);
withShoes(3,1) = 2;
withShoes(4,1) = 1;

withShoes(7,1) = 2;
withShoes(8,1) = 1;

%0: no shoes
%1: flexible soled shoes - running shoes
%2: stiff soled shoes - hiking shoes      

trialIsValid = ones(length(inputC3DFiles),1);

%trialIsValid(1,1)=0; %'Bar_1F_med.c3d',...  
%trialIsValid(2,1)=0; %'Bar_2F_far.c3d',...  
%trialIsValid(3,1)=0; %'Hik_1F_med.c3d',...  
%trialIsValid(4,1)=0; %'Run_1F_med.c3d',...

%trialIsValid(5,1)=0; %'Bar_2F_clo.c3d',...  
trialIsValid(6,1)=0; %'Bar_2F_med.c3d',...  
%trialIsValid(7,1)=0; %'Hik_2F_med.c3d',...  
%trialIsValid(8,1)=0; %'Run_2F_med.c3d'};


updateFootFrames = zeros(length(inputC3DFiles),3);
% update flag (0/1), index of c3d file to use, index of time sample to use

bareFootFrame     = [1,6,1];
hikingShoeFrame   = [1,7,1];
runningShoeFrame  = [1,8,1];


updateFootFrames(1,:)  = bareFootFrame; 
updateFootFrames(2,:)  = bareFootFrame; 
updateFootFrames(3,:)  = hikingShoeFrame; 
updateFootFrames(4,:)  = runningShoeFrame; 
updateFootFrames(5,:)  = bareFootFrame; 
updateFootFrames(6,:)  = bareFootFrame; 
updateFootFrames(7,:)  = hikingShoeFrame; 
updateFootFrames(8,:)  = runningShoeFrame; 


oneFoot = [1,3,4];

indexLeftForcePlate  = ones(length(inputC3DFiles),1).*2;
indexRightForcePlate = ones(length(inputC3DFiles),1).*1;

updateFootScale = zeros(length(inputC3DFiles),3);
updateFootScale(1,:) = bareFootFrame;


%inputC3DOffsetFile   = 'test0014_done.c3d';
%inputOffsetTimeIndex = 1;



processedOutliers = [];