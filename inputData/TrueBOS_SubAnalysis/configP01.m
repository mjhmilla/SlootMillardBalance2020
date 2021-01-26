subjectId = 'P01';
inputFolder = [subjectId];
outputFolder=[subjectId];

markerSet = {'FCC','TAM','FM1','FM2','FM5','FAL'};
markerLeftFootPrefix = 'L_';
markerRightFootPrefix = 'R_';

bodyWeight = 61*9.81;


inputC3DFiles = {'Bar_1F_med.c3d',... 
                 'Bar_2F_clo.c3d',...  
                 'Bar_2F_far.c3d',...  
                 'Bar_2F_med.c3d',...  
                 'Run_2F_med.c3d'};


withShoes = zeros(length(inputC3DFiles),1);
withShoes(end,1)=1;                
%0: no shoes
%1: flexible soled shoes - running shoes
%2: stiff soled shoes - hiking shoes                

trialIsValid = ones(length(inputC3DFiles),1);

%trialIsValid(1,1)=0; %'Bar_1F_med.c3d',...  
%trialIsValid(2,1)=0; %'Bar_2F_clo.c3d',...  
%trialIsValid(3,1)=0; %'Bar_2F_far.c3d',...  
%trialIsValid(4,1)=0; %'Bar_2F_med.c3d',...  
%trialIsValid(5,1)=0; %'Run_2F_med.c3d'};




updateFootFrames = zeros(length(inputC3DFiles),3);
% update flag (0/1), index of c3d file to use, index of time sample to use
updateFootFrames(1,:)   = [1,2, (2161-2161+1)]; %Note 2161 is the first data sample: earlier data was chopped
updateFootFrames(end,:) = [1,5,(463-417+1)]; %Note 417 is the first data sample


updateFootScale       = zeros(length(inputC3DFiles),3);
updateFootScale(1,:)  = [1,2,1];

%inputC3DOffsetFile   = 'test0014_done.c3d';
%inputOffsetTimeIndex = 1;

indexLeftForcePlate  = ones(length(inputC3DFiles),1).*2;
indexRightForcePlate = ones(length(inputC3DFiles),1).*1;



processedOutliers = [];