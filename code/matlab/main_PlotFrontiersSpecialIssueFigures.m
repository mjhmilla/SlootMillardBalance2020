
clc;
close all;
clear all;

%%
%
% Constants
%
%%

gravityVec                        = [0;0;-9.81];
flag_EventStart0Reference1End2    = 1;
flag_ConvexHullWithToes0Without1  = 0;

%%
%
% Plot Configuration
%
%%

%Names used here must match the names listed in trialTypeNames in
%FrontiersBalance2020/inputdata/processConfig near line 15
trialsToProcess = {'Chest','Conv','Leg','Side','Rob'};

numberOfFiguresPerPage = 3;
%assert(length(trialsToProcess)==1);
plotConfigFrontiers;
%%
%
% Subject/Grouping/Data
%
%%
subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE05','configE06',...
   'configE07','configE08'};

 
gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;

groups(2) = struct('index',[],'name',[],'color',[]);

indexGroupYoung   = 1;
indexGroupElderly = 2;

groups(indexGroupYoung  ).index = [1,2,3,4,5,6,7,8,9,10];
groups(indexGroupYoung  ).name = 'Y';
groups(indexGroupElderly).index = [11,12,13,14,15, 16,17];
groups(indexGroupElderly).name = 'E';

groups(indexGroupYoung  ).color   = [0,0,0];
groups(indexGroupElderly).color   = gorgeousGreen;




resultsData(length(subjectsToProcess),length(trialsToProcess)) =...
            struct('subjectId','',...
                   'trialId','',...
                   'com2edge',[],...
                   'com2cop' ,[],...
                   'comvel'  ,[],...
                   'fpe2edge',[],...
                   'fpelen'  ,[],...
                   'fpewidth',[]); 
              
              
%The rows store: min, 25%, mean, 75%, max
summaryData(length(subjectsToProcess)) =...
             struct('subjectId','',...
                    'com2edge',zeros(7,6),...
                    'com2cop' ,zeros(7,6),...
                    'comvel'  ,zeros(7,6),...
                    'fpe2edge',zeros(7,6),...
                    'fpelen'  ,zeros(7,6),...
                    'fpewidth',zeros(7,6)); 

%%
%
% Setup the input/output directory structure
%
%%
inputDirRelative    = '../../inputData'         ;
outputDirRelative   = '../../outputData'        ;
frontiersPlotDir    = '../../plots/Frontiers2020Pub/' ;
frontiersDataDir    = '../../outputData/Frontiers2020Pub/';

codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);



figComPlotsVec(length(trialsToProcess))=struct('h',[]);%figure;  
figFpePlotsVec(length(trialsToProcess))=struct('h',[]);%figure;

  axisXGroupWidth = 1+2.5;

  axisYLimComBos     = [-6,17];   
  yLabelOffsetComBos = -0.05;

  axisYLimComCop     = [ 0,16];
  yLabelOffsetComCop = -0.05;
  
  axisYLimComVel     = [0,65];
  yLabelOffsetComVel = -0.05;

  axisYLimFpeBos     = [-0.5,16];
  yLabelOffsetFpeBos = -0.05;

  axisYLimFpeCopWidth    = [-6,6];
  yLabelOffsetFpeCopWidth= -0.05;

  axisYLimFpeCopLength    = [-7,12];
  yLabelOffsetFpeCopLength= -0.05;



flag_firstCall=1;

for indexSubject = 1:1:length(subjectsToProcess)

  if(indexSubject > 1)
    flag_firstCall=0;
  end
    
  indexGroup = 0;
  for k=1:1:length(groups)
    if(any( groups(k).index == indexSubject ))
      indexGroup=k;
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1. Configure the list of input/output files for this subject
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
     
  numberOfTrials = length(inputC3DFiles); 

  disp(['Processing: ', subjectId]);  
  
  subjectIdOriginal = subjectId;
  subjectId = subjectId(2:3);
  if(contains(subjectId,'0'))
    idx = strfind(subjectId,'0');
    if(idx == 1)
      subjectId = subjectId(2);
    end
  end
    
  for indexTrialsToProcess = 1:1:length(trialsToProcess)

    indexTrial = 0;
    for k=1:1:length(trialTypeNames)
      if(contains(trialsToProcess{indexTrialsToProcess},trialTypeNames{k}))
        indexTrial = k;
      end
    end

    if( isempty(figComPlotsVec(indexTrialsToProcess).h)==1)
      figComPlotsVec(indexTrialsToProcess).h=figure;
    end
    if( isempty(figFpePlotsVec(indexTrialsToProcess).h)==1)
      figFpePlotsVec(indexTrialsToProcess).h=figure;
    end
    
    figComPlots = figComPlotsVec(indexTrialsToProcess).h;
    figFpePlots = figFpePlotsVec(indexTrialsToProcess).h;
    
    %Get the trial folder
    trialFolder = outputTrialFolders{indexTrial};    
        
    %Load the processed trial data
    
    c3DFileName       = updateFileExtension(inputC3DFiles{indexTrial},'mat');
    
    if( contains(c3DFileName,'static') == 0 && indexTrial <= 6)
      
      resultsData(indexSubject,indexTrialsToProcess).subjectId ...
         = subjectIdOriginal;
      resultsData(indexSubject,indexTrialsToProcess).trialId ...
         = trialsToProcess{indexTrialsToProcess};

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Fetch all of the data for this trial
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      load([trialFolder,c3DFileName]);
      
      disp(['  :',c3DFileName]);
      % c3dTime
      % c3dMarkers
      % c3dMarkerNames
      % c3dForcePlates
      % c3dForcePlateInfo
      % c3dGrf
      % c3dMarkerUnits
      
      wholeBodyFileName = updateFileExtension(...
                            inputWholeBodyFiles{indexTrial},'mat');                                               
      load([trialFolder,wholeBodyFileName]);
      % wholeBodyData
      % wholeBodyColNames

      anthroFileName    = updateFileExtension(...
                            inputAnthroFiles{indexTrial},'mat');          
      load([trialFolder,anthroFileName]);
      % anthroData
      % anthroColNames
      
      fpeFileName       = outputFpeFileNames{indexTrial};
      load([trialFolder,fpeFileName]);

      capFileName       = outputCapFileNames{indexTrial};  
      load([trialFolder,capFileName]);

      outputFpeToFootHullDistanceFileName ...
        = outputFpeToFootHullDistanceFileNames{indexTrial}   ;
      outputCapToFootHullDistanceFileName ...
        = outputCapToFootHullDistanceFileNames{indexTrial}   ;
      outputComGPToFootHullDistanceFileName ...
        = outputComGPToFootHullDistanceFileNames{indexTrial} ;
      outputCopToFootHullDistanceFileName ...
        = outputCopToFootHullDistanceFileNames{indexTrial}   ;    
      
      
      if( flag_ConvexHullWithToes0Without1==0)
        load([trialFolder, outputFpeToFootHullDistanceFileName ]);
        %fpe2FootConvexHullDist
        load([trialFolder, outputCapToFootHullDistanceFileName]);
        %cap2FootConvexHullDist
        load([trialFolder, outputComGPToFootHullDistanceFileName]);
        %comgp2FootConvexHullDist
        load([trialFolder, outputCopToFootHullDistanceFileName]);       
        %cop2FootConvexHullDist
      else
        load([trialFolder, outputFpeToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat' ]);
        %fpe2FootConvexHullDist
        load([trialFolder, outputCapToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat' ]);
        %cap2FootConvexHullDist
        load([trialFolder, outputComGPToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat' ]);
        %comgp2FootConvexHullDist
        load([trialFolder, outputCopToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat' ]);       
        %cop2FootConvexHullDist        
      end
      
      segmentFileName    = outputSegmentationFileNames{indexTrial};
      load([trialFolder,segmentFileName]);
      
      movementSequenceFileName = outputMovementSequenceFileNames{indexTrial};
      load([trialFolder,movementSequenceFileName]);
      
      indexSittingStatic      = 0;
      indexSittingDynamic     = 1;

      indexCrouchingStable    = 2;
      indexCrouchingUnstable  = 3;

      indexStandingStable     = 4;
      indexStandingUnstable   = 5;        
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Useful columns from the V3D data
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      headerRows          = 4; 
      textInRowBeforeData = 'ITEM';

      colMass     = getColumnIndex({'MASS';'METRIC';'PROCESSED';'X'},...
                                    headerRows,anthroColNames);
      mass        = anthroData(1,colMass);

      colHeight   = getColumnIndex({'HEIGHT';'METRIC';'PROCESSED';'X'},...
                                  headerRows,anthroColNames);
      height      = anthroData(1,colHeight);   
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Whole body quantities
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      colItem    = getColumnIndex({[];[];[];'ITEM'},headerRows,wholeBodyColNames);

      colComPos = zeros(1,3);
      colComPos(1,1)  = getColumnIndex(...
                    {'LBody_CoM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                    headerRows,wholeBodyColNames);
      colComPos(1,2) = colComPos(1,1)+1;
      colComPos(1,3) = colComPos(1,2)+1;

      colComVel = zeros(1,3);
      colComVel(1,1)  = getColumnIndex(...
                    {'LBody_CoM_vel';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                    headerRows,wholeBodyColNames);
      colComVel(1,2) = colComVel(1,1)+1;
      colComVel(1,3) = colComVel(1,2)+1;

      colHo = zeros(1,3);
      colHo(1,1)  = getColumnIndex(...
                    {'LBody_ANGMOM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                    headerRows,wholeBodyColNames);
      colHo(1,2) = colHo(1,1)+1;
      colHo(1,3) = colHo(1,2)+1;

      colJo       = zeros(1,3);
      colJo(1,1)  = getColumnIndex(...
                  {'LBody_MOMINERT';'LINK_MODEL_BASED';'PROCESSED_MATT';'0'},...
                  headerRows,wholeBodyColNames);
      for i=2:1:9
        colJo(1,i) = colJo(1,i-1)+1;
      end      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Movement Sequence
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      movementSequence = sitToStandSequence;
      movementSequenceName = 'Sit2Stand';


      timeLabel     = 'Time (s)';
      forceLabel    = 'Force (N)';
      angleLabel    = 'Angle (deg)';
      distanceLabel = 'Distance (cm)';
      velocityLabel = 'Velocity (cm/s)';
      scaleTime     = 1;
      scaleForce    = 1;
      scaleDistance = 100;
      scaleAngle    = 180/pi;
      scaleVelocity = 100;


      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot com-edge
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
      
      [row,col] = find(subPlotPanelIndex==1);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);
      
      data       = -comgp2FootConvexHullDist.distance;
      scaleData  = scaleDistance;                
      colorData  = groups(indexGroup).color;                
      yLabelData = distanceLabel;
      figTitle   = 'A. CoM to nearest BoS-edge'; 
      axisYLim          = axisYLimComBos;
      axisYLabelOffset  = yLabelOffsetComBos;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));
      
      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);
      
      flag_IntervalMode = 1;
      [figComPlots, dataSummary,dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figComPlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode, 1, flag_firstCall);
      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).com2edge(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events'];
                       
      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).com2edge = metric;

      %resultsData(indexSubject,indexTrialsToProcess)       
%                   'com2edge',metric,...
%                   'com2cop' ,metric,...
%                   'comvel'  ,metric,...
%                   'fpe2edge',metric,...
%                   'fpelen'  ,metric,...
%                   'fpewidth',metric);       
                       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot com-cop
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
      [row,col] = find(subPlotPanelIndex==2);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

      data = sum( (c3dGrf(index_FeetForcePlate).cop(:,1:2) ...
                           - wholeBodyData(:,colComPos(1,1:2))).^2 ,2).^0.5;
      scaleData   = scaleDistance;
      colorData   = groups(indexGroup).color;          
      yLabelData  = distanceLabel;
      figTitle    = 'B. CoM to CoP';
      axisYLim          = axisYLimComCop;
      axisYLabelOffset  = yLabelOffsetComCop;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));

      
      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);
      
      flag_IntervalMode = 1;
      [figComPlots, dataSummary,dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figComPlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode, 0, flag_firstCall);
      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).com2cop(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events']; 
      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).com2cop = metric;                       

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot com velocity
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
      [row,col] = find(subPlotPanelIndex==3);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

      data        = sum( wholeBodyData(:,colComVel).^2, 2).^0.5;
      scaleData   = scaleVelocity;
      colorData   = groups(indexGroup).color;          
      yLabelData  = velocityLabel;
      figTitle    = 'C. CoM Speed';      

      axisYLim          = axisYLimComVel;
      axisYLabelOffset  = yLabelOffsetComVel;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));      

      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);      
      
      flag_IntervalMode = 1;
      [figComPlots, dataSummary, dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figComPlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode,0,flag_firstCall);

      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).comvel(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events'];    

      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).comvel = metric; 
                       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot fpe-to-foot edge
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      [row,col] = find(subPlotPanelIndex==1);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

      data        = -fpe2FootConvexHullDist.distance;
      scaleData   = scaleDistance;
      colorData   = groups(indexGroup).color;          
      yLabelData  = distanceLabel;
      figTitle    = 'A. FPE to nearest BoS-edge';      
      
      axisYLim          = axisYLimFpeBos;
      axisYLabelOffset  = yLabelOffsetFpeBos;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));        

      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);      
      
      flag_IntervalMode = 1;
      [figFpePlots, dataSummary,dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figFpePlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode, 1, flag_firstCall);

      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).fpe2edge(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events'];  
  
      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).fpe2edge = metric;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot fpe-cop variation across width
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      [row,col] = find(subPlotPanelIndex==2);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

      flag_ModeBalancePointsVsCom0VsCop1   = 1;
      flag_ModeAnalyzeBalanceAlong0Across1 = 1;
      [data,dataName] = calcBalancePointDistance(fpeData.r0F0, ...
                    fpeData.u,fpeData.n,...
                    wholeBodyData(:,colComPos),...
                    c3dGrf(index_FeetForcePlate).cop,...
                    flag_ModeBalancePointsVsCom0VsCop1,...
                    flag_ModeAnalyzeBalanceAlong0Across1,...
                    'Fpe');

      scaleData   = scaleDistance;
      colorData   = groups(indexGroup).color;          
      yLabelData  = distanceLabel;
      figTitle    = ['B. ',dataName];      

      axisYLim          = axisYLimFpeCopWidth;
      axisYLabelOffset  = yLabelOffsetFpeCopWidth;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));        

      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);      
      
      flag_IntervalMode = 1;
      [figFpePlots, dataSummary,dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figFpePlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode, 0, flag_firstCall);

      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).fpewidth(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events']; 

      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).fpewidth = metric;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot fpe-cop variation across length
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      [row,col] = find(subPlotPanelIndex==3);          
      subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

      flag_ModeBalancePointsVsCom0VsCop1   = 1;
      flag_ModeAnalyzeBalanceAlong0Across1 = 0; 
      [data,dataName] = calcBalancePointDistance(fpeData.r0F0, ...
                    fpeData.u,fpeData.n,...
                    wholeBodyData(:,colComPos),...
                    c3dGrf(index_FeetForcePlate).cop,...
                    flag_ModeBalancePointsVsCom0VsCop1,...
                    flag_ModeAnalyzeBalanceAlong0Across1,...
                    'Fpe');

      scaleData   = scaleDistance;
      colorData   = groups(indexGroup).color;          
      yLabelData  = distanceLabel;
      figTitle    = ['C. ',dataName];      
      
      axisYLim          = axisYLimFpeCopLength;
      axisYLabelOffset  = yLabelOffsetFpeCopLength;
      yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1)); 

      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);      
      
      flag_IntervalMode = 1;
      [figFpePlots, dataSummary, dataRaw] = ...
        plotBoxWhiskerEventDataFrontiers(...
           figFpePlots, subplotPosition, ...
           movementSequence, ...
           indexSubject, 1,...              
           data, scaleData, ...
           c3dTime,...
           colorData,  ...
           0.33,...
           subjectId, yLabelPos,...
           '', yLabelData, ...
           [figTitle,': ',trialTypeNames{indexTrial}], axisLim,...
           flag_IntervalMode, 0, flag_firstCall);
      axisLim = zeros(1,4);
      axisLim(1)=0.5;
      axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
      axisLim(3)=axisYLim(1);
      axisLim(4)=axisYLim(2);         
      axis(axisLim);
      set(gca,'XTickLabel',[]); 
      
      summaryData(indexSubject).fpelen(:,indexTrial) = ...
                       [ dataSummary.min;...
                         dataSummary.p25;...
                         dataSummary.mean;...
                         dataSummary.p75;...
                         dataSummary.max;...
                         dataSummary.events'];

      metric = struct('min'      , dataSummary.min,...
                      'mean'     , dataSummary.mean,...
                      'max'      , dataSummary.max,...
                      'percent25', dataSummary.p25,...
                      'percent75', dataSummary.p75,...
                      'events', dataSummary.events,...
                      'raw',dataRaw);
                    
      resultsData(indexSubject,indexTrialsToProcess).fpelen = metric;

    end
      
      

    
  end
  
end

save([frontiersDataDir,'subjectTrialMetricData.mat'],'resultsData');

%%
%
% Add the group plot information
%
%%

for indexTrialsToProcess = 1:1:length(trialsToProcess)

  indexTrial = 0;
  for k=1:1:length(trialTypeNames)
    if(contains(trialsToProcess{indexTrialsToProcess},trialTypeNames{k}))
      indexTrial = k;
    end
  end
  
  figComPlots = figComPlotsVec(indexTrialsToProcess).h;
  figFpePlots = figFpePlotsVec(indexTrialsToProcess).h;
  
  if( contains(c3DFileName,'static') == 0 && indexTrial <= 6)

    for indexGroup = 1:1:length(groups)
      groupMembers = groups(indexGroup).index;

      xPositionLine = length(summaryData) + 1;
      xPositionData = length(summaryData) + 1 + indexGroup ;
      flag_drawLine = 0;
      if(indexGroup == 1)
        flag_drawLine = 1;
      end
      %%
      %Com-Bos
      %%
      [row,col] = find(subPlotPanelIndex==1);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);

      figComPlots = plotGroupAntsOnALog(figComPlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'com2edge', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetComBos,flag_drawLine);

      %%
      %Com-Cop
      %%
      [row,col] = find(subPlotPanelIndex==2);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);
      
      figComPlots = plotGroupAntsOnALog(figComPlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'com2cop', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetComCop,flag_drawLine); 
      %%
      %Com-Vel
      %%
      [row,col] = find(subPlotPanelIndex==3);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);
      
      figComPlots = plotGroupAntsOnALog(figComPlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'comvel', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetComVel,flag_drawLine); 

      %%
      %Fpe-Edge
      %%
      [row,col] = find(subPlotPanelIndex==1);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);
      
      figFpePlots = plotGroupAntsOnALog(figFpePlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'fpe2edge', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetFpeBos,flag_drawLine); 

      %%
      %Fpe-Width
      %%
      [row,col] = find(subPlotPanelIndex==2);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);
      
      figFpePlots = plotGroupAntsOnALog(figFpePlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'fpewidth', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetFpeCopWidth,flag_drawLine); 

      %%
      %Fpe-Length
      %%
      [row,col] = find(subPlotPanelIndex==3);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);
      
      figFpePlots = plotGroupAntsOnALog(figFpePlots, subPlotVec,...
                                      groups(indexGroup),...
                                      summaryData,...
                                      'fpelen', indexTrial,...
                                      xPositionData, xPositionLine,...
                                      yLabelOffsetFpeCopLength,flag_drawLine); 

    end

  end
end



for indexTrialsToProcess=1:1:length(trialsToProcess)

    figComPlots = figComPlotsVec(indexTrialsToProcess).h;
    figFpePlots = figFpePlotsVec(indexTrialsToProcess).h;
  
    figure(figComPlots);
    configPlotExporter;
    print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                    trialsToProcess{indexTrialsToProcess},'CoM','.pdf']);

    figure(figFpePlots);
    configPlotExporter;
    print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                    trialsToProcess{indexTrialsToProcess},'FPE','.pdf']);

end
