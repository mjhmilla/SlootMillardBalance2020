
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

numberOfPhases = 2;
indexPhaseStart2SeatOff=1;
indexPhaseSeatOff2Stand=2;
phaseNames = {'Phase1Sit2SeatOff','Phase2SeatOff2Stand'};
% 1st phase: start to seatoff
% 2nd phase: seatoff to standing;

%%
%
% Plot Configuration
%
%%

%Names used here must match the names listed in trialTypeNames in
%FrontiersBalance2020/inputdata/processConfig near line 15
trialsToProcess = {'Chest','Conv','Leg','Side','Rob'};

numberOfFiguresPerPage = 3;
flag_printGroupSummaryStats = 1;
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
   'configE07','configE08','configE09'};

 
gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;

groups(2) = struct('index',[],'name',[],'color',[],'label',[]);

indexGroupYoung   = 1;
indexGroupElderly = 2;

groups(indexGroupYoung  ).index = [1,2,3,4,5,6,7,8,9,10];
groups(indexGroupYoung  ).name = 'Y';
groups(indexGroupYoung  ).label = 'Young (Y)';
groups(indexGroupElderly).index = [11,12,13,14,15, 16,17,18];
groups(indexGroupElderly).name = 'E';
groups(indexGroupElderly).label = 'Elderly (E)';

groups(indexGroupYoung  ).color   = [0,0,0];
groups(indexGroupElderly).color   = gorgeousGreen;




resultsData(length(subjectsToProcess),length(trialsToProcess),numberOfPhases) =...
            struct('subjectId','',...
                   'trialId','',...
                   'com2edge',[],...
                   'com2cop' ,[],...
                   'comvel'  ,[],...
                   'fpe2edge',[],...
                   'fpelen'  ,[],...
                   'fpewidth',[],...
                   'duration',[]); 
              
              
%The rows store: min, 25%, mean, 75%, max
summaryData(length(subjectsToProcess),numberOfPhases) =...
             struct('subjectId','',...
                    'com2edge',zeros(8,6).*NaN,...
                    'com2cop' ,zeros(8,6).*NaN,...
                    'comvel'  ,zeros(8,6).*NaN,...
                    'fpe2edge',zeros(8,6).*NaN,...
                    'fpelen'  ,zeros(8,6).*NaN,...
                    'fpewidth',zeros(8,6).*NaN); 

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

  axisYLimComCop     = [ 0,13];
  yLabelOffsetComCop = -0.05;
  
  axisYLimComVel     = [0,60];
  yLabelOffsetComVel = -0.05;

  axisYLimFpeBos     = [-0.5,16];
  yLabelOffsetFpeBos = -0.05;

  axisYLimFpeCopWidth    = [-3,4.5];
  yLabelOffsetFpeCopWidth= -0.05;

  axisYLimFpeCopLength    = [-3,9];
  yLabelOffsetFpeCopLength= -0.05;



flag_firstCall=1;

for indexSubject = 1:1:length(subjectsToProcess)

  for indexPhase=1:1:numberOfPhases
    summaryData(indexSubject,indexPhase).com2edge = zeros(8,6).*NaN;
    summaryData(indexSubject,indexPhase).com2cop  = zeros(8,6).*NaN;
    summaryData(indexSubject,indexPhase).comvel   = zeros(8,6).*NaN;
    summaryData(indexSubject,indexPhase).fpe2edge = zeros(8,6).*NaN;
    summaryData(indexSubject,indexPhase).fpelen   = zeros(8,6).*NaN;
    summaryData(indexSubject,indexPhase).fpewidth = zeros(8,6).*NaN;  
  end
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
    if(isempty(inputC3DFiles{indexTrial})==0)

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

        for indexPhase=1:1:numberOfPhases
          resultsData(indexSubject,indexTrialsToProcess,indexPhase).subjectId ...
             = subjectIdOriginal;
          resultsData(indexSubject,indexTrialsToProcess,indexPhase).trialId ...
             = trialsToProcess{indexTrialsToProcess};

        end
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
        figTitle   = 'A. COM to nearest BOS-edge'; 
        axisYLim          = axisYLimComBos;
        axisYLabelOffset  = yLabelOffsetComBos;
        yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));

        axisLim = zeros(1,4);
        axisLim(1)=0.5;
        axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
        axisLim(3)=axisYLim(1);
        axisLim(4)=axisYLim(2);

    
        

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
             figTitle, axisLim,...
              1, flag_firstCall,numberOfPhases);
        axis(axisLim);
        set(gca,'XTickLabel',[]); 

        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).com2edge(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p,:));

          resultsData(indexSubject,indexTrialsToProcess,p).com2edge = metric;
          resultsData(indexSubject,indexTrialsToProcess,p).duration =...
            dataSummary(p).duration;
        end

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
        figTitle    = 'B. COM to COP';
        axisYLim          = axisYLimComCop;
        axisYLabelOffset  = yLabelOffsetComCop;
        yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));


        axisLim = zeros(1,4);
        axisLim(1)=0.5;
        axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
        axisLim(3)=axisYLim(1);
        axisLim(4)=axisYLim(2);


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
             figTitle, axisLim,...
              0, flag_firstCall,numberOfPhases);
        axis(axisLim);
        set(gca,'XTickLabel',[]); 

        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).com2cop(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p));

          resultsData(indexSubject,indexTrialsToProcess,p).com2cop = metric;

        end      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot com velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        [row,col] = find(subPlotPanelIndex==3);          
        subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

        data        = sum( wholeBodyData(:,colComVel).^2, 2).^0.5;
        scaleData   = scaleVelocity;
        colorData   = groups(indexGroup).color;          
        yLabelData  = velocityLabel;
        figTitle    = 'C. COM Speed';      

        axisYLim          = axisYLimComVel;
        axisYLabelOffset  = yLabelOffsetComVel;
        yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));      

        axisLim = zeros(1,4);
        axisLim(1)=0.5;
        axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
        axisLim(3)=axisYLim(1);
        axisLim(4)=axisYLim(2);      


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
             figTitle, axisLim,...
             0,flag_firstCall,numberOfPhases);

        axis(axisLim);
        set(gca,'XTickLabel',[]); 


        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).comvel(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p));

          resultsData(indexSubject,indexTrialsToProcess,p).comvel = metric;

        end      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot fpe-to-foot edge
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        [row,col] = find(subPlotPanelIndex==1);          
        subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

        data        = -fpe2FootConvexHullDist.distance;
        scaleData   = scaleDistance;
        colorData   = groups(indexGroup).color;          
        yLabelData  = distanceLabel;
        figTitle    = 'A. FPE to nearest BOS-edge';      

        axisYLim          = axisYLimFpeBos;
        axisYLabelOffset  = yLabelOffsetFpeBos;
        yLabelPos = axisYLim(1)+axisYLabelOffset*(axisYLim(2)-axisYLim(1));        

        axisLim = zeros(1,4);
        axisLim(1)=0.5;
        axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
        axisLim(3)=axisYLim(1);
        axisLim(4)=axisYLim(2);      


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
             figTitle, axisLim,...
              1, flag_firstCall,numberOfPhases);

        axis(axisLim);
        set(gca,'XTickLabel',[]); 


        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).fpe2edge(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p));

          resultsData(indexSubject,indexTrialsToProcess,p).fpe2edge = metric;

        end                     


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
                      'FPE');

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
             figTitle, axisLim,...
              0, flag_firstCall,numberOfPhases);

        axis(axisLim);
        set(gca,'XTickLabel',[]); 

        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).fpewidth(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p));

          resultsData(indexSubject,indexTrialsToProcess,p).fpewidth = metric;

        end           

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
                      'FPE');

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
             figTitle, axisLim,...
              0, flag_firstCall,numberOfPhases);
        axisLim = zeros(1,4);
        axisLim(1)=0.5;
        axisLim(2)=length(subjectsToProcess)+axisXGroupWidth;
        axisLim(3)=axisYLim(1);
        axisLim(4)=axisYLim(2);         
        axis(axisLim);
        set(gca,'XTickLabel',[]); 

        for p=1:1:numberOfPhases
          summaryData(indexSubject,p).fpelen(:,indexTrial) = ...
                           [ dataSummary(p).min;...
                             dataSummary(p).p25;...
                             dataSummary(p).mean;...
                             dataSummary(p).p75;...
                             dataSummary(p).max;...
                             dataSummary(p).start;...
                             dataSummary(p).end;...
                             dataSummary(p).duration];

          metric = struct('min'      , dataSummary(p).min,...
                          'mean'     , dataSummary(p).mean,...
                          'max'      , dataSummary(p).max,...
                          'percent25', dataSummary(p).p25,...
                          'percent75', dataSummary(p).p75,...
                          'meanValuePhaseStart', dataSummary(p).start,...
                          'meanValuePhaseEnd', dataSummary(p).end,...
                          'raw',dataRaw(p));

          resultsData(indexSubject,indexTrialsToProcess,p).fpelen = metric;

        end       

      end



    end
  end
  
end

save([frontiersDataDir,'subjectTrialMetricData.mat'],'resultsData');

%%
%
% Add the group plot information
%
%%
metricNameList = {'com2edge','com2cop','comvel',...
                  'fpe2edge','fpewidth','fpelen'};
metricSubPlotPanelIndexList = [1,2,3,1,2,3];
metricPlotIndexList = [1,1,1,2,2,2];
            
groupSummaryStats(length(trialsToProcess),length(groups),length(metricNameList),numberOfPhases)...
  = struct('trialName','','groupName','','metricName','','data',zeros(8,1));

emptyColumnNames = cell(1,length(groups)*length(metricNameList));
emptyRowNames    = cell(1,8);
emptyTrialData = zeros(8,length(groups)*length(metricNameList));


groupSummaryTables(length(trialsToProcess),numberOfPhases) =...
  struct( 'data',             emptyTrialData,...
          'rowLabels',        {emptyRowNames},...
          'columnLabelMetric',{emptyColumnNames},...
          'columnLabelGroup', {emptyColumnNames},...          
          'trialType',        ''               );

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

  disp(['Processing: ',trialsToProcess{indexTrialsToProcess}]);    
    
    for indexPhase=1:1:numberOfPhases
      groupSummaryTables(indexTrialsToProcess,indexPhase).rowLabels ...
        = {'mean(min)','25','mean','75','mean(max)',...
           'mean(seat-off)','mean(standing)','mean(duration)'};
      groupSummaryTables(indexTrialsToProcess,indexPhase).trialType ...
        = trialTypeNames{indexTrial};
    end  
  
    for indexGroup = 1:1:length(groups)
      
      
      disp(['  Group: ',groups(indexGroup).name]);
      groupMembers = groups(indexGroup).index;

      xPositionLine = size(summaryData,1) + 1;
      xPositionData = size(summaryData,1) + 1 + indexGroup ;
      flag_drawLine = 0;
      if(indexGroup == 1)
        flag_drawLine = 1;
      end
      
      for indexPhase=1:1:numberOfPhases
        for indexMetric = 1:1:length(metricNameList)
          indexTableColumn = (indexMetric-1)*length(groups) + indexGroup;

          groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelMetric{1,indexTableColumn} ...
            = metricNameList{indexMetric};
          groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelGroup{1,indexTableColumn} ...
            = groups(indexGroup).name;


          metricName              = metricNameList{indexMetric};
          metricPlotIndex         = metricPlotIndexList(1,indexMetric);
          metricSubPlotPanelIndex = metricSubPlotPanelIndexList(1,indexMetric);
          figH = [];
          switch(metricPlotIndex)
            case 1
              figH = figComPlots;
            case 2
              figH = figFpePlots;
            otherwise assert(0);
          end

          [row,col] = find(subPlotPanelIndex==metricSubPlotPanelIndex);          
          subPlotVec = reshape(subPlotPanel(row,col,:),1,4);

          groupSummaryStats(indexTrialsToProcess,indexGroup,indexMetric).trialName ...
            = trialTypeNames{indexTrial};
          groupSummaryStats(indexTrialsToProcess,indexGroup,indexMetric).groupName ...
            = groups(indexGroup).name;
          groupSummaryStats(indexTrialsToProcess,indexGroup,indexMetric).metricName ...
            = metricName;


          groupSummaryData = averageAcrossSubjects(groups(indexGroup),...
                                                summaryData(:,indexPhase),...
                                                metricName,...
                                                indexTrial,...
                                                flag_printGroupSummaryStats);      
          groupSummaryStats(indexTrialsToProcess,indexGroup,indexMetric,indexPhase).data = ...
            groupSummaryData;

          groupSummaryTables(indexTrialsToProcess,indexPhase).data(:,indexTableColumn) =...
            groupSummaryData;        
          if(indexPhase==indexPhaseSeatOff2Stand)
            figH = plotGroupAntsOnALog(figH, subPlotVec,...
                                            groups(indexGroup),...
                                            groupSummaryData,...
                                            xPositionData, xPositionLine,...
                                            yLabelOffsetComBos,...
                                            flag_drawLine);        
          end

        end
      end

    end

  end
end



for indexTrialsToProcess=1:1:length(trialsToProcess)
  for indexPhase=1:1:numberOfPhases
  
    tableName = ['../../outputData/Frontiers2020Pub/table',...
      groupSummaryTables(indexTrialsToProcess).trialType,...
      phaseNames{indexPhase},'Group.csv'];
  
% groupSummaryTables(length(trialTypeNames)) = ...
%   struct('data',zeros(5, length(groups)*length(metricList)),...
%          'rowLabels',{'max','75','mean','25','min'},...
%          'columnLabelMetric',cell(1, length(groups)*length(metricList)),...
%          'columnLabelGroup', cell(1, length(groups)*length(metricList)));

    
    fid =fopen(tableName,'w');
      fprintf(fid,',%s',groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelMetric{1});    
      for i=2:1:length(groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelMetric)
        fprintf(fid,',%s',groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelMetric{i});
      end
      fprintf(fid,'\n');
      fprintf(fid,',%s',groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelGroup{1});    
      for i=2:1:length(groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelGroup)
        fprintf(fid,',%s',groupSummaryTables(indexTrialsToProcess,indexPhase).columnLabelGroup{i});
      end
      fprintf(fid,'\n');
    
      for i=1:1:length(groupSummaryTables(indexTrialsToProcess,indexPhase).rowLabels)
        fprintf(fid,'%s',groupSummaryTables(indexTrialsToProcess,indexPhase).rowLabels{i});
        for j=1:1:size(groupSummaryTables(indexTrialsToProcess,indexPhase).data,2)
          fprintf(fid,',%1.2f',groupSummaryTables(indexTrialsToProcess,indexPhase).data(i,j));          
        end
        fprintf(fid,'\n');       
      end      
    fclose(fid);
  end  
  
    figComPlots = figComPlotsVec(indexTrialsToProcess).h;
    figFpePlots = figFpePlotsVec(indexTrialsToProcess).h;
    
    legendFontSize = get(groot,'defaultAxesFontSize');
    
    figure(figComPlots);
    [row,col] = find(subPlotPanelIndex==1);          
    subplotPosition = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subplotPosition);
    axisLim = axis;
    yMin = axisLim(3);
    yMax = axisLim(4);
    x0 = 12;
    y0 = yMin + 0.1*(yMax-yMin);
    dy = 0.1*(yMax-yMin);
    dx = 1;    
    for indexGroup=1:1:length(groups)    
      plot([x0,x0+dx],[y0,y0],'Color',groups(indexGroup).color);
      hold on;
      text(x0+2*dx,y0,groups(indexGroup).label,'FontSize',legendFontSize,...
           'Interpreter','latex','HorizontalAlignment','left');
      hold on;
      y0 = y0+dy;
    end
    
    configPlotExporter;
    print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                    trialsToProcess{indexTrialsToProcess},'CoM','.pdf']);

    figure(figFpePlots);
    [row,col] = find(subPlotPanelIndex==1);          
    subplotPosition = reshape(subPlotPanel(row,col,:),1,4);    
    subplot('Position',subplotPosition);

    axisLim = axis;
    yMin = axisLim(3);
    yMax = axisLim(4);
    x0 = 12;
    y0 = yMin + 0.1*(yMax-yMin);
    dy = 0.1*(yMax-yMin);
    dx = 1;    
    for indexGroup=1:1:length(groups)    
      plot([x0,x0+dx],[y0,y0],'Color',groups(indexGroup).color);
      hold on;
      text(x0+2*dx,y0,groups(indexGroup).label,'FontSize',legendFontSize,...
           'Interpreter','latex','HorizontalAlignment','left');
      hold on;
      y0 = y0+dy;      
    end
    
    configPlotExporter;
    print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                    trialsToProcess{indexTrialsToProcess},'FPE','.pdf']);

end
