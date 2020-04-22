
clc;
close all;
clear all;


flag_loadPrecomputedSubjectStatistics = 0;
flag_loadPrecomputedGroupStatistics   = 0;
flag_writeCSVGroupTables              = 1;
%%
%
% Constants
%
%%
flag_ConvexHullWithToes0Without1  = 0;
nameToeTag = '';
if(flag_ConvexHullWithToes0Without1 == 1)
  nameToeTag = '_NoToes';
end


metricNameList = {'duration',...
                  'com2edge','com2cop','comvel',...
                  'comvelx','comvely','comvelz',...
                  'angvel','angvelx','angvely','angvelz','angveln',...
                  'fpe2edge','fpewidth','fpelen','fpeasmhkz','timemaxfpe'};
                
% duration: the time taken to complete the phase
%
% com2edge: distance between the com-ground-projection (comgp) and the 
%           nearest edge of the base of support (bos). Positive means the  
%           comgp is inside the bos, negative outside.
%
% com2cop:  the distance between the com-ground-projection and the
%           center-of-pressure.
%
% comvel: the velocity of the com
%
% fpe2edge: distance between the foot-placement-estimator and the 
%           nearest edge of the base of support (bos). Positive means the  
%           comgp is inside the bos, negative outside.
%
% fpewidth: distance between the fpe and the cop in the n direction
%
% fpelength: distance between the fpe and the cop in the u direction
%
% fpeasmhkz: fpe assumption that the vertical component of the momentum
%            vector computed about the comgp point is small. What is
%            stored in this component is the percentage of Hgp that
%            is in the vertical direction.

metricSubFields =  {'phase','start','end'};

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

%%
%
% Subject/Grouping/Data
%
%%
subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE04','configE05',...
   'configE06','configE07','configE08','configE09'};

 

groups(2) = struct('index',[],'name',[],'color',[],'label',[]);

gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;

indexGroupYoung   = 1;
indexGroupElderly = 2;

groups(indexGroupYoung  ).index = [1,2,3,4,5,6,7,8,9,10];
groups(indexGroupYoung  ).name = 'Y';
groups(indexGroupYoung  ).label = 'Young (Y)';
groups(indexGroupElderly).index = [11,12,13,14,15, 16,17,18,19];
groups(indexGroupElderly).name = 'E';
groups(indexGroupElderly).label = 'Elderly (E)';

groups(indexGroupYoung  ).color   = [0,0,0];
groups(indexGroupElderly).color   = gorgeousGreen;

assert( (length(groups(indexGroupElderly).index)...
        +length(groups(indexGroupYoung  ).index)) ...
        == length(subjectsToProcess)); 


numberOfInterpolationPoints = 20; %Number of points to use to interpolate the raw data


basicStatsStruct=struct('p05',[],'p25',[],'median',[],'p75',[],'p95',[],...
                        'min',[],'mean',[],'max',[],...
                        'n',0,'data',[]);      
                        
phaseStatsStruct = struct('phase',basicStatsStruct,...
                          'start',basicStatsStruct,...
                          'end',basicStatsStruct);
                        
                        
subjectData(length(subjectsToProcess),length(trialsToProcess),numberOfPhases) =...
            struct('subjectId','',...
                   'trialId','',...
                   'duration'       ,basicStatsStruct,...
                   'com2edge'      ,phaseStatsStruct,...
                   'com2cop'       ,phaseStatsStruct,...
                   'comvel'        ,phaseStatsStruct,...
                   'comvelx'       ,phaseStatsStruct,...
                   'comvely'       ,phaseStatsStruct,...
                   'comvelz'       ,phaseStatsStruct,...
                   'angvel'        ,phaseStatsStruct,...                   
                   'angvelx'       ,phaseStatsStruct,...
                   'angvely'       ,phaseStatsStruct,...
                   'angvelz'       ,phaseStatsStruct,...  
                   'angveln'       ,phaseStatsStruct,...
                   'fpe2edge'      ,phaseStatsStruct,...
                   'fpelen'        ,phaseStatsStruct,...
                   'fpewidth'      ,phaseStatsStruct,...
                   'fpeasmhkz'     ,phaseStatsStruct,...
                   'timemaxfpe'    ,phaseStatsStruct); 
              
groupData(length(groups),length(trialsToProcess),numberOfPhases) = ...
      struct(  'name','',...
               'label','',...
               'subjectId',{''},...
               'subjectIndex',[0],...
               'trialId','',...
               'trialTypeIndex',0,...
               'duration'      ,phaseStatsStruct,...
               'com2edge'      ,phaseStatsStruct,...
               'com2cop'       ,phaseStatsStruct,...
               'comvel'        ,phaseStatsStruct,...
               'comvelx'       ,phaseStatsStruct,...
               'comvely'       ,phaseStatsStruct,...
               'comvelz'       ,phaseStatsStruct,...
               'angvel'        ,phaseStatsStruct,...                   
               'angvelx'       ,phaseStatsStruct,...
               'angvely'       ,phaseStatsStruct,...
               'angvelz'       ,phaseStatsStruct,...
               'angveln'       ,phaseStatsStruct,...               
               'fpe2edge'      ,phaseStatsStruct,...
               'fpelen'        ,phaseStatsStruct,...
               'fpewidth'      ,phaseStatsStruct,...
               'fpeasmhkz'     ,phaseStatsStruct,...
               'timemaxfpe'    ,phaseStatsStruct);

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
  processConfigCommon;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);




if(flag_loadPrecomputedSubjectStatistics == 0)

  for indexSubject = 1:1:length(subjectsToProcess)


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

        %Get the trial folder
        trialFolder = outputTrialFolders{indexTrial};    

        %Load the processed trial data
        c3DFileName       = updateFileExtension(inputC3DFiles{indexTrial},'mat');

        if( contains(c3DFileName,'static') == 0 && indexTrial <= 6)

          if(indexSubject==18 && indexTrial==5)
            here=1;
          end
          
          for indexPhase=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,indexPhase).subjectId ...
               = subjectIdOriginal;
            subjectData(indexSubject,indexTrialsToProcess,indexPhase).trialId ...
               = trialsToProcess{indexTrialsToProcess};
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fetch all of the data for this trial
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          c3dGrfDataAvailable = -1;
          
          load([trialFolder,c3DFileName]);
          
          %Make sure that this flag has been overwritten: c3dGrfDataAvailable
          assert(c3dGrfDataAvailable ~= -1);
         
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
          angularVelocityLabel = 'Angular Velocity (deg/s)';
          scaleTime     = 1;
          scaleForce    = 1;
          scaleDistance = 100;
          scaleAngle    = 180/pi;
          scaleVelocity = 100;
          scaleAngularVelocity = scaleAngle;
          


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com-edge
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

          data       = -comgp2FootConvexHullDist.distance;
          scaleData  = scaleDistance;                

          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).com2edge = phaseStats(p);
            subjectData(indexSubject,indexTrialsToProcess,p).duration = timingStats(p);          
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com-cop
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          
          if(c3dGrfDataAvailable==1)
            data = sum( (c3dGrf(index_FeetForcePlate).cop(:,1:2) ...
                                 - wholeBodyData(:,colComPos(1,1:2))).^2 ,2).^0.5;
            scaleData   = scaleDistance;
            [phaseStats, timingStats] = ...
              calcMovementSequenceDataSummary(movementSequence, ...
                  c3dTime, data, scaleData,numberOfPhases,...
                  numberOfInterpolationPoints); 

            for p=1:1:numberOfPhases
              subjectData(indexSubject,indexTrialsToProcess,p).com2cop = phaseStats(p);    
            end
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com velocity
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = sum( wholeBodyData(:,colComVel).^2, 2).^0.5;
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).comvel = phaseStats(p);    
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com velocity x
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = wholeBodyData(:,colComVel(1,1));
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).comvelx = phaseStats(p);    
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com velocity y
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = wholeBodyData(:,colComVel(1,2));
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).comvely = phaseStats(p);    
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com velocity z
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = wholeBodyData(:,colComVel(1,3));
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).comvelz = phaseStats(p);    
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com angular velocity
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = sum( fpeData.w0C0.^2, 2).^0.5;
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).angvel = phaseStats(p);    
          end          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com angular velocity x
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = fpeData.w0C0(:,1);
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).angvelx = phaseStats(p);    
          end          

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com angular velocity y
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = fpeData.w0C0(:,2);
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).angvely = phaseStats(p);    
          end          

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com angular velocity z
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = fpeData.w0C0(:,3);
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).angvelz = phaseStats(p);    
          end                    
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Com angular velocity n
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          data        = fpeData.w0C0n;
          scaleData   = scaleVelocity;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).angveln = phaseStats(p);    
          end                    
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe-to-foot edge
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          data        = -fpe2FootConvexHullDist.distance;
          scaleData   = scaleDistance;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, data, scaleData,numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).fpe2edge = phaseStats(p);    
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe-cop variation across width
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          
          if(c3dGrfDataAvailable ==1)
            flag_ModeBalancePointsVsCom0VsCop1 =   1;
            flag_ModeAnalyzeBalanceAlong0Across1 = 1;
            [data,dataName] = calcBalancePointDistance(fpeData.r0F0, ...
                          fpeData.u,fpeData.n,...
                          wholeBodyData(:,colComPos),...
                          c3dGrf(index_FeetForcePlate).cop,...
                          flag_ModeBalancePointsVsCom0VsCop1,...
                          flag_ModeAnalyzeBalanceAlong0Across1,...
                          'Fpe');
            scaleData   = scaleDistance;
            [phaseStats, timingStats] = ...
              calcMovementSequenceDataSummary(movementSequence, ...
                  c3dTime, data, scaleData,numberOfPhases,...
                  numberOfInterpolationPoints); 

            for p=1:1:numberOfPhases
              subjectData(indexSubject,indexTrialsToProcess,p).fpewidth = phaseStats(p);    
            end
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe-cop variation across length
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

          if(c3dGrfDataAvailable==1)
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
           [phaseStats, timingStats] = ...
              calcMovementSequenceDataSummary(movementSequence, ...
                  c3dTime, data, scaleData,numberOfPhases,...
                  numberOfInterpolationPoints); 

            for p=1:1:numberOfPhases
              subjectData(indexSubject,indexTrialsToProcess,p).fpelen = phaseStats(p);    
            end
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe assumption that the vertical component of Hgp is negligible
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          
          scaleOne = 1;
          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
                c3dTime, fpeData.projectionError, scaleOne, numberOfPhases,...
                numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).fpeasmhkz = phaseStats(p);    
          end
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Time relative to seat-off at which the distance (COM-FPE)u is
          % maximized
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
          flag_ModeBalancePointsVsCom0VsCop1   = 0;
          flag_ModeAnalyzeBalanceAlong0Across1 = 0; 
          [fpelength,dataName] = calcBalancePointDistance(fpeData.r0F0, ...
                        fpeData.u,fpeData.n,...
                        wholeBodyData(:,colComPos),...
                        [],...
                        flag_ModeBalancePointsVsCom0VsCop1,...
                        flag_ModeAnalyzeBalanceAlong0Across1,...
                        'Fpe');          
          
          timemaxfpe=calcTimeToMaximum(movementSequence,c3dTime,fpelength);

          [phaseStats, timingStats] = ...
            calcMovementSequenceDataSummary(movementSequence, ...
              c3dTime, fpelength, scaleOne,numberOfPhases,...
              numberOfInterpolationPoints); 

          for p=1:1:numberOfPhases
            subjectData(indexSubject,indexTrialsToProcess,p).timemaxfpe ...
              = phaseStats(p);    
          end
          
          
        end



      end
    end

  end

  save([frontiersDataDir,'subjectTrialPhaseMetricData',nameToeTag,'.mat'],'subjectData');
else
  load([frontiersDataDir,'subjectTrialPhaseMetricData',nameToeTag,'.mat']);
end

if(flag_loadPrecomputedGroupStatistics == 0)

  for indexTrialsToProcess = 1:1:length(trialsToProcess)

    indexTrial = 0;
    for k=1:1:length(trialTypeNames)
      if(contains(trialsToProcess{indexTrialsToProcess},trialTypeNames{k}))
        indexTrial = k;
      end
    end

    if( contains(trialTypeNames{indexTrial},'static') == 0 && indexTrial <= 6)

    disp(['Processing: ',trialsToProcess{indexTrialsToProcess}]);    
      for indexGroup = 1:1:length(groups)

        disp(['  Group: ',groups(indexGroup).name]);


        for indexPhase=1:1:numberOfPhases

          groupData(indexGroup,indexTrialsToProcess,indexPhase).name ...
            = groups(indexGroup).name;
          groupData(indexGroup,indexTrialsToProcess,indexPhase).label ...
            = groups(indexGroup).label;
          groupData(indexGroup,indexTrialsToProcess,indexPhase).subjectIndex ...
            = groups(indexGroup).index;
          groupData(indexGroup,indexTrialsToProcess,indexPhase).trialId ...
           = trialTypeNames{indexTrial};
          groupData(indexGroup,indexTrialsToProcess,indexPhase).trialTypeIndex ...
           = indexTrial;

          for indexMetric = 1:1:length(metricNameList)

            
            if(indexTrialsToProcess == 4 ...
                && indexGroup == indexGroupElderly ...
                && indexPhase == indexPhaseSeatOff2Stand ...
                && indexMetric == length(metricNameList))
              here=1;
            end

            flag_verbose=1;
            groupData(indexGroup,indexTrialsToProcess,indexPhase) ...
              = updateGroupStats(metricNameList{indexMetric},...
                  groups(indexGroup).index,...
                  subjectData(:,indexTrialsToProcess,indexPhase),...
                  groupData(indexGroup,indexTrialsToProcess,indexPhase),...
                  flag_verbose);       
            here=1;
          end
        end

      end

    end
  end
  
  save([frontiersDataDir,'groupTrialPhaseMetricData',nameToeTag,'.mat'],'groupData');
else
  load([frontiersDataDir,'groupTrialPhaseMetricData',nameToeTag,'.mat']);  
end


if(flag_writeCSVGroupTables==1)
  tableFolderName = ['../../outputData/Frontiers2020Pub/'];
  success = writeGroupMetricTable(tableFolderName, groupData, groups, ...
                                  metricNameList, metricSubFields,...
                                  trialsToProcess, trialTypeNames,...
                                  phaseNames,nameToeTag);
end                              

