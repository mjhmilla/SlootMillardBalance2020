
if(exist('flag_outerLoopMode','var') ==0)
    clc;
    close all;
    clear all;

    flag_ConvexHullWithToes0Without1          = 0;
    flag_ModePoint                            = 3;
    %0 fpe
    %1 cap
    %2 cop
    %3 com
      flag_pointToNearestFootEdge             = 0;
      flag_pointToNearestFootEdgeSubject      = 0;


    flag_ModeBalancePortrait                  = 0;
    % 0: fpe
    % 1: cap
    %      
      flag_balancePortraitSubject             = 0;        
      flag_balancePortraitSeatOff             = 0;

    flag_ModeBalancePointsVsCom0VsCop1        = 1;
    flag_ModeAnalyzeBalanceAlong0Across1      = 0;

      flag_fpeTimeSeriesPlotsSubject          = 0;
      flag_fpeTimeSeriesPlotsEnsemble         = 0;
      flag_fpeSeatOffEnsemble                 = 0;

      flag_capTimeSeriesPlotsSubject          = 0;
      flag_capTimeSeriesPlotsEnsemble         = 0;
      flag_capSeatOffEnsemble                 = 0;

    flag_ModeErrorDistance0Angle1             = 0;
      flag_fpeCapErrorTimeSeriesPlotsSubject  = 0;
      flag_fpeCapErrorTimeSeriesPlotsEnsemble = 0;
    %flag_fpeCapErrorSeatOffEnsemble      
      
    flag_ModeComVelX0VelY1VelZ2Speed3ComGpVsCop4 = 3;
            
      flag_comKinematicsTimeSeriesSubject   = 1;
      flag_comKinematicsTimeSeriesEnsemble  = 0;
      flag_comKinematicsSeatOffEnsemble     = 0;

    flag_GrfzTimeSeriesPlotsSubject         = 0;
    flag_GrfzTimeSeriesPlotsEnsemble        = 0;
    flag_GrfzSeatOffEnsemble                = 0;

end



gravityVec = [0;0;-9.81];
flag_EventStart0Reference1End2 = 1;

plotConfig;





subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE05','configE06',...
   'configE07','configE08'};

 
groups(2) = struct('index',[],'name',[],'color',[]);
groups(1).index = [1,2,3,4,5,6,7,8,9,10];
groups(1).name = 'Y';
groups(2).index = [11,12,13,14,15, 16,17];
groups(2).name = 'E';

%The rows store: min, 25%, mean, 75%, max
subjectData(length(subjectsToProcess)) =...
             struct('com2edge',zeros(7,6),...
                    'com2cop' ,zeros(7,6),...
                    'comvel'  ,zeros(7,6),...
                    'fpe2edge',zeros(7,6),...
                    'fpelen'  ,zeros(7,6),...
                    'fpewidth',zeros(7,6)); 
%'configE06',...

%%
% Setup the input/output directory structure
%%
inputDirRelative = '../../inputData';
outputDirRelative = '../../outputData';

codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);

figSubjectPointToFootEdge(length(subjectsToProcess)) = struct('h',[]);
figTrialPointToFootEdge = [];
if(flag_pointToNearestFootEdge == 1)
  figTrialPointToFootEdge = figure;
  flag_firstCall = 1;
end

figSubjectBalancePortrait(length(subjectsToProcess)) = struct('h',[]);

figBalancePortraitSeatOff = [];
if(flag_balancePortraitSeatOff == 1)
  figBalancePortraitSeatOff = figure;
end

figSubjectCom(length(subjectsToProcess)) = struct('h',[]);
figAllCom = [];
if(flag_comKinematicsTimeSeriesEnsemble==1)
  figAllCom = figure;
end

figTrailSeatOffCom = [];
if(flag_comKinematicsSeatOffEnsemble==1)
  figTrialSeatOffCom = figure;
end

figSubjectFpe(length(subjectsToProcess)) = struct('h',[]);
figAllFpe = [];
if(flag_fpeTimeSeriesPlotsEnsemble == 1)
  figAllFpe = figure;
end

 
%figTrialSeatOffFpe(numberOfDataTrials) = struct('h',[]);
figTrialSeatOffFpe =[];
if( flag_fpeSeatOffEnsemble == 1)
  figTrialSeatOffFpe = figure;
end


figSubjectCap(length(subjectsToProcess)) = struct('h',[]);
figAllCap = [];
if(flag_capTimeSeriesPlotsEnsemble == 1)
  figAllCap = figure;
end
figTrialSeatOffCap =[];
if(flag_capSeatOffEnsemble==1)
  figTrialSeatOffCap = figure;
end


figSubjectFpeCapErr(length(subjectsToProcess)) = struct('h',[]);
figAllFpeCapErr = [];
if(flag_fpeCapErrorTimeSeriesPlotsEnsemble == 1)
  figAllFpeCapErr = figure;
end

figSubjectGrfZ(length(subjectsToProcess)) = struct('h',[]);
figAllGrfZ = [];
if(flag_GrfzTimeSeriesPlotsEnsemble == 1)
  figAllGrfZ = figure;
end

figTrialGrfZSeatOff = [];
if(flag_GrfzSeatOffEnsemble==1)
  figTrialGrfZSeatOff = figure;
end

green  = [102 204 0]./255; 
blue   = [51 153 255]./255; 
orange = [255 128 0]./255;

colorElderly = ones(4,3).*green;
colorYoung   = ones(4,3).*[0,0,0];

groups(1).color = colorYoung(1,:);
groups(2).color = colorElderly(1,:);

for indexSubject = 1:1:length(subjectsToProcess)

  if(indexSubject > 1)
    flag_firstCall=0;
  end
    
  if(flag_comKinematicsTimeSeriesSubject == 1)
    figSubjectCom(      indexSubject).h = figure;
  end
  if(flag_fpeTimeSeriesPlotsSubject==1)
    figSubjectFpe(      indexSubject).h = figure;
  end
  if(flag_capTimeSeriesPlotsSubject==1)
    figSubjectCap(      indexSubject).h = figure;
  end  
  if(flag_GrfzTimeSeriesPlotsSubject==1)
    figSubjectGrfZ(     indexSubject).h = figure;
  end
  if(flag_fpeCapErrorTimeSeriesPlotsSubject==1)
    figSubjectFpeCapErr(      indexSubject).h = figure;
  end
  if(flag_balancePortraitSubject == 1)
    figSubjectBalancePortrait(indexSubject).h = figure;
  end
  if(flag_pointToNearestFootEdgeSubject==1)
    figSubjectPointToFootEdge(indexSubject).h = figure;
  end

  colorFpe    = [0,0,0];
  colorFpeCap = [0,0,0];
  colorCap    = [0,0,0];
  colorCom    = [0,0,0];
  colorCop    = [0,0,0];  
  colorGrfz   = [0,0,0];
  flag_young = 0;
  if( isempty(strfind(subjectsToProcess{indexSubject},'E'))==0)
    colorFpe    = colorElderly(1,:);
    colorFpeCap = colorElderly(1,:);  
    colorCap    = colorElderly(1,:);    
    colorCom    = colorElderly(1,:);
    colorCop    = colorElderly(1,:);    
    colorGrfz   = colorElderly(1,:);
  else
    colorFpe    = colorYoung(1,:);    
    colorFpeCap = colorYoung(1,:);  
    colorCap    = colorYoung(1,:);    
    colorCom    = colorYoung(1,:);
    colorCop    = colorYoung(1,:);    
    colorGrfz   = colorYoung(1,:);
    flag_young = 1;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1. Configure the list of input/output files for this subject
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
     
  numberOfTrials = length(inputC3DFiles); 
  
  subjectId = subjectId(2:3);
  if(contains(subjectId,'0'))
    idx = strfind(subjectId,'0');
    if(idx == 1)
      subjectId = subjectId(2);
    end
  end
  
  disp(['Processing: ', subjectId]);  
  
  for indexTrial = 1:1:numberOfTrials
    

    %Get the trial folder
    trialFolder = outputTrialFolders{indexTrial};    
        
    %Load the processed trial data
    
    c3DFileName       = updateFileExtension(inputC3DFiles{indexTrial},'mat');
    
    if( contains(c3DFileName,'static') == 0)
      
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

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Plot the FPE data at the sit-to-stand transition
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if(indexTrial <= 6)


          timeLabel = 'Time (s)';
          forceLabel= 'Force (N)';
          angleLabel= 'Angle (deg)';
          distanceLabel = 'Distance (cm)';
          velocityLabel = 'Velocity (cm/s)';
          scaleTime = 1;
          scaleForce = 1;
          scaleDistance = 100;
          scaleAngle    = 180/pi;
          scaleVelocity = 100;


          [row,col] = find(subPlotPanelIndex == (indexTrial-offsetTrialIndex));          
          subplotPosition = reshape(subPlotPanel(row,col,:),1,4);
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Balance point to foot edge
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

          if(   flag_pointToNearestFootEdge         == 1 ...
             || flag_pointToNearestFootEdgeSubject  == 1)

            data      = [];
            scaleData = [];
            colorData = [];
            yLabelData =[];
            figTitle = '';
            axisYLim = [];
            %0 fpe
            %1 cap
            %2 cop
            %3 com
            switch flag_ModePoint
              case 0
                data  = -fpe2FootConvexHullDist.distance;
                scaleData = scaleDistance;
                colorData = colorFpe;
                yLabelData = distanceLabel;
                figTitle = 'Fpe--Foot-Edge';
                axisYLim = [-7, 16];
              case 1
                data  = -cap2FootConvexHullDist.distance;
                scaleData = scaleDistance;                
                colorData = colorCap;
                yLabelData = distanceLabel;
                figTitle = 'Cap--Foot-Edge';    
                axisYLim = [-7,16];
              case 2
                data  = -cop2FootConvexHullDist.distance;
                scaleData = scaleDistance;     
                colorData = colorCop;
                yLabelData = distanceLabel; 
                figTitle = 'Cop--Foot-Edge';      
                axisYLim = [0,16];                
              case 3
                data  = -comgp2FootConvexHullDist.distance;
                scaleData = scaleDistance;                
                colorData = colorCom;                
                yLabelData = distanceLabel;
                figTitle = 'Com--Foot-Edge'; 
                axisYLim = [-10,17];
              otherwise assert(0);
            end

            if(flag_pointToNearestFootEdge == 1)
              yLabelPos = axisYLim(1,1)+1.5;
              %if(flag_young == 0)
              %  yLabelPos = 14;
              %end
              
              flag_IntervalMode = 1;
              [figTrialPointToFootEdge, dataSummary] = ...
                plotBoxWhiskerEventData(...
                   figTrialPointToFootEdge, subplotPosition, ...
                   movementSequence, ...
                   indexSubject, 1,...              
                   data, scaleData, ...
                   colorData,  ...
                   0.33,...
                   subjectId, yLabelPos,...
                   '', yLabelData, ...
                   [figTitle,': ',trialTypeNames{indexTrial}], ...
                   flag_IntervalMode, flag_firstCall);
              axisLim = axis;
              axisLim(1)=-0.5;
              axisLim(2)=length(subjectsToProcess)+2.5;
              axisLim(3)=axisYLim(1);
              axisLim(4)=axisYLim(2);
              axis(axisLim);
              set(gca,'XTickLabel',[]); 
              
              dataSummaryVec = [ dataSummary.min;...
                                 dataSummary.p25;...
                                 dataSummary.mean;...
                                 dataSummary.p75;...
                                 dataSummary.max;...
                                 dataSummary.events'];
              switch flag_ModePoint
                case 0
                  subjectData(indexSubject).fpe2edge(:,indexTrial) = ...
                    dataSummaryVec;
                case 1
                  subjectData(indexSubject).cap2edge(:,indexTrial) = ...
                    dataSummaryVec;
                case 2
                  subjectData(indexSubject).cop2edge(:,indexTrial) = ...
                    dataSummaryVec;         
                case 3
                  subjectData(indexSubject).com2edge(:,indexTrial) = ...
                    dataSummaryVec;
                otherwise assert(0);              
              end
            end

            %figSubjectPointToFootEdge
            if(flag_pointToNearestFootEdgeSubject ==  1)
              figSubjectPointToFootEdge(indexSubject).h = ...
                plotTimeSeriesData(...
                  figSubjectPointToFootEdge(indexSubject).h, subplotPosition, ...
                  movementSequence, ...
                  c3dTime, scaleTime, ...
                  data, scaleData, ...
                  colorData,  ...
                  timeLabel, ...
                  yLabelData, ...
                  [figTitle,'(',subjectId,'): ',trialTypeNames{indexTrial}]);
            end
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Balance Portrait
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(   flag_balancePortraitSubject == 1 ...
             || flag_balancePortraitSeatOff == 1)

            balancePoint = [];
            balanceTangentDir = [];
            balanceNormalDir = [];
            lineColor = [];
            figTitle = '';
            switch flag_ModeBalancePortrait
              case 0
                balancePoint = fpeData.r0F0;
                balanceTangentDir = fpeData.u;
                balanceNormalDir  = fpeData.n;              
                lineColor = colorFpe;
                figTitle = 'Fpe BP ';
              case 1
                balancePoint = capData.r0F0;
                balanceTangentDir = capData.u;
                balanceNormalDir  = capData.n;              
                lineColor = colorCap;
                figTitle = 'Cap BP ';
            end
            c3dFootMarkersLeftNames = {'L_FAL','L_FM1','L_FM2','L_FM5','L_TAM','L_FCC'};
            c3dFootMarkersRightNames = {'R_FAL','R_FM1','R_FM2','R_FM5','R_TAM','R_FCC'};
            
            lineColorBalancePoint = [1,1,1/0.25].*0.25;
            lineColorBalanceDir   = [1,1,1/0.75].*0.75;
            lineColorCop          = [1,0,0];
            lineColorCom          = [0.5,0.5,0.5];
            footPatchColor        = [1,1,1].*0.9;
            
            if(flag_balancePortraitSubject == 1)
              figSubjectBalancePortrait(indexSubject).h ...
              = plotBalancePortrait(...
                  figSubjectBalancePortrait(indexSubject).h,...
                  subplotPosition,...
                  indexSubject, subjectId,...
                  c3dFootMarkersRightNames, c3dFootMarkersLeftNames, ...                  
                  movementSequence, ...
                  c3dTime, c3dMarkers, ...
                  wholeBodyData(:,colComPos), ...
                  balancePoint, balanceTangentDir, ...
                  c3dGrf(index_FeetForcePlate),...
                  0,...
                  lineColorCom,lineColorBalancePoint,lineColorBalanceDir,...
                  lineColorCop,...
                  footPatchColor,...
                  [figTitle, '(',subjectId,'): ', trialTypeNames{indexTrial}]);
            end
            
            if(flag_balancePortraitSeatOff == 1)
              figBalancePortraitSeatOff ...
                = plotBalancePortraitAtPointInTime(...
                    figBalancePortraitSeatOff,...
                    subplotPosition,...
                    indexSubject, subjectId,...
                    c3dFootMarkersRightNames, c3dFootMarkersLeftNames, ...                  
                    movementSequence, ...
                    scaleTime,scaleDistance,scaleAngle,scaleForce,...
                    c3dTime, c3dMarkers, ...
                    wholeBodyData(:,colComPos), ...
                    balancePoint, balanceTangentDir, balanceNormalDir,...
                    c3dGrf(index_FeetForcePlate),...
                    lineColorCom,lineColorBalancePoint,lineColorBalanceDir,...
                    lineColorCop,...
                    footPatchColor,...
                    '',distanceLabel,...
                    [figTitle,': ', trialTypeNames{indexTrial}]);
              axisLim = axis;
              axisLim(1) = 0;
              axisLim(2) = (length(subjectsToProcess)*0.05 +0.025).*scaleDistance;
              axisLim(3) = (-0.1)*scaleDistance;
              axisLim(4) = (0.2)*scaleDistance;
              axis(axisLim);
            end
           
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Whole Body Kinematic Portrait
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          if(    flag_comKinematicsTimeSeriesSubject   == 1 ...
              || flag_comKinematicsTimeSeriesEnsemble  == 1 ...
              || flag_comKinematicsSeatOffEnsemble     == 1)

            errorVec = [];
            errorLabel = [];
            errorScale = [];
            errorName = [];
            figTitle = '';
            switch flag_ModeComVelX0VelY1VelZ2Speed3ComGpVsCop4
              case 0
                errorVec = wholeBodyData(:,colComVel(1,1));
                errorLabel= velocityLabel;
                errorScale= scaleVelocity;
                figTitle = 'Com Vel-X ';
              case 1
                errorVec = wholeBodyData(:,colComVel(1,2));
                errorLabel= velocityLabel;
                errorScale= scaleVelocity;
                figTitle = 'Com Vel-Y ';                
              case 2
                errorVec = wholeBodyData(:,colComVel(1,3));
                errorLabel= velocityLabel;
                errorScale= scaleVelocity;
                figTitle = 'Com Vel-Z ';
              case 3
                errorVec = sum( wholeBodyData(:,colComVel).^2, 2).^0.5;
                errorLabel= velocityLabel;
                errorScale= scaleVelocity;
                figTitle = 'Com Vel ';                
              case 4
                errorVec = sum( (c3dGrf(index_FeetForcePlate).cop(:,1:2) ...
                                     - wholeBodyData(:,colComPos(1,1:2))).^2 ,2).^0.5;
                errorLabel= distanceLabel;
                errorScale= scaleDistance;
                figTitle = 'Cop-Com ';                
              otherwise assert(0);
            end
            
            if(flag_comKinematicsTimeSeriesSubject==1)
              figSubjectCom(indexSubject).h =...
                plotTimeSeriesData(...
                  figSubjectCom(indexSubject).h, ...
                  subplotPosition,...
                  movementSequence,...
                  c3dTime, scaleTime,...
                  errorVec, errorScale,...
                  colorCom,...
                  timeLabel,errorLabel,...
                  [figTitle,'(',subjectId,'): ',trialTypeNames{indexTrial}]);
              axis tight;              
            end
            if(flag_comKinematicsTimeSeriesEnsemble==1)
              figAllCom = plotTimeSeriesData(...
                  figAllCom, ...
                  subplotPosition,...
                  movementSequence,...
                  c3dTime, scaleTime,...
                  errorVec, errorScale,...
                  colorCom,...
                  timeLabel,errorLabel,...
                  [figTitle,trialTypeNames{indexTrial}]);
              axis tight              
            end
            if(flag_comKinematicsSeatOffEnsemble==1)
              figTrialSeatOffCom = plotEventData( ...
                figTrialSeatOffCom, subplotPosition,...
                movementSequence,...
                indexSubject, 1,...
                errorVec, errorScale,...
                colorCom,...
                subjectId, 1,...
                '',errorLabel,...
                [figTitle,' At Seat Off',trialTypeNames{indexTrial}],...
                flag_EventStart0Reference1End2);

              axis tight;
              axisLim = axis;            
              axisLim(1) = 0;
              axisLim(2) = length(subjectsToProcess)+1;
              switch flag_ModeComVelX0VelY1VelZ2Speed3ComGpVsCop4
                case 0
                  axisLim(3)=0;
                  axisLim(4)=60;
                case 1
                  axisLim(3)=-10;
                  axisLim(4)=10;                
                case 2
                  axisLim(3)=-5;
                  axisLim(4)=40;                
                case 3
                  axisLim(3)=0;
                  axisLim(4)=70;                
                case 4
                  axisLim(3)=0;
                  axisLim(4)=15;                
                otherwise assert(0);                
              end
              axis(axisLim);
              set(gca,'XTickLabel',[]);    
            end
    
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Whole Body Kinematic Portrait : Seat off
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    

          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe Time Series
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          
          if(    flag_fpeTimeSeriesPlotsSubject        == 1 ...
              || flag_fpeTimeSeriesPlotsEnsemble       == 1 ...
              || flag_fpeSeatOffEnsemble               == 1 )     

            [errorVec,errorName] = calcBalancePointDistance(fpeData.r0F0, ...
                          fpeData.u,fpeData.n,...
                          wholeBodyData(:,colComPos),...
                          c3dGrf(index_FeetForcePlate).cop,...
                          flag_ModeBalancePointsVsCom0VsCop1,...
                          flag_ModeAnalyzeBalanceAlong0Across1,...
                          'Fpe');
                   
            %flag_ModeBalancePointsVsCom0VsCop1: scale/label is the same
            errorScale= scaleDistance;
            errorLabel=distanceLabel;

            if(flag_fpeTimeSeriesPlotsSubject==1)
              figSubjectFpe(indexSubject).h = ...
                plotTimeSeriesData(...
                          figSubjectFpe(indexSubject).h,...
                          subplotPosition,...
                          movementSequence,...
                          c3dTime, scaleTime,...
                          errorVec, errorScale,...
                          colorFpe, ...
                          timeLabel, errorLabel,...
                          [errorName,': (',subjectId,'): ',trialTypeNames{indexTrial}]);              
            end 

            if(flag_fpeTimeSeriesPlotsEnsemble==1)
              figAllFpe = plotTimeSeriesData(...
                          figAllFpe,subplotPosition,...
                          movementSequence,...
                          c3dTime, scaleTime,...
                          errorVec, errorScale,...
                          colorFpe, ...
                          timeLabel, errorLabel,...
                          [errorName,': ',trialTypeNames{indexTrial}]);   
            end 

            if(flag_fpeSeatOffEnsemble ==1 )
              figTrialSeatOffFpe = plotEventData(...
                          figTrialSeatOffFpe,subplotPosition,...
                          movementSequence,...
                          indexSubject, 1,...
                          errorVec, errorScale,...
                          colorFpe, subjectId, 1,...
                          '',errorLabel,...
                          [errorName,' Seat-Off: ',trialTypeNames{indexTrial}],...
                          flag_EventStart0Reference1End2);             
              set(gca,'XTickLabel',[]);                
              axis tight;
              axisLim = axis;
              axisLim(1) = 0;
              axisLim(2) = length(subjectsToProcess)+1;
              switch(flag_ModeBalancePointsVsCom0VsCop1)
                case 0
                  switch(flag_ModeAnalyzeBalanceAlong0Across1)
                    case 0
                      axisLim(3) =-0.5;
                      axisLim(4) = 20.;
                    case 1
                      axisLim(3) =-5;
                      axisLim(4) = 5;                      
                    otherwise assert(0)                      
                  end
                case 1
                  switch(flag_ModeAnalyzeBalanceAlong0Across1)
                    case 0
                      axisLim(3) =-12;
                      axisLim(4) = 12;
                    case 1
                      axisLim(3) =-5;
                      axisLim(4) = 5;                      
                    otherwise assert(0)                      
                  end                  
                otherwise assert(0)                  
              end

              axis(axisLim);
              set(gca,'XTickLabel',[]);
            end

          end
          
          
 
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Cap - Seat-Off
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          
          if(    flag_capTimeSeriesPlotsSubject   == 1 ...
              || flag_capTimeSeriesPlotsEnsemble  == 1 ...
              || flag_capSeatOffEnsemble          == 1)
             
            [errorVec,errorName] = calcBalancePointDistance(capData.r0F0, ...
                          capData.u,capData.n,...
                          wholeBodyData(:,colComPos),...
                          c3dGrf(index_FeetForcePlate).cop,...
                          flag_ModeBalancePointsVsCom0VsCop1,...
                          flag_ModeAnalyzeBalanceAlong0Across1,...
                          'Cap');
                   
            %flag_ModeBalancePointsVsCom0VsCop1: scale/label is the same
            errorScale= scaleDistance;
            errorLabel=distanceLabel;

            if(flag_capTimeSeriesPlotsSubject==1)
              figSubjectCap(indexSubject).h = ...
                plotTimeSeriesData(...
                          figSubjectCap(indexSubject).h,...
                          subplotPosition,...
                          movementSequence,...
                          c3dTime, scaleTime,...
                          errorVec, errorScale,...
                          colorCap, ...
                          timeLabel, errorLabel,...
                          [errorName,': (',subjectId,'): ',trialTypeNames{indexTrial}]);              
            end 

            if(flag_capTimeSeriesPlotsEnsemble==1)
              figAllCap = plotTimeSeriesData(...
                          figAllCap,subplotPosition,...
                          movementSequence,...
                          c3dTime, scaleTime,...
                          errorVec, errorScale,...
                          colorCap, ...
                          timeLabel, errorLabel,...
                          [errorName,': ',trialTypeNames{indexTrial}]);   
            end 

            if(flag_capSeatOffEnsemble ==1 )
              figTrialSeatOffCap = plotEventData(...
                          figTrialSeatOffCap,subplotPosition,...
                          movementSequence,...
                          indexSubject, 1,...
                          errorVec, errorScale,...
                          colorCap, subjectId, 1,...
                          '',errorLabel,...
                          [errorName,' Seat-Off: ',trialTypeNames{indexTrial}],...
                          flag_EventStart0Reference1End2);             
              set(gca,'XTickLabel',[]);                
              axis tight;
              axisLim = axis;
              axisLim(1) = 0;
              axisLim(2) = length(subjectsToProcess)+1;
              switch(flag_ModeBalancePointsVsCom0VsCop1)
                case 0
                  switch(flag_ModeAnalyzeBalanceAlong0Across1)
                    case 0
                      axisLim(3) =-0.5;
                      axisLim(4) = 16.5;
                    case 1
                      axisLim(3) =-11;
                      axisLim(4) = 11;                      
                    otherwise assert(0)                      
                  end
                case 1
                  switch(flag_ModeAnalyzeBalanceAlong0Across1)
                    case 0
                      axisLim(3) =-11;
                      axisLim(4) = 11;
                    case 1
                      axisLim(3) =-5;
                      axisLim(4) = 5;                      
                    otherwise assert(0)                      
                  end                  
                otherwise assert(0)                  
              end

              axis(axisLim);
              set(gca,'XTickLabel',[]);
            end
          end          
          

          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Fpe-Cap Time Series
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          
          if(     flag_fpeCapErrorTimeSeriesPlotsSubject == 1 ...
              ||  flag_fpeCapErrorTimeSeriesPlotsEnsemble == 1)
            errorVec = [];
            errorScale = [];
            errorLabel = '';
            switch flag_ModeErrorDistance0Angle1
              case 0
                errorVec =  sum( (fpeData.r0F0-capData.r0F0).*fpeData.u, 2);
                errorScale = scaleDistance;
                errorLabel = distanceLabel;
              case 1
                errorVec = acos( sum( fpeData.n .* capData.n, 2) );
                errorScale=scaleAngle;
                errorLabel=angleLabel;
              otherwise assert(0); 
            end
            if(flag_fpeCapErrorTimeSeriesPlotsSubject == 1)
              figSubjectFpeCapErr(indexSubject).h = ...
                plotTimeSeriesData(...
                  figSubjectFpeCapErr(indexSubject).h, ...
                  subplotPosition,...
                  movementSequence,...
                  c3dTime, scaleTime,...
                  errorVec,errorScale,...
                  colorFpeCap,...
                  timeLabel, errorLabel, ...
                  ['Fpe-Cap Err. ','(',subjectId,'): ',trialTypeNames{indexTrial}]);
            end
            if(flag_fpeCapErrorTimeSeriesPlotsEnsemble == 1)
               figAllFpeCapErr = ...
                plotTimeSeriesData( ...
                  figAllFpeCapErr, ...
                  subplotPosition,...
                  movementSequence,...
                  c3dTime, scaleTime,...
                  errorVec,errorScale,...
                  colorFpeCap,...
                  timeLabel, errorLabel, ...
                  ['Fpe-Cap Err.: ', trialTypeNames{indexTrial}]);  
            end
          end 
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Ground Forces Time Series
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

          if(flag_GrfzTimeSeriesPlotsSubject == 1)
            figSubjectGrfZ(indexSubject).h = plotTimeSeriesData(...
                        figSubjectGrfZ(indexSubject).h, ...
                        subplotPosition,... 
                        movementSequence, ...
                        c3dTime, scaleTime, ...
                        c3dGrf(index_ChairForcePlate).force(:,3), scaleForce,...
                        colorGrfz, ...
                        timeLabel, forceLabel,figGrfZTitle);  
          end

          if(flag_GrfzTimeSeriesPlotsEnsemble == 1)
            figAllGrfZ = plotTimeSeriesData(figAllGrfZ,...
                        subplotPosition,... 
                        movementSequence, ...
                        c3dTime, scaleTime, ...
                        c3dGrf(index_ChairForcePlate).force(:,3), scaleForce,...
                        colorGrfz, ...
                        timeLabel, forceLabel,...
                        ['GRFz: ', trialTypeNames{indexTrial}]);  
          end

                    
          if(flag_GrfzSeatOffEnsemble==1)
            figTrialGrfZSeatOff = plotEventData(...
                    figTrialGrfZSeatOff, subplotPosition,... 
                    movementSequence,...
                    indexSubject, 1,...
                    c3dGrf(index_ChairForcePlate).force(:,3),scaleForce,...
                    colorGrfz, subjectId, 2.5,...
                    '',forceLabel, ...
                    ['GRFz Seat-Off:', trialTypeNames{indexTrial}],...
                    flag_EventStart0Reference1End2); 
            set(gca,'XTickLabel',[]);
          end


      end
      
      
    end
    
  end
  
  
  if(flag_pointToNearestFootEdgeSubject == 1)
    analysisType = '';
    switch flag_ModePoint
      case 0
        analysisType = 'FpeFootEdge';
      case 1
        analysisType = 'CapFootEdge';
      case 2
        analysisType = 'CopFootEdge';      
      case 3
        analysisType = 'ComFootEdge';      
      otherwise assert(0);
    end  
    if flag_ConvexHullWithToes0Without1 == 1
      analysisType = [analysisType,'NoToes'];
    end    
    figure(figSubjectPointToFootEdge(indexSubject).h);
    configPlotExporter;
    print('-dpdf',['../../plots/footEdge/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  
  if(flag_balancePortraitSubject==1 )

    analysisType='';
    switch flag_ModeBalancePortrait      
      case 0
        analysisType = ['FpeBalancePortrait'];
      case 1
        analysisType = ['CapBalancePortrait'];        
    end
    
    
    figure(figSubjectBalancePortrait(indexSubject).h);
    configPlotExporter;
    print('-dpdf',['../../plots/balancePortrait/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);  
    
  end



  if(flag_comKinematicsTimeSeriesSubject==1)
    figure(figSubjectCom(indexSubject).h);
    configPlotExporter;    
    analysisType = '';
    switch flag_ModeComVelX0VelY1VelZ2Speed3ComGpVsCop4
      case 0
        analysisType = 'ComVelX';
      case 1
        analysisType = 'ComVelY';
      case 2
        analysisType = 'ComVelZ';
      case 3
        analysisType = 'ComVel';
      case 4
        analysisType = 'ComGPCop';
      otherwise assert(0);
    end                
    print('-dpdf',['../../plots/com/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);
  end  
  
  if(flag_fpeTimeSeriesPlotsSubject==1)
    figure(figSubjectFpe(indexSubject).h);
    configPlotExporter;    
    analysisType = 'FpeVs';
    switch flag_ModeBalancePointsVsCom0VsCop1
      case 0
        analysisType = [analysisType,'Com'];
      case 1
        analysisType = [analysisType,'Cop'];
    end    
    switch flag_ModeAnalyzeBalanceAlong0Across1
      case 0
        analysisType = [analysisType,'Length'];
      case 1
        analysisType = [analysisType,'Width'];
    end            
    print('-dpdf',['../../plots/fpe/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  if(flag_capTimeSeriesPlotsSubject==1)
    figure(figSubjectCap(indexSubject).h);
    configPlotExporter;  
    analysisType = 'CapVs';
    switch flag_ModeBalancePointsVsCom0VsCop1
      case 0
        analysisType = [analysisType,'Com'];
      case 1
        analysisType = [analysisType,'Cop'];
    end    
    switch flag_ModeAnalyzeBalanceAlong0Across1
      case 0
        analysisType = [analysisType,'Length'];
      case 1
        analysisType = [analysisType,'Width'];
    end            
    print('-dpdf',['../../plots/cap/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  
  if(flag_fpeCapErrorTimeSeriesPlotsSubject==1)
    figure(figSubjectFpeCapErr(indexSubject).h);
    configPlotExporter;    
    errorType = '';
    switch flag_ModeErrorDistance0Angle1
      case 0
        errorType = 'DistanceError';
      case 1
        errorType = 'AngleError';
      otherwise assert(0); 
    end
    
    print('-dpdf',['../../plots/fpeCap/fig_FpeCap',errorType,'S2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  
  if(flag_GrfzTimeSeriesPlotsSubject==1)
    figure(figSubjectGrfZ(indexSubject).h);
    configPlotExporter;
    print('-dpdf',['../../plots/grfz/fig_GrfZS2SQuietTimeSeries_',subjectId,'.pdf']);
  end
end

%%
%
% Add in the group statistics to each of the plots
%
%%
% subjectData(length(subjectsToProcess)) =...
%              struct('com2edge',zeros(5,6),...
%                     'com2cop' ,zeros(5,6),...
%                     'comvel'  ,zeros(5,6),...
%                     'fpe2edge',zeros(5,6),...
%                     'fpelen'  ,zeros(5,6),...
%                     'fpewidth',zeros(5,6)); 


if(flag_pointToNearestFootEdge == 1)
  for indexTrial = 1:1:numberOfTrials

    if(indexTrial > offsetTrialIndex && indexTrial <= 6) 


    %0 fpe
    %1 cap
    %2 cop
    %3 com
    
    
      groupDataSummary(length(groups)) = struct('data',zeros(7,1),'n',0);
      for indexGroup = 1:1:length(groups)
        groupDataSummary(indexGroup).n = length(groups(indexGroup).index);
        groupDataSummary(indexGroup).data = zeros(7,1);
        for groupMembers = 1:1:length(groups(indexGroup).index)
          
          indexSubject = groups(indexGroup).index(1,groupMembers);
          
          switch flag_ModePoint
            case 0
              groupDataSummary(indexGroup).data=groupDataSummary(indexGroup).data...
                                    + subjectData(indexSubject).fpe2edge(:,indexTrial); 
            case 1
              groupDataSummary(indexGroup).data=groupDataSummary(indexGroup).data...
                                    + subjectData(indexSubject).cap2edge(:,indexTrial); 
              
            case 2
              groupDataSummary(indexGroup).data=groupDataSummary(indexGroup).data...
                                    + subjectData(indexSubject).cop2edge(:,indexTrial); 
              
            case 3
              groupDataSummary(indexGroup).data=groupDataSummary(indexGroup).data...
                                    + subjectData(indexSubject).com2edge(:,indexTrial); 
              
            otherwise assert(0);
          end
        end
        groupDataSummary(indexGroup).data = ...
          groupDataSummary(indexGroup).data./groupDataSummary(indexGroup).n;        

        [row,col] = find(subPlotPanelIndex == (indexTrial-offsetTrialIndex));          
        subplotPosition = reshape(subPlotPanel(row,col,:),1,4);

        axisLim = axis;
        figTrialPointToFootEdge = plotAntsOnALog(...
          figTrialPointToFootEdge, subplotPosition,...
            length(subjectsToProcess)+indexGroup, ...
            groupDataSummary(indexGroup).data(1,1),...
            groupDataSummary(indexGroup).data(3,1),...
            groupDataSummary(indexGroup).data(5,1),...
            groupDataSummary(indexGroup).data(2,1),...
            groupDataSummary(indexGroup).data(4,1),...
            [groupDataSummary(indexGroup).data(6,1),...
             groupDataSummary(indexGroup).data(7,1)],...
            {'o','o'},[1,1,1;groups(indexGroup).color],...
            0.33,groups(indexGroup).color);
        
          xDataMid = length(subjectsToProcess)+indexGroup;
          axis(axisLim);
          yLabelPos = axisLim(3)+1.5;
          
          text( xDataMid, yLabelPos, groups(indexGroup).name,...
            'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');  
          hold on;
      end   
      axisLim = axis;
      plot([1;1].*(length(subjectsToProcess)+0.5),...
           [axisLim(3);axisLim(4)],'-','Color',[1,1,1],'LineWidth',1.5);
      hold on;
      plot([1;1].*(length(subjectsToProcess)+0.5),...
           [axisLim(3);axisLim(4)],'-','Color',[1,1,1].*0,'LineWidth',0.5);
      hold on;

    end
  end
end


if(flag_balancePortraitSeatOff==1 )

  analysisType='';
  switch flag_ModeBalancePortrait      
    case 0
      analysisType = ['FpeBalancePortraitSeatOff'];
    case 1
      analysisType = ['CapBalancePortraitSeatOff'];        
  end
  figure(figBalancePortraitSeatOff);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);  

end


if(flag_pointToNearestFootEdge == 1)
  analysisType = '';
  switch flag_ModePoint
    case 0
      analysisType = 'FpeFootEdge';
    case 1
      analysisType = 'CapFootEdge';
    case 2
      analysisType = 'CopFootEdge';      
    case 3
      analysisType = 'ComFootEdge';      
    otherwise assert(0);
  end  
  if flag_ConvexHullWithToes0Without1 == 1
    analysisType = [analysisType,'NoToes'];
  end
  
  figure(figTrialPointToFootEdge);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);
  figure(figTrialPointToFootEdge);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);
end

if(flag_comKinematicsTimeSeriesEnsemble==1 ...
    || flag_comKinematicsSeatOffEnsemble==1)  
  analysisType = '';
  switch flag_ModeComVelX0VelY1VelZ2Speed3ComGpVsCop4
    case 0
      analysisType = 'ComVelX';
    case 1
      analysisType = 'ComVelY';
    case 2
      analysisType = 'ComVelZ';
    case 3
      analysisType = 'ComVel';
    case 4
      analysisType = 'ComGPCop';
    otherwise assert(0);
  end  
  if(flag_comKinematicsTimeSeriesEnsemble==1)
    figure(figAllCom);
    configPlotExporter;      
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);
  end
  if(flag_comKinematicsSeatOffEnsemble==1)
    figure(figTrialSeatOffCom)
    configPlotExporter;        
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);  
  end
end 


if(flag_fpeTimeSeriesPlotsEnsemble==1 || flag_fpeSeatOffEnsemble ==1)

  analysisType = 'FpeVs';
  switch flag_ModeBalancePointsVsCom0VsCop1
    case 0
      analysisType = [analysisType,'Com'];
    case 1
      analysisType = [analysisType,'Cop'];
  end    
  switch flag_ModeAnalyzeBalanceAlong0Across1
    case 0
      analysisType = [analysisType,'Length'];
    case 1
      analysisType = [analysisType,'Width'];
  end        
  if(flag_fpeTimeSeriesPlotsEnsemble ==1)
    figure(figAllFpe);
    configPlotExporter;  
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);  
  end
  if(flag_fpeSeatOffEnsemble ==1)
    figure(figTrialSeatOffFpe);    
    configPlotExporter;  
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);        
  end
end

if(flag_capSeatOffEnsemble == 1 || flag_capTimeSeriesPlotsEnsemble==1)
  analysisType = 'CapVs';
  switch flag_ModeBalancePointsVsCom0VsCop1
    case 0
      analysisType = [analysisType,'Com'];
    case 1
      analysisType = [analysisType,'Cop'];
  end    
  switch flag_ModeAnalyzeBalanceAlong0Across1
    case 0
      analysisType = [analysisType,'Length'];
    case 1
      analysisType = [analysisType,'Width'];
  end   
  
  if(flag_capSeatOffEnsemble == 1)
    figure(figTrialSeatOffCap);
    configPlotExporter;  
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);  
  end
  if(flag_capTimeSeriesPlotsEnsemble==1)
    figure(figAllCap);
    configPlotExporter;  
    print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);  
  end
end

if(flag_fpeCapErrorTimeSeriesPlotsEnsemble==1)
  figure(figAllFpeCapErr);
  configPlotExporter;
    errorType = '';
    switch flag_ModeErrorDistance0Angle1
      case 0
        errorType = 'DistanceError';
      case 1
        errorType = 'AngleError';
      otherwise assert(0); 
    end    
    print('-dpdf',['../../plots/fig_FpeCap',errorType,'S2SQuietTimeSeries.pdf'])    
end  
if(flag_GrfzTimeSeriesPlotsEnsemble==1)
  figure(figAllGrfZ);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_GrfZS2SQuietTimeSeries','.pdf']);
end
if(flag_GrfzSeatOffEnsemble==1)
  figure(figTrialGrfZSeatOff);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_GrfZS2SQuietSeatOff','.pdf']);
end
