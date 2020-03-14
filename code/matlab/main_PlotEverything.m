
if(exist('flag_outerLoopMode','var') ==0)
    clc;
    close all;
    clear all;

    flag_ModeBalancePortrait = 0;
    % 0: fpe
    % 1: cap
    %
    
    flag_balancePortraitSubject         = 1;    
    flag_balancePortraitEnsemble        = 0;
    flag_balancePortraitSeatOffEnsemble = 0;
    
    flag_ModeBalancePointsVsCom0VsCop1   = 1;
    flag_ModeAnalyzeBalanceAlong0Across1 = 0;

      flag_fpeTimeSeriesPlotsSubject   = 0;
      flag_fpeTimeSeriesPlotsEnsemble  = 0;
      flag_fpeSeatOffEnsemble          = 0;

      flag_capTimeSeriesPlotsSubject   = 0;
      flag_capTimeSeriesPlotsEnsemble  = 0;
      flag_capSeatOffEnsemble          = 0;

    flag_fpeCapErrorTimeSeriesPlotsSubject   = 0;
    flag_fpeCapErrorTimeSeriesPlotsEnsemble  = 0;
    %flag_fpeCapErrorSeatOffEnsemble      
      
    flag_ModeComVel0ComGpVsCop1              = 0;
            
      flag_comKinematicsTimeSeriesSubject  = 0;
      flag_comKinematicsTimeSeriesEnsemble = 0;
      flag_comKinematicsSeatOffEnsemble    = 0;

    flag_GrfzTimeSeriesPlotsSubject  = 0;
    flag_GrfzTimeSeriesPlotsEnsemble = 0;
    flag_GrfzSeatOffEnsemble         = 0;

end



gravityVec = [0;0;-9.81];

plotConfig;

subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE05', ...
   'configE06','configE07','configE08'};

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


figSubjectBalancePortrait(length(subjectsToProcess)) = struct('h',[]);
figAllBalancePortrait = [];
if(flag_balancePortraitEnsemble==1)
  figAllBalancePortrait = figure;
end

figTrialSeatOffBalancePortrait = [];
if(flag_balancePortraitSeatOffEnsemble==1)
  figTrialSeatOffBalancePortrait = figure;
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

colorElderly = [0,0,0;...
                1,0,0;...
                1,0,1;...
                0,0,1];
colorYoung = [0.5,0.5,0.5;...
              1.0,0.5,0.5;...
              1.0,0.5,1.0;...
              0.5,0.5,1.0];


for indexSubject = 1:1:length(subjectsToProcess)


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
  
  colorFpe    = [0,0,0];
  colorFpeCap = [0,0,0];
  colorCap    = [0,0,0];
  colorCom    = [0,0,0];
  colorGrfz   = [0,0,0];
  
  if( isempty(strfind(subjectsToProcess{indexSubject},'E'))==0)
    colorFpe    = colorElderly(1,:);
    colorFpeCap = colorElderly(1,:);  
    colorCap    = colorElderly(1,:);    
    colorCom    = colorElderly(1,:);
    colorGrfz   = colorElderly(1,:);
  else
    colorFpe    = colorYoung(1,:);    
    colorFpeCap = colorYoung(1,:);  
    colorCap    = colorYoung(1,:);    
    colorCom    = colorYoung(1,:);
    colorGrfz   = colorYoung(1,:);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1. Configure the list of input/output files for this subject
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
     
  numberOfTrials = length(inputC3DFiles); 
  
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
      %Plot the FPE data at the sit-to-stand transition
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if(indexTrial <= 6)
        for figType = 1:1:2
          figComH = [];
          figFpeH = [];
          figCapH = [];
          figFpeCapH = [];
          figGrfZH = [];
          figBPH = [];

          figComTitle = [];
          figFpeTitle = [];
          figCapTitle = [];
          figFpeCapTitle = [];
          figGrfZTitle = [];
          figBPTitle = [];

          [row,col] = find(subPlotPanelIndex == (indexTrial-offsetTrialIndex));          
          subplotPosition = reshape(subPlotPanel(row,col,:),1,4);
          
          flag_zeroAtSeatOff = 1;
          
          if(figType==1)
            if(flag_balancePortraitSubject == 1)
              figBPH = figSubjectBalancePortrait(indexSubject).h;
              analysisType = '';
              switch flag_ModeBalancePortrait
                case 0
                  analysisType = 'Fpe BP ';
                case 1
                  analysisType = 'Cap BP ';
              end
              figBPTitle = [analysisType,'(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];

            end

            if(flag_comKinematicsTimeSeriesSubject==1)
              figComH = figSubjectCom(indexSubject).h;
              analysisType = '';
              switch flag_ModeComVel0ComGpVsCop1
                case 0
                  analysisType = 'Com-Vel';
                case 1
                  analysisType = 'Com-Cop';
              end
              figComTitle = [analysisType,'(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end
            
            if(flag_fpeTimeSeriesPlotsSubject==1)
              figFpeH = figSubjectFpe(indexSubject).h;
              
              figTitle = '';
              switch flag_ModeBalancePointsVsCom0VsCop1
                case 0
                  figTitle = 'Fpe-Com ';
                case 1
                  figTitle = 'Fpe-Cop ';
              end
              
              figFpeTitle = [figTitle,'(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_capTimeSeriesPlotsSubject==1)              
              figCapH = figSubjectCap(indexSubject).h;
              figTitle = '';
              switch flag_ModeBalancePointsVsCom0VsCop1
                case 0
                  figTitle = 'Cap-Com ';
                case 1
                  figTitle = 'Cap-Cop ';
              end

              figCapTitle = [figTitle,'(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end  
            
            if(flag_fpeCapErrorTimeSeriesPlotsSubject==1)
              figFpeCapH = figSubjectFpeCapErr(indexSubject).h;
              figFpeCapTitle = ['Fpe-Cap Err. ','(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_GrfzTimeSeriesPlotsSubject==1 )
              figGrfZH = figSubjectGrfZ(indexSubject).h;
              figGrfZTitle = ['GRFz ','(',subjectId,'): ',...
                               trialTypeNames{indexTrial}];
            end
            
          else
            if(flag_balancePortraitEnsemble == 1)
              figBPH = figAllBalancePortrait;
              analysisType = '';
              switch flag_ModeBalancePortrait
                case 0
                  analysisType = 'Fpe BP: ';
                case 1
                  analysisType = 'Cap BP: ';
              end
              figBPTitle = [analysisType, ...
                             trialTypeNames{indexTrial}];
            end

            if(flag_comKinematicsTimeSeriesEnsemble==1)
              figComH = figAllCom;
              analysisType = '';
              switch flag_ModeComVel0ComGpVsCop1
                case 0
                  analysisType = 'Com-Vel';
                case 1
                  analysisType = 'Com-Cop';
              end
              figComTitle = [analysisType,':', trialTypeNames{indexTrial}];
            end            
            if(flag_fpeTimeSeriesPlotsEnsemble==1)
              figFpeH = figAllFpe;
              switch flag_ModeBalancePointsVsCom0VsCop1
                case 0
                  figTitle = 'Fpe-Com: ';
                case 1
                  figTitle = 'Fpe-Cop: ';
              end
              figFpeTitle = [figTitle,...
                             trialTypeNames{indexTrial}];
            end  
            
            if(flag_capTimeSeriesPlotsEnsemble==1)
              figCapH = figAllCap;
              figTitle = '';
              switch flag_ModeBalancePointsVsCom0VsCop1
                case 0
                  figTitle = 'Cap-Com: ';
                case 1
                  figTitle = 'Cap-Cop: ';
              end
              figCapTitle = [figTitle,...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_fpeCapErrorTimeSeriesPlotsEnsemble==1)
              figFpeCapH = figAllFpeCapErr;
              figFpeCapTitle = ['Fpe-Cap Err.: ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_GrfzTimeSeriesPlotsEnsemble==1 )
              figGrfZH = figAllGrfZ;
              figGrfZTitle = ['GRFz: ',...
                               trialTypeNames{indexTrial}];
            end
            
          end
          
          if(isempty(figBPTitle)==0)
            balancePoint = [];
            balanceStepDir = [];
            lineColor = [];
            switch flag_ModeBalancePortrait
              case 0
                balancePoint = fpeData.r0F0;
                balanceStepDir = zeros(size(balancePoint));
                eK = -gravityVec./norm(gravityVec);
                for k=1:1:size(balancePoint,1)
                  eN = fpeData.n(k,:)';
                  balanceStepDir(k,:) = (getCrossProductMatrix(eN)*eK)';
                end
                lineColor = colorFpe;
              case 1
                balancePoint = capData.r0F0;
                balanceStepDir = capData.u;
                lineColor = colorCap;
            end
            c3dFootMarkersLeftNames = {'L_FAL','L_FM1','L_FM2','L_FM5','L_TAM','L_FCC'};
            c3dFootMarkersRightNames = {'R_FAL','R_FM1','R_FM2','R_FM5','R_TAM','R_FCC'};

            allowedFootMarkerMovement = 0.02;
            
            lineColorBalancePoint = [1,1,1/0.25].*0.25;
            lineColorBalanceDir =[1,1,1/0.75].*0.75;
            lineColorCop     = [1,0,0];
            lineColorCom     = [0.5,0.5,0.5];
            footPatchColor   = [1,1,1].*0.75;
            

            figBPH = plotBalancePortrait(figBPH, subplotPosition,...
                      indexSubject, subjectId,...
                      c3dTime, c3dMarkers, ...
                      c3dFootMarkersRightNames, c3dFootMarkersLeftNames, ...
                      sitToStandQuietSequence, segmentedData,...
                      c3dGrf(index_ChairForcePlate), ...
                      c3dGrf(index_FeetForcePlate),...
                      wholeBodyData(:,colComPos), ...
                      balanceStepDir, balancePoint,...
                      indexSittingDynamic, indexCrouchingStable,...
                      lineColorBalanceDir,lineColorBalancePoint,...
                      lineColorCop,lineColorCom, ...
                      footPatchColor,figBPTitle, ...
                      flag_zeroAtSeatOff,...
                      allowedFootMarkerMovement);
                  
            here=1;
          end

          if(isempty(figComTitle)==0)
            figComH = plotComKinematicsTimeSeries(...
                   figComH, subplotPosition, ...
                   c3dTime, sitToStandQuietSequence, ...
                   c3dGrf(index_ChairForcePlate),...
                   c3dGrf(index_FeetForcePlate),...
                   wholeBodyData(:,colComPos),...
                   wholeBodyData(:,colComVel),...
                   gravityVec,...
                   indexSittingDynamic, indexCrouchingStable,...
                   colorCom,...
                   figComTitle,...
                   flag_zeroAtSeatOff,...
                   flag_ModeComVel0ComGpVsCop1);
          end
          
          if(flag_comKinematicsSeatOffEnsemble==1)
            flag_pointInTimeSelector = 1;
            figTitle = '';
            switch flag_ModeComVel0ComGpVsCop1
              case 0
                figTitle = ['Com Vel At Seat-Off: ',...
                              trialTypeNames{indexTrial}];
              case 1
                figTitle = ['Com-Cop Dist. At Seat-Off: ',...
                              trialTypeNames{indexTrial}];
            end
            figTrialSeatOffCom = ...
              plotComKinematicsAtPointInTime(...
                   figTrialSeatOffCom, subplotPosition, ...
                   indexSubject, subjectId,...
                   c3dTime, sitToStandQuietSequence, ...
                   c3dGrf(index_ChairForcePlate),...
                   c3dGrf(index_FeetForcePlate),...
                   wholeBodyData(:,colComPos),...
                   wholeBodyData(:,colComVel),...
                   gravityVec,...
                   indexSittingDynamic, indexCrouchingStable,...
                   colorCom,...
                   figTitle,...
                   flag_pointInTimeSelector,...
                   flag_ModeComVel0ComGpVsCop1);

          end
          
          
          if( isempty(figFpeTitle)==0)              
              figFpeH = ...
                plotFpeStepLengthTimeSeries(...
                     figFpeH, subplotPosition,... 
                     c3dTime, sitToStandQuietSequence, ...
                     c3dGrf(index_ChairForcePlate),...
                     c3dGrf(index_FeetForcePlate), fpeData, gravityVec,...
                     indexSittingDynamic, indexCrouchingStable, colorFpe,...
                     figFpeTitle,...
                     flag_zeroAtSeatOff,...
                     flag_ModeBalancePointsVsCom0VsCop1,...
                     flag_ModeAnalyzeBalanceAlong0Across1);
          end
          
          
          
          if(flag_fpeSeatOffEnsemble==1)
            figTitle = '';
            switch flag_ModeBalancePointsVsCom0VsCop1
              case 0
                figTitle = 'Fpe-Com ';
              case 1
                figTitle = 'Fpe-Cop ';
            end
            switch flag_ModeAnalyzeBalanceAlong0Across1
              case 0
                figTitle = [figTitle,'at Seat-Off: '];
              case 1
                figTitle = [figTitle,'at Seat-Off: '];
            end
            flag_pointInTimeSelector = 1;
            figTrialSeatOffFpe = plotFpeStepLengthAtPointInTime(...
                    figTrialSeatOffFpe, subplotPosition,... 
                   indexSubject, subjectId, sitToStandQuietSequence, ...
                   c3dGrf(index_ChairForcePlate),...
                   c3dGrf(index_FeetForcePlate),fpeData,gravityVec,...
                   indexSittingDynamic, indexCrouchingStable,...
                   colorFpe,...
                   [figTitle,trialTypeNames{indexTrial}],...
                   flag_pointInTimeSelector,...
                   flag_ModeBalancePointsVsCom0VsCop1,...
                     flag_ModeAnalyzeBalanceAlong0Across1);
            axis tight;
            axisLim = axis;            
            axisLim(1) = 0;
            axisLim(2) = length(subjectsToProcess)+1;
            if(flag_ModeBalancePointsVsCom0VsCop1 == 0)
              axisLim(3) = -1;
              axisLim(4) = 20;
            else
              if(flag_ModeAnalyzeBalanceAlong0Across1==0)
                axisLim(3) = -14;
                axisLim(4) = 8;
              else
                axisLim(3) = -5;
                axisLim(4) = 5;
              end
            end
            axis(axisLim);
            set(gca,'XTickLabel',[]);
          end
          
          if(flag_capSeatOffEnsemble==1)
            figTitle = '';
            switch flag_ModeBalancePointsVsCom0VsCop1
              case 0
                figTitle = 'Cap-Com ';
              case 1
                figTitle = 'Cap-Cop ';
            end
            switch flag_ModeAnalyzeBalanceAlong0Across1
              case 0
                figTitle = [figTitle,'at Seat-Off: '];
              case 1
                figTitle = [figTitle,'at Seat-Off: '];
            end
            flag_pointInTimeSelector = 1;
            figTrialSeatOffCap = plotCapStepLengthAtPointInTime(...
                    figTrialSeatOffCap, subplotPosition,... 
                   indexSubject, subjectId, sitToStandQuietSequence,...
                   c3dGrf(index_ChairForcePlate),...
                   c3dGrf(index_FeetForcePlate),capData,...
                   indexSittingDynamic, indexCrouchingStable,...
                   colorCap,...
                   [figTitle,trialTypeNames{indexTrial}],...
                   flag_pointInTimeSelector,...
                   flag_ModeBalancePointsVsCom0VsCop1,...
                   flag_ModeAnalyzeBalanceAlong0Across1);
            axis tight;
            axisLim = axis;
            axisLim(1) = 0;
            axisLim(2) = length(subjectsToProcess)+1;
            if(flag_ModeBalancePointsVsCom0VsCop1 == 0)
              axisLim(3) = -1;
              axisLim(4) = 20;
            else
              if(flag_ModeAnalyzeBalanceAlong0Across1==0)
                axisLim(3) = -14;
                axisLim(4) = 8;
              else
                axisLim(3) = -5;
                axisLim(4) = 5;
              end
            end
            axis(axisLim);
            set(gca,'XTickLabel',[]);
          end          
          
          if( isempty(figCapTitle)==0)              
              figCapH = ...
                plotCapStepLengthTimeSeries(...
                     figCapH, subplotPosition,... 
                     c3dTime, sitToStandQuietSequence, ...
                     c3dGrf(index_ChairForcePlate),...
                     c3dGrf(index_FeetForcePlate),capData,...
                     indexSittingDynamic, indexCrouchingStable, colorCap,...
                     figCapTitle,...
                     flag_zeroAtSeatOff,...
                     flag_ModeBalancePointsVsCom0VsCop1,...
                   flag_ModeAnalyzeBalanceAlong0Across1);
          end
          
          if(isempty(figFpeCapTitle)==0)
            figFpeCapH = ...
              plotFpeCapErrorTimeSeries(...
                   figFpeCapH, subplotPosition,... 
                   c3dTime, sitToStandQuietSequence, ...
                   c3dGrf(index_ChairForcePlate),...
                   fpeData, capData,...
                   indexSittingDynamic, indexCrouchingStable, colorFpeCap,...
                   figFpeCapTitle,...
                   flag_zeroAtSeatOff);
          end 
          
          if(isempty(figGrfZH)==0)
            figGrfZH = ...
                plotGrfZTimeSeries(...
                     figGrfZH, subplotPosition,... 
                     c3dTime, sitToStandQuietSequence, ...
                     c3dGrf(index_ChairForcePlate),...
                     indexSittingDynamic, indexCrouchingStable, colorGrfz,...
                     figGrfZTitle,...
                     flag_zeroAtSeatOff);  
          end
          
          if(flag_GrfzSeatOffEnsemble==1)
            flag_pointInTimeSelector = 1;
            figTrialGrfZSeatOff = ...
                plotGrfZAtPointInTime(...
                     figTrialGrfZSeatOff, subplotPosition,... 
                     indexSubject, subjectId, c3dTime, ...
                     sitToStandQuietSequence, c3dGrf(index_ChairForcePlate),...
                     indexSittingDynamic, indexCrouchingStable, colorGrfz,...
                     figGrfZTitle,...
                     flag_pointInTimeSelector); 
          end

        end
      end
      
      
    end
    
  end
  
  if(flag_balancePortraitSubject==1)
    figure(figSubjectBalancePortrait(indexSubject).h);
    configPlotExporter;
    analysisType='';
    switch flag_ModeBalancePortrait      
      case 0
        analysisType = ['FpeBalancePortrait'];
      case 1
        analysisType = ['CapBalancePortrait'];        
    end
    print('-dpdf',['../../plots/balancePortrait/fig_',analysisType,'_S2SQuietTimeSeries_',subjectId,'.pdf']);  
    
  end


  if(flag_comKinematicsTimeSeriesSubject==1)
    figure(figSubjectCom(indexSubject).h);
    configPlotExporter;    
    analysisType = '';
    switch flag_ModeComVel0ComGpVsCop1
      case 0
        analysisType = ['ComVel'];
      case 1
        analysisType = ['ComVsCop'];
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
    print('-dpdf',['../../plots/fpeCap/fig_FpeCapErrS2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  
  if(flag_GrfzTimeSeriesPlotsSubject==1)
    figure(figSubjectGrfZ(indexSubject).h);
    configPlotExporter;
    print('-dpdf',['../../plots/grfz/fig_GrfZS2SQuietTimeSeries_',subjectId,'.pdf']);
  end
end


if(flag_comKinematicsTimeSeriesEnsemble==1)  
  figure(figAllCom);
  configPlotExporter;    
  analysisType = '';
  switch flag_ModeComVel0ComGpVsCop1
    case 0
      analysisType = ['ComVel'];
    case 1
      analysisType = ['ComVsCop'];
  end               
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);
end 

if(flag_comKinematicsSeatOffEnsemble==1)
  figure(figTrialSeatOffCom)
  configPlotExporter;    
  analysisType = '';
  switch flag_ModeComVel0ComGpVsCop1
    case 0
      analysisType = ['ComVel'];
    case 1
      analysisType = ['ComVsCop'];
  end               
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);
end

if(flag_fpeTimeSeriesPlotsEnsemble==1)
  figure(figAllFpe);
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
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);  
end

if(flag_fpeSeatOffEnsemble)
  figure(figTrialSeatOffFpe);    
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
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);    
  
end
if(flag_capSeatOffEnsemble)
  figure(figTrialSeatOffCap);
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
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietSeatOff.pdf']);  
end

if(flag_capTimeSeriesPlotsEnsemble==1)
  figure(figAllCap);
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
  print('-dpdf',['../../plots/fig_',analysisType,'_S2SQuietTimeSeries.pdf']);  
end

if(flag_fpeCapErrorTimeSeriesPlotsEnsemble==1)
  figure(figAllFpeCapErr);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_FpeCapErrS2SQuietTimeSeries','.pdf']);
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
