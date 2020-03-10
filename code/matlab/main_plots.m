clc;
close all;
clear all;

subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE05','configE08'};



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

figSubject(length(subjectsToProcess)) = struct('h',[]);
figAll = figure;

colorElderly = [0,0,0];
colorYoung = [0.5,0.5,1];


for indexSubject = 1:1:length(subjectsToProcess)
  figSubject(indexSubject).h = figure;
  
  colorTrial = [1,0,1];
  if( isempty(strfind(subjectsToProcess{indexSubject},'E'))==0)
    colorTrial = colorElderly;
  else
    colorTrial = colorYoung;    
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
          if(figType==1)
            figure(figSubject(indexSubject).h);          
          else
            figure(figAll);
          end        
          subplot(2,3,indexTrial-1);      
          for z=1:1:length(sitToStandSequence)
            if( sum(isnan(sitToStandSequence(z).phaseTransitions))==0)

              idx0 = sitToStandSequence(z).phaseTransitionIndices(1);

              idx1 = idx0;
              numberOfTransitions = ...
                  size(sitToStandSequence(z).phaseTransitionIndices,1);

              for indexPhase=1:1:numberOfTransitions
                if(sitToStandSequence(z).phaseTransitions(indexPhase,1)...
                    ==indexSittingDynamic)
                  idx1 = sitToStandSequence(z).phaseTransitionIndices(indexPhase);
                end
              end


              idx2 = sitToStandSequence(z).phaseTransitionIndices(end);

              timeBias = c3dTime(idx1);

              rGF0 = fpeData.r0F0(idx0:1:idx2,:) ...
                    -fpeData.r0G0(idx0:1:idx2,:);
              distanceFromStability = sum(rGF0.^2,2).^0.5;

              plot( c3dTime(idx0:1:idx2,1)-timeBias, distanceFromStability.*100,...
                    'Color',colorTrial);
              hold on;
              xlabel('Time (s)');
              ylabel('Distance from Stability (cm)');
              box off;
              %subjectId
              if(figType==1)
                title([subjectId,': ',trialTypeNames{indexTrial}]);        
              else
                title(['All: ',trialTypeNames{indexTrial}]);
              end           

            end

          end

        end
      end
    end
    
  end
  
  figure(figSubject(indexSubject).h);
  saveas(gcf,['../../plots/fig_FpeLenthAtSeatOff_',subjectId],'pdf');
  
end
  figure(figAll);
  saveas(gcf,['../../plots/fig_FpeLenthAtSeatOff'],'pdf');

