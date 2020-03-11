clc;
close all;
clear all;

flag_fpeTimeSeriesPlotsSubject   = 1;
flag_fpeTimeSeriesPlotsEnsemble  = 1;
flag_fpeSeatOffEnsemble          = 1;

flag_capTimeSeriesPlotsSubject   = 1;
flag_capTimeSeriesPlotsEnsemble  = 1;
%flag_capSeatOffEnsemble

%flag_comVelTimeSeriesPlotsSubject   = 0;
%flag_comVelTimeSeriesPlotsEnsemble  = 0;
%flag_comVelSeatOffEnsemble

%flag_comVelTimeSeriesPlotsSubject   = 0;
%flag_comVelTimeSeriesPlotsEnsemble  = 0;
%flag_comVelSeatOffEnsemble

flag_fpeCapErrorTimeSeriesPlotsSubject   = 1;
flag_fpeCapErrorTimeSeriesPlotsEnsemble  = 1;
%flag_fpeCapErrorSeatOffEnsemble

flag_GrfzTimeSeriesPlotsSubject  = 1;
flag_GrfzTimeSeriesPlotsEnsemble = 1;
%flag_GrfzSeatOffEnsemble


plotConfig;

subjectsToProcess =  ...
  {'configH01','configH02','configH03','configH04','configH05',...
   'configH06','configH07','configH08','configH09','configH10',...
   'configE01','configE02','configE03','configE05', ...
   'configE07','configE08'};

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


colorElderly = [0,0,0;...
                1,0,0;...
                1,0,1;...
                0,0,1];
colorYoung = [0.5,0.5,0.5;...
              1.0,0.5,0.5;...
              1.0,0.5,1.0;...
              0.5,0.5,1.0];


for indexSubject = 1:1:length(subjectsToProcess)
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
  %figSubjectFpeCapErr(indexSubject).h = figure;
  %figSubjectComVel(   indexSubject).h = figure;
  %figSubjectComCop(   indexSubject).h = figure;
  
  colorFpe    = [0,0,0];
  colorFpeCap = [0,0,0];
  colorCap    = [0,0,0];
  colorComVel = [0,0,0];
  colorComCop = [0,0,0];
  colorGrfz   = [0,0,0];
  
  
  
  if( isempty(strfind(subjectsToProcess{indexSubject},'E'))==0)
    colorFpe    = colorElderly(1,:);
    colorFpeCap = colorElderly(1,:);  
    colorCap    = colorElderly(1,:);    
    colorComVel = colorElderly(1,:);
    colorComCop = colorElderly(1,:);  
    colorGrfz   = colorElderly(1,:);
  else
    colorFpe    = colorYoung(1,:);    
    colorFpeCap = colorYoung(1,:);  
    colorCap    = colorYoung(1,:);    
    colorComVel = colorYoung(1,:);
    colorComCop = colorYoung(1,:);    
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
          figFpeH = [];
          figCapH = [];
          figFpeCapH = [];
          figGrfZH = [];
          
          figFpeTitle = [];
          figCapTitle = [];
          figFpeCapTitle = [];
          figGrfZTitle = [];
          
          [row,col] = find(subPlotPanelIndex == (indexTrial-offsetTrialIndex));          
          subplotPosition = reshape(subPlotPanel(row,col,:),1,4);
          
          flag_zeroAtSeatOff = 1;
          
          if(figType==1)
            
            if(flag_fpeTimeSeriesPlotsSubject==1)
              figFpeH = figSubjectFpe(indexSubject).h;
              figFpeTitle = ['FPE ','(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_capTimeSeriesPlotsSubject==1)
              figCapH = figSubjectCap(indexSubject).h;
              figCapTitle = ['Cap ','(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end  
            
            if(flag_fpeCapErrorTimeSeriesPlotsSubject==1)
              figFpeCapH = figSubjectFpeCapErr(indexSubject).h;
              figFpeCapTitle = ['FPE-Cap Err. ','(',subjectId,'): ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_GrfzTimeSeriesPlotsSubject==1 )
              figGrfZH = figSubjectGrfZ(indexSubject).h;
              figGrfZTitle = ['GRFz ','(',subjectId,'): ',...
                               trialTypeNames{indexTrial}];
            end
            
          else
            
            if(flag_fpeTimeSeriesPlotsEnsemble==1)
              figFpeH = figAllFpe;
              figFpeTitle = ['FPE: ',...
                             trialTypeNames{indexTrial}];
            end  
            
            if(flag_capTimeSeriesPlotsEnsemble==1)
              figCapH = figAllCap;
              figCapTitle = ['Cap: ',...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_fpeCapErrorTimeSeriesPlotsEnsemble==1)
              figFpeCapH = figAllFpeCapErr;
              figFpeCapTitle = ['FPE-Cap Err.: ', ...
                             trialTypeNames{indexTrial}];
            end 
            
            if(flag_GrfzTimeSeriesPlotsEnsemble==1 )
              figGrfZH = figAllGrfZ;
              figGrfZTitle = ['GRFz: ',...
                               trialTypeNames{indexTrial}];
            end
            
          end
          
          if( isempty(figFpeTitle)==0)              
              figFpeH = ...
                plotFpeStepLengthTimeSeries(...
                     figFpeH, subplotPosition,... %[2,3,indexTrial-offsetTrialIndex], ...
                     c3dTime, sitToStandQuietSequence, fpeData,...
                     indexSittingDynamic, indexCrouchingStable, colorFpe,...
                     figFpeTitle,...
                     flag_zeroAtSeatOff);
          end
             
          if(flag_fpeSeatOffEnsemble==1)
            flag_pointInTimeSelector = 1;
            figTrialSeatOffFpe = plotFpeStepLengthAtPointInTime(...
                    figTrialSeatOffFpe, subplotPosition,... %[2,3,indexTrial-offsetTrialIndex], ...
                   indexSubject, subjectId, sitToStandQuietSequence, fpeData,...
                   indexSittingDynamic, indexCrouchingStable,...
                   colorFpe,...
                   ['FPE-Length At Seat-Off ',trialTypeNames{indexTrial}],...
                   flag_pointInTimeSelector);
            axis tight;
            axisLim = axis;
            axisLim(1) = 0;
            axisLim(2) = length(subjectsToProcess)+1;
            axisLim(3) = -1;
            axisLim(4) = 20;
            axis(axisLim);
            set(gca,'XTickLabel',[]);
          end
          
          if( isempty(figCapTitle)==0)              
              figCapH = ...
                plotCapStepLengthTimeSeries(...
                     figCapH, subplotPosition,... %[2,3,indexTrial-offsetTrialIndex], ...
                     c3dTime, sitToStandQuietSequence, capData,...
                     indexSittingDynamic, indexCrouchingStable, colorFpe,...
                     figFpeTitle,...
                     flag_zeroAtSeatOff);
          end
          
          if(isempty(figFpeCapTitle)==0)
            figFpeCapH = ...
              plotFpeCapErrorTimeSeries(...
                   figFpeCapH, subplotPosition,... %[2,3,indexTrial-offsetTrialIndex],  ...
                   c3dTime, sitToStandQuietSequence, fpeData, capData,...
                   indexSittingDynamic, indexCrouchingStable, colorFpeCap,...
                   figFpeCapTitle,...
                   flag_zeroAtSeatOff);
          end 
          
          if(isempty(figGrfZH)==0)
            figGrfZH = ...
                plotGrfZTimeSeries(...
                     figGrfZH, subplotPosition,... %[2,3,indexTrial-offsetTrialIndex], ...
                     c3dTime, sitToStandQuietSequence, c3dGrf(index_ChairForcePlate),...
                     indexSittingDynamic, indexCrouchingStable, colorFpe,...
                     figGrfZTitle,...
                     flag_zeroAtSeatOff);  
          end
            

        end
      end
      
      
    end
    
  end
  

  
  if(flag_fpeTimeSeriesPlotsSubject==1)
    figure(figSubjectFpe(indexSubject).h);
    configPlotExporter;
    print('-dpdf',['../../plots/fpe/fig_FpeS2SQuietTimeSeries_',subjectId,'.pdf']);
  end
  if(flag_capTimeSeriesPlotsSubject==1)
    figure(figSubjectCap(indexSubject).h);
    configPlotExporter;    
    print('-dpdf',['../../plots/cap/fig_CapS2SQuietTimeSeries_',subjectId,'.pdf']);
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

if(flag_fpeTimeSeriesPlotsEnsemble==1)
  figure(figAllFpe);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_FpeS2SQuietTimeSeries','.pdf']);
end
if(flag_fpeSeatOffEnsemble)
  figure(figTrialSeatOffFpe);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_FpeS2SQuietSeatOff','.pdf']);
  
end

if(flag_capTimeSeriesPlotsEnsemble==1)
  figure(figAllCap);
  configPlotExporter;
  print('-dpdf',['../../plots/fig_capS2SQuietTimeSeries','.pdf']);
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

