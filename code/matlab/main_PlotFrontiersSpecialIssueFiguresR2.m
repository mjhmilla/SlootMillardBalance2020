
clc;
close all;
clear all;

%%
%
% Constants
%
%%

flag_PaperPlots0Presentation1 = 0;

flag_useBasicBosModel = 0;
bosModelTag = '';
if(flag_useBasicBosModel == 1)
  bosModelTag = '_SimpleBos';  
end

flag_ConvexHullWithToes0Without1  = 0;
nameToeTag = '';
if(flag_ConvexHullWithToes0Without1 == 1)
  nameToeTag = '_NoToes';
end


flag_excludeE06 = 0;
nameExcludeE06Tag = '';
if(flag_excludeE06 == 1)
  nameExcludeE06Tag = '_NoE06';
end

numberOfPhases = 2;
indexPhaseToAnalyze=NaN;
phaseNames = {''};

flag_motionSequence = 0; %0. Sit to stand
                         %1. Stand to sit
                         
flag_motionSegmentEmphasis = NaN;
%0. Start
%1. Whole-movement
%2. End
                         
motionSequenceTag = '';
outputPlotFolder  = '';
outputDataFolder  = '';
phaseLabels = {};
switch (flag_motionSequence)
  case 0
    motionSequenceTag = '_sit2Stand';
    outputPlotFolder='../../plots/Frontiers2020Pub/';
    outputDataFolder='../../outputData/Frontiers2020Pub/';
    phaseNames = {'Phase1Sit2SeatOff','Phase2SeatOff2Stand'};
    numberOfPhases = length(phaseNames);
    indexPhaseToAnalyze=2;
    flag_motionSegmentEmphasis = 0;
    phaseLabels = {'seat-off','sit-to-stand','standing'};
  case 1
    motionSequenceTag = '_stand2Sit';    
    outputPlotFolder='../../plots/StandToSitPub/';  
    outputDataFolder='../../outputData/StandToSitPub/'; 
    phaseNames = {'Phase1Stand2SeatOn','Phase2SeatOn2Sit'};
    numberOfPhases = length(phaseNames);
    indexPhaseToAnalyze=1;
    flag_motionSegmentEmphasis = 2;
    phaseLabels = {'standing','stand-to-sit','seat-on'};
  otherwise
    assert(0);
end


% 1st phase: start to seatoff
% 2nd phase: seatoff to standing;

%%
%
% Plot Configuration
%
%%
flagPlotStart =1;
flagPlotPhase =1;
flagPlotEnd   =1;
flag_useCustomYTicks = 0;

indexPhase = indexPhaseToAnalyze;

trialsToProcess = {'Side','Chest','Conv','Leg','Side','Rob'};

quietStandingFile = 'QuietStanding_SubAnalysis/quietStanding_NoE06.mat';


%%
%
% Setup the input/output directory structure
%
%%
inputDirRelative    = '../../inputData'         ;
outputDirRelative   = '../../outputData'        ;
%frontiersPlotDir    = '../../plots/Frontiers2020Pub/' ;
%frontiersDataDir    = '../../outputData/Frontiers2020Pub/';

codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  processConfigCommon;  
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd; 
  quietStandingFilePath = [outputPath,'/',quietStandingFile];
cd(codeDir);

%%
%
% Load the data
%
%%

tmp = load(quietStandingFilePath);
quietStandingData = tmp.dataStruct;

load([outputDataFolder,'groupTrialPhaseMetricData',...
      bosModelTag,nameToeTag,nameExcludeE06Tag,...
      motionSequenceTag,'.mat']);
load([outputDataFolder,'subjectTrialPhaseMetricData',...
      bosModelTag,nameToeTag,nameExcludeE06Tag,...
      motionSequenceTag,'.mat']);

fpeCopInTOneSided = mean( [mean(abs(quietStandingData.fpeCopInT(3,:))), ...
                          mean(abs(quietStandingData.fpeCopInT(2,:)))]);                             
fpeCopInSOneSided = mean( [mean(abs(quietStandingData.fpeCopInS(3,:))), ...
                          mean(abs(quietStandingData.fpeCopInS(2,:)))] );

angZSpeedOneSided = mean( [abs(mean(quietStandingData.angZSpeed(3,:))),...
                          abs(mean(quietStandingData.angZSpeed(2,:)))]);   
    
%%
% Configure the plots
%%

flag_plotConfig = 0; 
% 0. Frontiers balance metrics
% 1. Com velocity
% 2. Com angular velocity
% 3. Switching Conditions: Time since peak (fpe-com).u, fz err at seat-off
% 4. Summary plot
plotConfigName = '';

flag_plotSubjectData = 0;
flag_plotGroupData   = 0;


switch(flag_PaperPlots0Presentation1)
  case 0
    lineWidth  = 0.75;
    boxWidth   = 0.33;
    panelHeight = 5;
    panelWidth  = 16;
  case 1
    lineWidth  = 0.5;
    boxWidth   = 0.33;
    panelWidth = 12.8;
    panelHeight  = 9.6;
  otherwise assert(0);
end    


metricYTick = [];
metricOrder = [];

metricUpperLeftNotes = {};
metricLowerLeftNotes = {};

switch flag_plotConfig
  case 0
    
    if(flag_PaperPlots0Presentation1)
      panelHeight = panelHeight*0.5;
    end
        
    metricOrder = [1,2,3,7, 8,4,5,6];
    
    metricNameList = {'com2edge','com2cop','comvel',...
                      'fpe2edge','fpewidth','fpelen','angvel','angvelz'};

    %metricYLim = [-9,12; 0,15.1; 0,65;...
    %              -2.,12.0; -3.6,3.6; -4,11;...
    %              0,102; -35, 35];
    metricYLim = [ -11, 11;...
                   0.0, 15.1; ...
                   0.0, 65;...
                  -2.0, 12.0; ...
                  -3.6,  3.6;...
                    -8, 10.1;...
                     0,102;...
                   -35, 35];

         
    if(flag_motionSequence==1)
      metricYLim = [ -11, 11;...
                     0.0, 15.1; ...
                     0.0, 65;...
                    -20.0, 10.0; ...
                    -3.6,  3.6;...
                      -5.1, 30.;...
                       0,130;...
                     -40.5, 15];
    end                
                
                
    metricYTick = [5;3;10; ...
                   3;1;3;...
                   20;10];
                
    metricYAxisLabel = {'Displacement (cm)','Distance (cm)','Speed (cm/s)',...
                        'Displacement (cm)','Displacement (cm)','Displacement (cm)',...
                        'Angular Speed ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)'};

    metricYAxisScale = [1,1,1,1,...
                         1,1,1,1];


                       
    metricTitle = {'A. Displacement from nearest BOS edge to COM$_{\mathrm{GP}}$ ($_{B(C)}r{_C}$)',...
                   'B. Distance from COM$_{\mathrm{GP}}$ to COP ($|{_C}r{_P}|_2$)',...
                   'C. COM Speed ($|v{_C}|_2$)',...
                   'B. Displacement from BOS edge to FPE ($_{B(F)}r{_F}$)',...
                   'C. Displacement from COP to FPE in $\hat{t}$ (${_P}r{_F} \cdot \hat{t}$)',...
                   'D. Displacement from COP to FPE in $\hat{s}$ (${_P}r{_F} \cdot \hat{s}$)',...
                   'D. Angular speed ($|\omega|_2$)',...
                   'A. Angular velocity about the vertical axis ($\omega_Z$)'};

    metricPlotHalfPlaneBox = [1,1,1,1,1,1,1,1];   
    
    metricPlotZeroLine = [0,0,0,0,1,1,0,1];
    
      
    
                   
    metricPlotHalfPlaneYmax = [metricYLim(1,2),...
                               mean(quietStandingData.comCopDist(3,:)),...
                               mean(quietStandingData.comSpeed(3,:)),...
                               metricYLim(4,2),...
                               fpeCopInTOneSided,...
                               fpeCopInSOneSided,...
                               mean(quietStandingData.angSpeed(3,:)),...
                               angZSpeedOneSided];
                             
    metricPlotHalfPlaneYmin = [0,...
                               0,...
                               0,...
                               0,...
                               -fpeCopInTOneSided,...
                               -fpeCopInSOneSided,...
                               0,...
                               -angZSpeedOneSided];
                            
                             
    metricUpperLeftNotes = {'(+) Inside BOS','','',...
                            '(+) Inside BOS','(+) Left turn','(+) Accelerate Forwards',...
                            '','(+) CCW'};
    metricLowerLeftNotes = {'(-) Outside BOS','','',...
                            '(-) Outside BOS','(-) Right turn','(-) Accelerate Backwards',...
                            '','(-) CW'};

                            
    metricPlotHalfPlaneText = ...
      {'Stat. Balanced:',...
         '${_C}r{_B} > 0\,\mathrm{cm}$';...
       'Quiet standing:', ...
          sprintf('$|{_C}r{_P}|_2 < %1.1f\\,\\mathrm{cm}$',mean(quietStandingData.comCopDist(3,:)));...
       'Quiet standing:', ...
          sprintf('$|v{_C}|_2 < %1.1f\\,\\mathrm{cm}/\\mathrm{s}$',...
                 mean(quietStandingData.comSpeed(3,:)));...
       'Dyn. balanced','standing ${_F}r{_B} > 0$';...
       'Quiet standing:',sprintf('$|{_F}r{_C} \\cdot \\hat{t}| < %1.1f\\,\\mathrm{cm}$', fpeCopInTOneSided);...
       'Quiet standing:',sprintf('$|{_F}r{_C} \\cdot \\hat{s}| < %1.1f\\,\\mathrm{cm}$', fpeCopInSOneSided);...
       'Quiet standing:',sprintf('$|\\omega|_2 < %1.1f\\,^\\circ/s$',mean(quietStandingData.angSpeed(3,:)));...
       'Quiet standing:',sprintf('$|\\omega_Z| < %1.1f\\,^\\circ/s$',angZSpeedOneSided),...
      };

    metricAnnotationLines = [0,0,0,0,0,0,0,0];
    
    plotConfigName = 'ComFpeBalanceMetrics';   
    
    flag_plotSubjectData = 1;
    flag_plotGroupData   = 1;

    % Plot Config    

    numberOfFiguresPerPage        = length(metricNameList);
    numberOfVerticalPlotRows      = numberOfFiguresPerPage;
    numberOfHorizontalPlotColumns = 1;    
    assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
             >= numberOfFiguresPerPage);

    plotHorizMarginCm = 1.5;
    plotVertMarginCm  = 1.5;           
    pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
    pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;           
    plotHeight  = panelHeight;
    plotWidth   = panelWidth ;
    
    plotConfigGeneric;
 
  case 1
    metricNameList = {'comvel','comvelx','comvely','comvelz'};

    metricYLim = [ 0,65; -20,65; -10,10; -40,65];
         
    if(flag_motionSequence==1)
      metricYLim(2,:) =[-50,5];
      metricYLim(4,:) =[-50,5];
    end
    
    metricYAxisLabel = {'Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)'};

    metricYAxisScale = [1,1,1,1];

    metricTitle = {'COM Speed',...
                   'COM Velocity in X',...
                   'COM Velocity in Y',...
                   'COM Velocity in Z'};

    metricPlotHalfPlaneBox = [0,0,0,0];             
    metricPlotHalfPlaneYmax = [0,0,0,0];                             
    metricPlotHalfPlaneYmin = metricYLim(:,1);    
    metricPlotZeroLine = [0,0,0,0];
    
    metricAnnotationLines = [0,0,0,0];

    plotConfigName = 'ComVelocity';    

    flag_plotSubjectData = 1;
    flag_plotGroupData   = 1;
    
    % Plot Config    
    numberOfFiguresPerPage        = length(metricNameList);
    numberOfVerticalPlotRows      = numberOfFiguresPerPage;
    numberOfHorizontalPlotColumns = 1;    
    assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
             >= numberOfFiguresPerPage);
           
    plotHorizMarginCm = 1;
    plotVertMarginCm  = 1.5;           
    pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
    pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;              
    plotHeight  = panelHeight;
    plotWidth   = panelWidth ;
    
    plotConfigGeneric;
    
  case 2
    metricNameList = {'angvel','angvelx','angvely','angvelz','angveln'};

    metricYLim = [ 0, 120; -25,25; -40, 120; -90,60; -30, 120];

    if(flag_motionSequence==1)
      metricYLim(3,:) =[-130,20];
      metricYLim(5,:) =[-30,130];
    end

    metricYAxisLabel = {'Angular Velocity ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)',...
                        'Angular Velocity ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)',...
                        'Angular Velocity ($$^\circ$$/s)'};

    metricYAxisScale = [1,1,1,1,1];

    metricTitle = {'Angular Speed',...
                   'Angular Velocity about X',...
                   'Angular Velocity about Y',...
                   'Angular Velocity about Z',...
                   'Angular Velocity about $$\hat{t}$$'};

    metricPlotHalfPlaneBox  = [0,0,0,0,0];             
    metricPlotHalfPlaneYmax = [0,0,0,0,0];                             
    metricPlotHalfPlaneYmin = metricYLim(:,1);            
    metricAnnotationLines = [0,0,0,0,0];
    metricPlotZeroLine = [0,0,0,0,0];
    plotConfigName = 'ComAngularVelocity';  
        
    flag_plotSubjectData = 1;
    flag_plotGroupData   = 1;
    
    % Plot Config    
    numberOfFiguresPerPage        = length(metricNameList);
    numberOfVerticalPlotRows      = numberOfFiguresPerPage;
    numberOfHorizontalPlotColumns = 1;    
    assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
             >= numberOfFiguresPerPage);
           
    plotHorizMarginCm = 1;
    plotVertMarginCm  = 1.5;           
    pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
    pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;          
    plotHeight  = panelHeight;
    plotWidth   = panelWidth ;
    
    plotConfigGeneric;
    
  case 3
    metricNameList = {'timemaxfpe','fzerrseatoff'};
    metricYLim = [-0.05,0.25; -5, 5];
    
    
    if(flag_motionSequence==1)
      metricNameList = {'timemaxfpe','fzerrseatoff'};
      metricYLim = [-1,0; 0, 85];      
    end
    
    metricYTickSpacing = [0.05,1];    
    metricYAxisLabel = {'Time (s)','Force (N)'};
    metricYAxisScale = [1, 1];
    metricTitle = {'Time since max((FPE-COM).u)','Chair Fz Err at Seat-Off'};
    metricPlotHalfPlaneBox = [0,0];
    metricPlotHalfPlaneYmax = [0,0];
    metricPlotHalfPlaneYmin = metricYLim(:,1);    
    metricPlotZeroLine = [0,0];
    
    metricAnnotationLines = [0,1];
    plotConfigName = 'SwitchingConditions';
    
    flag_plotSubjectData = 1;
    flag_plotGroupData   = 1;
        
    % Plot Config    
    numberOfFiguresPerPage        = length(metricNameList);
    numberOfVerticalPlotRows      = numberOfFiguresPerPage;
    numberOfHorizontalPlotColumns = 1;    
    assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
             >= numberOfFiguresPerPage);
           
    plotHorizMarginCm = 1;
    plotVertMarginCm  = 1.5;           
    pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
    pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;  
                  
    plotHeight  = panelHeight;
    plotWidth   = panelWidth ;
    
    plotConfigGeneric;
    
  case 4
    
    tmp = panelWidth;
    panelWidth = panelHeight;
    panelHeight = tmp;
    
    metricNameList = {'com2edge','com2cop','comvel',...
                      'angvel','angvelz','fpe2edge'};

    metricYLim = [-10.1,7.1;...
                   0,15.1; ...
                   0,50.1;...
                  0,85; ...
                  -20, 25; ...
                  -1,12.1];

    if(flag_motionSequence==1)
      metricYLim = [-10.1,7.1;...
                     0,15.1; ...
                     0,50.1;...
                    0,110; ...
                    -30, 10; ...
                    -18.1,10.1];
    end
                
    metricYAxisLabel = {'Displacement (cm)','Distance (cm)','Speed (cm/s)',...
                        'Angular Speed ($$^\circ$$/s)',...
                        'Angular Velocity ($$^\circ$$/s)','Displacement (cm)'};

    metricYAxisScale = [1,1,1,...
                        1,1,1];

    metricYTick = [5;5;10;20;10;5];                      
                      
%     metricTitle = {'A. COM$_{GP}$ inside BOS: $_{B(C)}d{_C}$',...
%                    'B. COM$_{GP}$-COP Dist.: $|{_C}r{_P}|_2$',...
%                    'C. COM Speed: $|v{_C}|_2$',...
%                    'D. Ang. Speed: $|\omega|_2$',...
%                    'E. FPE Asmpt.: $\omega_Z \approx 0$',...
%                    'F. FPE Dyn. Stab. Margin: $_{B(F)}d{_F}$'};

    metricTitle = {'COM$_{GP}$ to nearest BOS-edge',...
                   'COM$_{GP}$-COP',...
                   'COM Speed',...
                   'Wholebody Angular Speed',...
                   'Assumption: $\omega_Z \approx 0$',...
                   'Dynamic Stability Margin'};


    metricUpperLeftNotes = {'(+) Inside BOS','','',...
                            '','(+) CCW','(+) Inside BOS'};
    metricLowerLeftNotes = {'(-) Outside BOS','','',...
                            '','(-) CW','(-) Outside BOS'};

                            
    metricPlotHalfPlaneText = ...
      {'',...
         '';...
       'Quiet standing:', ...
          sprintf('$|{_C}r{_P}|_2 < %1.1f\\,\\mathrm{cm}$',mean(quietStandingData.comCopDist(3,:)));...
       'Quiet standing:', ...
          sprintf('$|v{_C}|_2 < %1.1f\\,\\mathrm{cm}/\\mathrm{s}$',...
                 mean(quietStandingData.comSpeed(3,:)));...
       'Quiet standing:',sprintf('$|\\omega|_2 < %1.1f\\,^\\circ/s$',mean(quietStandingData.angSpeed(3,:)));...                 
       'Quiet standing:',sprintf('$|\\omega_Z| < %1.1f\\,^\\circ/s$',angZSpeedOneSided);...
       'Dyn. balanced','standing $_{B(F)}d{_F} > 0$'};

%  metricTitle = {'A. Displacement from nearest BOS edge to COM$_{\mathrm{GP}}$ ($_{B(C)}r{_C}$)',...
%                    'B. Distance from COM$_{\mathrm{GP}}$ to COP ($|{_C}r{_P}|_2$)',...
%                    'C. COM Speed ($|v{_C}|_2$)',...
%                    'B. Displacement from BOS edge to FPE ($_{B(F)}r{_F}$)',...
%                    'C. Displacement from COP to FPE in $\hat{t}$ (${_P}r{_F} \cdot \hat{t}$)',...
%                    'D. Displacement from COP to FPE in $\hat{s}$ (${_P}r{_F} \cdot \hat{s}$)',...
%                    'D. Angular speed ($|\omega|_2$)',...
%                    'A. Angular velocity about the vertical axis ($\omega_Z$)'};                 




    metricPlotHalfPlaneBox = [1,1,1,1,1,1];   
    
    metricPlotZeroLine = [0,0,0,0,0,0];
    
      
    
                   
    metricPlotHalfPlaneYmax = [metricYLim(1,2),...
                               mean(quietStandingData.comCopDist(3,:)),...
                               mean(quietStandingData.comSpeed(3,:)),...
                               mean(quietStandingData.angSpeed(3,:)),...
                               angZSpeedOneSided,...
                               metricYLim(6,2)];
                             
    metricPlotHalfPlaneYmin = [0,...
                               0,...
                               0,...
                               0,...
                               -angZSpeedOneSided,...
                               0];


                            
    metricAnnotationLines = [0,0,0,0,0,0];
    
    plotConfigName = 'ComFpeBalanceMetricsSummary';   
       
    flag_plotSubjectData = 0;
    flag_plotGroupData   = 1;
    
    
    % Plot Config
    numberOfFiguresPerPage        = length(metricNameList);
    numberOfVerticalPlotRows      = 3;
    numberOfHorizontalPlotColumns = 2;    
    assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
             >= numberOfFiguresPerPage);

    plotHorizMarginCm = 1.5;
    plotVertMarginCm  = 1.5;    
    
    plotHeight  =  panelHeight / numberOfVerticalPlotRows;
    plotWidth   = (panelWidth)  / numberOfHorizontalPlotColumns;
    
    pageHeight  = numberOfVerticalPlotRows*(plotHeight) + 2;
    pageWidth   = numberOfHorizontalPlotColumns*(plotWidth) + 2;  
    
    plotConfigGeneric;
    
  otherwise
    assert(0);
    
end
    




metricYTickSpacingStruct(size(metricYLim,1))=struct('default',[],'custom',[]);
for z=1:1:size(metricYLim,1)
  dy = metricYLim(z,2)-metricYLim(z,1);
  deltaY = floor(dy/5); 

  deltaYCases = [1,5,10,20];
  deltaYErr = Inf;
  kBest=0;
  deltaY = 0;
  if(isempty(metricYTick)==1)
    for k=1:1:length(deltaYCases)
      if(abs(deltaY-deltaYCases(1,k)) < deltaYErr)
        deltaYErr = abs(deltaY-deltaYCases(1,k));
        kBest=k;
      end
    end
    deltaY = deltaYCases(1,kBest);
  else
    deltaY = metricYTick(z,1);
  end
 
  tickYN = [0];
  if(metricYLim(z,1) < 0)
    tickYN = [-floor( -metricYLim(z,1)/deltaY)*deltaY: deltaY: 0];
  end
  tickYP = [deltaY:deltaY:floor(metricYLim(z,2)/deltaY)*deltaY];

  metricYTickSpacingStruct(z).default = [tickYN,tickYP];
end


%%
%
% Subject/Grouping/Data
%
%%

gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;

groups(2) = struct('name',[],'color',[],'label',[],'x',0);

indexGroupYoung   = 1;
indexGroupElderly = 2;

groups(indexGroupYoung  ).name  = 'Y';
groups(indexGroupYoung  ).label = 'Younger Adult (Y)';
groups(indexGroupElderly).name  = 'O';
groups(indexGroupElderly).label = 'Older Adult (O)';

colorMod = 0.66;

colorA = [0,0,0];
colorB = colorA.*(1-colorMod) + [1,1,1].*colorMod;
groups(indexGroupYoung  ).color   = [colorA;colorB];

colorA = gorgeousGreen;
colorB = colorA.*(1-colorMod) + [1,1,1].*colorMod;
groups(indexGroupElderly).color   = [colorA;...
                                     colorB];





%%
%Generate the plots
%%

figPlotMatrix(length(trialsToProcess) )= struct('h',[]);
  
xDeltaSubjectGroup = 0.5; 
marginWidthForSignificanceStats = 2;

xPlotLim = 0;
xPlotLimSub = 0;
xPlotLimGroup = 0;

if(flag_plotSubjectData == 1 && flag_plotGroupData == 1)
  xPlotLim = size(subjectData,1) ...
           + size(groupData,1) ...
           + xDeltaSubjectGroup ...
           + marginWidthForSignificanceStats;
  xPlotLimSub = size(subjectData,1);
  xPlotLimGroup = size(subjectData,1) ...
           + size(groupData,1) ...
           + xDeltaSubjectGroup;
         
end

if(flag_plotSubjectData == 0 && flag_plotGroupData == 1)
  marginWidthForSignificanceStats = 2.0;
  xPlotLim = size(groupData,1) ...
           + xDeltaSubjectGroup ...
           + marginWidthForSignificanceStats;
  xPlotLimSub = NaN;
  xPlotLimGroup = size(groupData,1) ...
           + xDeltaSubjectGroup;         
end

if(flag_motionSequence==1)
  xPlotLimGroup = xPlotLimGroup+xDeltaSubjectGroup;
end




for indexTrialsToProcess = 1:1:length(trialsToProcess)
  indexTrial = 0;
  for k=1:1:size(subjectData,2)
    if(contains(trialsToProcess{indexTrialsToProcess},subjectData(1,k,1).trialId))
      indexTrial = k;
      assert( contains(trialsToProcess{indexTrialsToProcess},groupData(1,k,1).trialId))
    end
  end  
  assert(indexTrial ~= 0);
  
  

  figPlotMatrix(indexTrialsToProcess).h = figure;
  figH = figPlotMatrix(indexTrialsToProcess).h;
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  
  for indexMetricCount =1:1:length(metricNameList)
    indexMetric = indexMetricCount;
    if(isempty(metricOrder)==0)
      indexMetric = metricOrder(1,indexMetricCount);
    end
    
    metricName = metricNameList{indexMetric};

    [row,col] = find(subPlotPanelIndex==indexMetricCount);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    
    figure(figH);
    if(length(subPlotVec) == 3)
      subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
    end  
    if(length(subPlotVec) == 4)
      subplot('Position',subPlotVec);
    end    
        
    if(metricPlotHalfPlaneBox(1,indexMetric)==1)
      %Plot the grey background box
      x0 = -0.5;
      x1 = xPlotLimGroup;
      y0 = metricYLim(indexMetric,1);
      y1 = metricYLim(indexMetric,2);
      dataBox = [x0,y0;x1,y0;x1,y1;x0,y1;x0,y0];
      
      fill(dataBox(:,1),dataBox(:,2),[1,1,1],'EdgeColor','none');
      hold on;
      
      %Plot the stable white box
      x0 = -0.5;
      x1 = xPlotLimGroup;
      y0 = metricPlotHalfPlaneYmin(1,indexMetric);
      y1 = metricPlotHalfPlaneYmax(1,indexMetric);
      
      val = metricPlotHalfPlaneYmin(1,indexMetric);
      flag_boxUp = 1;
      %Valid if Ymax is the same as the plot limit.
      %Invalid if Ymin is the same as the plot limt
      valOptionError = abs(metricPlotHalfPlaneYmin(1,indexMetric)...
                   -metricYLim(indexMetric,1));
      if(  valOptionError < 1e-3)
        val = metricPlotHalfPlaneYmax(1,indexMetric);      
        flag_boxUp=0;
      end
      
      
      dataBox = [x0,y0;x1,y0;x1,y1;x0,y1;x0,y0];
      fill(dataBox(:,1),dataBox(:,2),[1,1,1].*0.9,'EdgeColor','none');
      hold on; 
      
 
      
      
      dy = 0.1*(metricYLim(indexMetric,2)-metricYLim(indexMetric,1));
      
      dy6pt = (0.3/(panelHeight))*(metricYLim(indexMetric,2)-metricYLim(indexMetric,1));
      
      %text(xPlotLimGroup,y0+dy6pt,'Stable:',...
      %     'HorizontalAlignment','left',...
      %     'VerticalAlignment','bottom','FontSize',6); 
      %hold on;
      
      ly1 = y1;
      ly0 = ly1;
      
      if( ly0 > metricPlotHalfPlaneYmax(1,indexMetric))
        ly0 = metricPlotHalfPlaneYmax(1,indexMetric);
      end
      
      plot([x1;(x1+0.25)],[ly0;ly1],'-','LineWidth',1.,'Color',[1,1,1].*0.9);
      hold on;
      text(xPlotLimGroup+0.25,y1+dy*0.5,...
            metricPlotHalfPlaneText{indexMetric,1},...
           'HorizontalAlignment','left',...
           'VerticalAlignment','top','FontSize',6,...
           'interpreter','latex');
      hold on;
      text(xPlotLimGroup+0.25,y1,...
            metricPlotHalfPlaneText{indexMetric,2},...
           'HorizontalAlignment','left',...
           'VerticalAlignment','top','FontSize',6,...
           'interpreter','latex');
      hold on;
      
      %Get the group summary statistics to generate the custom y ticks
      [phaseResults, startResults, endResults] = ...
        calcGroupComparison(...
        groupData(indexGroupYoung,indexTrial,indexPhase).(metricName),...
        groupData(indexGroupElderly,indexTrial,indexPhase).(metricName));
   
      medY = groupData(indexGroupYoung,indexTrial,indexPhase).(metricName).start.median;      
      medE = groupData(indexGroupElderly,indexTrial,indexPhase).(metricName).start.median;
      
      minY = groupData(indexGroupYoung,indexTrial,indexPhase).(metricName).start.min;      
      minE = groupData(indexGroupElderly,indexTrial,indexPhase).(metricName).start.min;
      
      maxY = groupData(indexGroupYoung,indexTrial,indexPhase).(metricName).start.max;      
      maxE = groupData(indexGroupElderly,indexTrial,indexPhase).(metricName).start.max;
      
      minYE = min([minY,minE]);
      maxYE = max([maxY,maxE]);
      medYE = 0.5*(medY+medE);

      %metricYTickSpacing(indexMetric).custom = sort([minY,medY,maxY,minE,medE,maxE]);
      metricYTickSpacingStruct(indexMetric).custom = sort([minYE,medYE,maxYE]);
      ySpan = metricYLim(indexMetric,2)-metricYLim(indexMetric,1);

%       for z=2:1:(length(metricYTickSpacing(indexMetric).custom)-1)
%         dL = metricYTickSpacing(indexMetric).custom(1,z) ...
%            - metricYTickSpacing(indexMetric).custom(1,z-1);
%         dR = metricYTickSpacing(indexMetric).custom(1,z+1) ...
%            - metricYTickSpacing(indexMetric).custom(1,z);
% 
%         if(dL/ySpan < 0.05 && dR/ySpan < 0.05)
%           metricYTickSpacing(indexMetric).custom(1,z-1) = NaN;
%           metricYTickSpacing(indexMetric).custom(1,z+1) = NaN;
%         end
%         if(dL/ySpan < 0.05 && dR/ySpan > 0.05)
%           metricYTickSpacing(indexMetric).custom(1,z-1) = NaN;
%         end
%         if(dL/ySpan > 0.05 && dR/ySpan < 0.05)
%           metricYTickSpacing(indexMetric).custom(1,z+1) = NaN;
%         end
%          
%       end      

      isValueCloseToZero = 0;
      for z=1:1:length(metricYTickSpacingStruct(indexMetric).custom)
        if( abs(metricYTickSpacingStruct(indexMetric).custom(1,z)/ySpan) < 0.05)
          %metricYTickSpacing(indexMetric).custom(1,z) = NaN;
          isValueCloseToZero = 1;
        end                  
      end
      if(isValueCloseToZero == 0)
        metricYTickSpacingStruct(indexMetric).custom = ...
          sort([0,metricYTickSpacingStruct(indexMetric).custom]);
      end
      
      idx = find( isnan(metricYTickSpacingStruct(indexMetric).custom) == 0) ;
      metricYTickSpacingStruct(indexMetric).custom = ...
        metricYTickSpacingStruct(indexMetric).custom(idx);

      %metricYTickSpacing(indexMetric).custom;

      
      if(flag_useCustomYTicks == 1)
        for z=1:1:length(metricYTickSpacingStruct(indexMetric).custom)
          l0 = metricYTickSpacingStruct(indexMetric).custom(1,z);
          lineColor = [1,1,1].*0.85;
          if(l0 >= y0 && l0 <= y1)
            lineColor = [1,1,1].*0.75;
          end

          lineWidth = 0.1;
          lineStyle = '-';
          if( abs(l0-medYE) < 1e-3)
            lineWidth = 0.5;
            lineStyle = '-';
          end


          plot([-0.5,xPlotLimGroup],[l0,l0],lineStyle,'Color',lineColor,...
               'LineWidth',lineWidth);
          hold on;
        end
      end

      if(flag_useCustomYTicks == 1)
        metricYTickSpacingStruct(indexMetric).custom = ...
          round(metricYTickSpacingStruct(indexMetric).custom,2,'significant');      
        yticks([metricYTickSpacingStruct(indexMetric).custom]);
      else
        metricYTickSpacingStruct(indexMetric).default = ...
          round(metricYTickSpacingStruct(indexMetric).default,2,'significant');      
        yticks([metricYTickSpacingStruct(indexMetric).default]);
      end
      
    else
      x0 = -0.5;
      x1 = xPlotLim;
      plot([x0;x1],[0;0],'Color',[1,1,1].*0.9,'LineWidth',0.5);
      hold on;
    end
    
        
    
    for indexGroup=1:1:size(groupData,1)
 
      groupXPosition = indexGroup;
      axisLimits = [0,xPlotLim, ...
                    metricYLim(indexMetric,1),...
                    metricYLim(indexMetric,2)];      
                  
      if(flag_plotSubjectData == 1)
        subjectsInGroup = groupData(indexGroup,indexTrial,indexPhase).subjectIndex;
                
        for indexGroupMember=1:1:length(subjectsInGroup)
          indexSubject = subjectsInGroup(1,indexGroupMember);

          if( isempty(subjectData(indexSubject,indexTrial,indexPhase).subjectId) == 0)



            subjectId = subjectData(indexSubject,indexTrial,indexPhase).subjectId(2:3);
            if(contains(subjectId,'0'))
              idx = strfind(subjectId,'0');
              if(idx == 1)
                subjectId = subjectId(2);
              end
            end                    

            fontColor = groups(indexGroup).color(1,:);

            flagEnableMinMaxLabels=1;
            [figH] = plotMetricDistributionEventData2(figH, subPlotVec, indexSubject,...
                  subjectData(indexSubject,indexTrial,indexPhase).(metricName),...
                  metricYAxisScale(1,indexMetric),...
                  subjectId, groups(indexGroup).color, axisLimits, ...
                  boxWidth, lineWidth, plotFontName, fontColor,...
                  flagPlotStart,flagPlotPhase,flagPlotEnd,...
                  flag_motionSegmentEmphasis,...
                  flagEnableMinMaxLabels);
            hold on;

            
            
          end
        end
        
        groupXPosition =   size(subjectData,1) + indexGroup+xDeltaSubjectGroup;
      end      
             
      if(flag_plotGroupData == 1)
      
        groups(indexGroup).x = groupXPosition;

        fontColor = groups(indexGroup).color(1,:);        
        flagEnableMinMaxLabels = 0;
        
        flagPlotStartTemp = flag_motionSegmentEmphasis == 0;
        flagPlotPhaseTemp = flag_motionSegmentEmphasis == 1;
        flagPlotEndTemp   = flag_motionSegmentEmphasis == 2;
        
        figH = plotMetricDistributionEventData2(figH, subPlotVec, ...
                groupXPosition,...
                groupData(indexGroup,indexTrial,indexPhase).(metricName),...
                metricYAxisScale(1,indexMetric),...
                groups(indexGroup).name, groups(indexGroup).color, ...
                axisLimits,boxWidth,lineWidth,plotFontName,fontColor,...
                flagPlotStartTemp,flagPlotPhaseTemp,flagPlotEndTemp,...
                flag_motionSegmentEmphasis,...
                flagEnableMinMaxLabels);
          hold on;

        %Compute the Wilcoxin rank sum test  
        if(indexGroup == size(groupData,1))
          assert(indexGroup == 2)
          [phaseResults, startResults, endResults] = ...
            calcGroupComparison(...
            groupData(indexGroupYoung,indexTrial,indexPhase).(metricName),...
            groupData(indexGroupElderly,indexTrial,indexPhase).(metricName));

          %disp(metricName);
          %disp(startResults);
          %disp('Young');
          %disp(groupData(indexGroupYoung,indexTrial,indexPhase).(metricName).start)
          %disp('Elderly');          
          %disp(groupData(indexGroupElderly,indexTrial,indexPhase).(metricName).start)
          
          startResultsPlot  = [];
          phaseResultsPlot  = [];
          endResultsPlot    = [];
          
          fontColor = [0,0,0];
          switch(flag_motionSegmentEmphasis)
            case 0
              startResultsPlot = startResults;
            case 1
              phaseResultsPlot = phaseResults;
            case 2
              endResultsPlot = endResults;
            otherwise
              assert(0);
          end
          
          figH = plotGroupComparisons(figH, subPlotVec,...
                   groups(1).x, groupData(1,indexTrial,indexPhase).(metricName),...
                   groups(2).x, groupData(2,indexTrial,indexPhase).(metricName),...
                   phaseResultsPlot, startResultsPlot, endResultsPlot,...
                   phaseLabels{2}, phaseLabels{1}, phaseLabels{3},...
                   metricAnnotationLines(1,indexMetric),boxWidth,...
                   plotFontName,fontColor);          

        end
      end
      
      if(indexGroup == 1)          
        axisLimitsScaled = axis;
        ylabel(metricYAxisLabel{indexMetric});
        titleFontSize =  get(groot,'defaultAxesFontSize')...
                        *get(groot,'defaultAxesTitleFontSizeMultiplier');

        xTitle = axisLimitsScaled(1,1); 
        

        if(flag_plotConfig==4)
          xTitle = xTitle - 0.1*( axisLimitsScaled(1,2)-axisLimitsScaled(1,1) );
        end
        yTitle = axisLimitsScaled(1,4) ...
                + 0.1*(axisLimitsScaled(1,4)-axisLimitsScaled(1,3));          

        %'Interpreter','latex',      
        text(xTitle,yTitle,metricTitle{indexMetric},...
             'FontSize',titleFontSize,...
             'HorizontalAlignment','left',...
             'fontname',plotFontName);
        hold on;

        maxSubject = size(subjectData,1);
        x0 = maxSubject+xDeltaSubjectGroup*1.5;
        y0 = axisLimitsScaled(1,3);
        y1 = axisLimitsScaled(1,4);
        dy =y1-y0;
        plot([x0;x0],[(y0-0.05*dy);(y1+0.05*dy)],...
             '-','Color',[1,1,1],'LineWidth',2);
        hold on;
        
        plot([x0;x0],[(y0-0.05*dy);(y1+0.05*dy)],...
             '-','Color',[1,1,1].*0.5,'LineWidth',0.5);
        hold on;
        
        if(isempty(metricUpperLeftNotes)==0)
          text(0.125, y1, metricUpperLeftNotes{1,indexMetric},...
               'FontSize',6,...
                'VerticalAlignment','top',...
                'HorizontalAlignment','left',...
                'interpreter','latex');        
          hold on;
        end
        if(isempty(metricLowerLeftNotes)==0)        
          text(0.125, y0, metricLowerLeftNotes{1,indexMetric},...
                'FontSize',6,...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','left',...
                'interpreter','latex');
          hold on;
        end
      end      
      

      axisLimits = [0,xPlotLim, ...
                    metricYLim(indexMetric,1),...
                    metricYLim(indexMetric,2)]; 
                  
      axis(axisLimits);
      
      tickLength = get(gca,'TickLength');
      set(gca,'TickLength',tickLength.*0.5);
      set(gca,'XTickLabel',[]); 
      set(gca, 'XTick',[]);
      box off;
      
    end 
  end
end

for indexTrialsToProcess=1:1:length(trialsToProcess)
  formatTag = '';
  if(flag_PaperPlots0Presentation1 == 1)
    formatTag = 'Beamer';
  end
  
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  configPlotExporter;
  print('-dpdf',[outputPlotFolder,'fig2_Results',...
                  formatTag,...
                  trialsToProcess{indexTrialsToProcess},...
                  plotConfigName,... 
                  bosModelTag,nameToeTag,nameExcludeE06Tag,...
                  motionSequenceTag,'.pdf']);  
end  
  
        

