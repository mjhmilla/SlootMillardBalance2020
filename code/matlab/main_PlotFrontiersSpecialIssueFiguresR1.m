
clc;
close all;
clear all;

%%
%
% Constants
%
%%

flag_PaperPlots0Presentation1 = 0;

flag_ConvexHullWithToes0Without1  = 0;
nameToeTag = '';
if(flag_ConvexHullWithToes0Without1 == 1)
  nameToeTag = '_NoToes';
end


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



indexPhase = indexPhaseSeatOff2Stand;

trialsToProcess = {'Side','Chest','Conv','Leg','Side','Rob'};

flag_plotConfig = 0; 
% 0. Frontiers balance metrics
% 1. Com velocity
% 2. Com angular velocity
% 3. Time since peak (fpe-com).u
plotConfigName = '';

switch flag_plotConfig
  case 0
    metricNameList = {'com2edge','com2cop','comvel',...
                      'fpe2edge','fpewidth','fpelen','angvel','angvelz'};

    metricYLim = [-6,17; 0,17; 0,65;...
                  -6,17; -5,10; -5,10.5; 0,100; -90, 60];
    

    metricYAxisLabel = {'Distance (cm)','Distance (cm)','Velocity (cm/s)',...
                        'Distance (cm)','Distance (cm)','Distance (cm)',...
                        'Angular Velocity ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)'};

    metricYAxisScale = [1,1,1,1,...
                         1,1,1,1];

    metricTitle = {'A. COM$_{\mathrm{GP}}$ to nearest BOS edge',...
                   'B. COM$_{\mathrm{GP}}$ to COP',...
                   'C. COM Speed',...
                   'A. FPE to nearest BOS edge',...
                   'B. FPE-COP in $\hat{t}$',...
                   'C. FPE-COP in $\hat{s}$',...
                   'A. Angular speed',...
                   'B. Angular velocity about the vertical axis'};

    metricPlotHalfPlaneBox = [1,0,0,1,0,0,0,0];             

    metricAnnotationLines = [1,0,0,1,0,0,1,0];
    
    plotConfigName = 'ComFpeBalanceMetrics';
  case 1
    metricNameList = {'comvel','comvelx','comvely','comvelz'};

    metricYLim = [ 0,65; -5,65; -5,5; -5,65];
         

    metricYAxisLabel = {'Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)'};

    metricYAxisScale = [1,1,1,1];

    metricTitle = {'COM Speed',...
                   'COM Velocity in X',...
                   'COM Velocity in Y',...
                   'COM Velocity in Z'};

    metricPlotHalfPlaneBox = [0,0,0,0];             

    metricAnnotationLines = [0,0,0,0];

    plotConfigName = 'ComVelocity';    
    
  case 2
    metricNameList = {'angvel','angvelx','angvely','angvelz','angveln'};

    metricYLim = [ 0, 100; -25,25; -40, 100; -90,60; -30, 110];
         

    metricYAxisLabel = {'Angular Velocity ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)',...
                        'Angular Velocity ($$^\circ$$/s)','Angular Velocity ($$^\circ$$/s)',...
                        'Angular Velocity ($$^\circ$$/s)'};

    metricYAxisScale = [1,1,1,1,1];

    metricTitle = {'Angular Speed',...
                   'Angular Velocity about X',...
                   'Angular Velocity about Y',...
                   'Angular Velocity about Z',...
                   'Angular Velocity about $$\hat{t}$$'};

    metricPlotHalfPlaneBox = [0,0,0,0,0];             

    metricAnnotationLines = [0,0,0,0,0];

    plotConfigName = 'ComAngularVelocity';  
    
  case 3
      metricNameList = {'timemaxfpe'};
      metricYLim = [-0.05,0.150];
      metricYAxisLabel = {'Time'};
      metricYAxisScale = [1];
      metricTitle = {'Time since max((FPE-COM).u)'};
      metricPlotHalfPlaneBox = [0];
      metricAnnotationLines = [0];
      plotConfigName = 'MaxFpeTimeOffset';
  otherwise assert(0);
    
end
    
numberOfFiguresPerPage = length(metricNameList);

boxWidth = 0.33;
      
switch flag_PaperPlots0Presentation1
  case 0
    plotConfigFrontiers;
  case 1
    plotConfigPresentation;
  otherwise assert(0);
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
groups(indexGroupYoung  ).label = 'Young (Y)';
groups(indexGroupElderly).name  = 'E';
groups(indexGroupElderly).label = 'Elderly (E)';

colorMod = 0.66;

colorA = [0,0,0];
colorB = colorA.*(1-colorMod) + [1,1,1].*colorMod;
groups(indexGroupYoung  ).color   = [colorA;colorB];

colorA = gorgeousGreen;
colorB = colorA.*(1-colorMod) + [1,1,1].*colorMod;
groups(indexGroupElderly).color   = [colorA;...
                                     colorB];




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

%%
%
% Load the data
%
%%

load([frontiersDataDir,'groupTrialPhaseMetricData',nameToeTag,'.mat']);
load([frontiersDataDir,'subjectTrialPhaseMetricData',nameToeTag,'.mat']);

  

figPlotMatrix(length(trialsToProcess) )= struct('h',[]);
  
xDeltaSubjectGroup = 0.5; 
marginWidthForSignificanceStats = 2;

xPlotLim = size(subjectData,1) ...
         + size(groupData,1) ...
         + xDeltaSubjectGroup ...
         + marginWidthForSignificanceStats;





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
  
  for indexMetric =1:1:length(metricNameList)
    metricName = metricNameList{indexMetric};

    [row,col] = find(subPlotPanelIndex==indexMetric);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    
    figure(figH);
    if(length(subPlotVec) == 3)
      subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
    end  
    if(length(subPlotVec) == 4)
      subplot('Position',subPlotVec);
    end    
    
    if(metricPlotHalfPlaneBox(1,indexMetric)==1)
      x0 = -0.5;
      x1 = xPlotLim;
      y0 = -10;
      y1 = 0;
      dataBox = [x0,y0;x1,y0;x1,y1;x0,y1;x0,y0];
      fill(dataBox(:,1),dataBox(:,2),[1,1,1].*0.9,'EdgeColor','none');
      hold on;      
    else
      x0 = -0.5;
      x1 = xPlotLim;
      plot([x0;x1],[0;0],'Color',[1,1,1].*0.75,'LineWidth',0.5);
      hold on;
    end
    
        
    
    for indexGroup=1:1:size(groupData,1)
 
      subjectsInGroup = groupData(indexGroup,indexTrial,indexPhase).subjectIndex;
      for indexGroupMember=1:1:length(subjectsInGroup)
        indexSubject = subjectsInGroup(1,indexGroupMember);
        
        if( isempty(subjectData(indexSubject,indexTrial,indexPhase).subjectId) == 0)
                
          axisLimits = [0,xPlotLim, ...
                        metricYLim(indexMetric,1),...
                        metricYLim(indexMetric,2)];

          subjectId = subjectData(indexSubject,indexTrial,indexPhase).subjectId(2:3);
          if(contains(subjectId,'0'))
            idx = strfind(subjectId,'0');
            if(idx == 1)
              subjectId = subjectId(2);
            end
          end                    


          figH = plotMetricDistributionEventData(figH, subPlotVec, indexSubject,...
                subjectData(indexSubject,indexTrial,indexPhase).(metricName),...
                metricYAxisScale(1,indexMetric),...
                subjectId, groups(indexGroup).color, axisLimits, boxWidth, plotFontName);
          hold on;

          if(indexSubject == 1)          
            axisLimitsScaled = axis;
            ylabel(metricYAxisLabel{indexMetric});
            titleFontSize =  get(groot,'defaultAxesFontSize')...
                            *get(groot,'defaultAxesTitleFontSizeMultiplier');

            xTitle = axisLimitsScaled(1,1);
            yTitle = axisLimitsScaled(1,4) ...
                    + 0.1*(axisLimitsScaled(1,4)-axisLimitsScaled(1,3));          

            %'Interpreter','latex',      
            text(xTitle,yTitle,metricTitle{indexMetric},...
                 'FontSize',titleFontSize,...
                 'HorizontalAlignment','left',...
                 'fontname',plotFontName);
            hold on;

            maxSubject = size(subjectData,1);
            x0 = maxSubject+xDeltaSubjectGroup;
            y0 = axisLimitsScaled(1,3);
            y1 = axisLimitsScaled(1,4);
            dy =y1-y0;
            plot([x0;x0],[(y0+0.05*dy);(y1-0.05*dy)],...
                 '-','Color',[1,1,1].*0.5,'LineWidth',0.5)
            hold on;
          end


          set(gca,'TickLength',[0 0]);
          set(gca,'XTickLabel',[]); 
          box off;

        end
      end
      
      xPosition =   size(subjectData,1)...
                   + indexGroup;
               
      groups(indexGroup).x = xPosition;
      
      figH = plotMetricDistributionEventData(figH, subPlotVec, xPosition,...
              groupData(indexGroup,indexTrial,indexPhase).(metricName),...
              metricYAxisScale(1,indexMetric),...
              groups(indexGroup).name, groups(indexGroup).color, ...
              axisLimits,boxWidth,plotFontName);
        hold on;
        
      %Compute the Wilcoxin rank sum test  
      if(indexGroup == size(groupData,1))
        assert(indexGroup == 2)
        [phaseResults, startResults, endResults] = ...
          calcGroupComparison(...
          groupData(indexGroupYoung,indexTrial,indexPhase).(metricName),...
          groupData(indexGroupElderly,indexTrial,indexPhase).(metricName));
        
        figH = plotGroupComparisons(figH, subPlotVec,...
                 groups(1).x, groupData(1,indexTrial,indexPhase).(metricName),...
                 groups(2).x, groupData(2,indexTrial,indexPhase).(metricName),...
                 phaseResults, startResults, [],...
                 'all', 'seat-off', 'standing',...
                 metricAnnotationLines(1,indexMetric),boxWidth,plotFontName);
        
        
     
        
      end
    end 
  end
end

for indexTrialsToProcess=1:1:length(trialsToProcess)
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  configPlotExporter;
  print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                  trialsToProcess{indexTrialsToProcess},...
                  plotConfigName,... 
                  nameToeTag,'.pdf']);  
end  
  
        

