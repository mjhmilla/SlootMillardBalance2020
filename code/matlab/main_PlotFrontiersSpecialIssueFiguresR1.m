
clc;
close all;
clear all;

%%
%
% Constants
%
%%

flag_PaperPlots0Presentation1 = 1;

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

flag_plotConfig = 4; 
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


switch flag_plotConfig
  case 0
    
    if(flag_PaperPlots0Presentation1)
      panelHeight = panelHeight*0.5;
    end
        
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
 
  case 1
    metricNameList = {'comvel','comvelx','comvely','comvelz'};

    metricYLim = [ 0,65; -20,65; -10,10; -40,65];
         

    metricYAxisLabel = {'Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)','Velocity (cm/s)'};

    metricYAxisScale = [1,1,1,1];

    metricTitle = {'COM Speed',...
                   'COM Velocity in X',...
                   'COM Velocity in Y',...
                   'COM Velocity in Z'};

    metricPlotHalfPlaneBox = [0,0,0,0];             

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
    metricYAxisLabel = {'Time (s)','Force (N)'};
    metricYAxisScale = [1, 1];
    metricTitle = {'Time since max((FPE-COM).u)','Chair Fz Err at Seat-Off'};
    metricPlotHalfPlaneBox = [0,0];
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

    metricYLim = [-1,17;    0,15;  0,65;...
                  0,115; -25, 25; -1,17];
    

    metricYAxisLabel = {'Distance (cm)','Distance (cm)','Velocity (cm/s)',...
                        'Angular Speed ($$^\circ$$/s)',...
                        'Angular Speed ($$^\circ$$/s)','Distance (cm)'};

    metricYAxisScale = [1,1,1,...
                        1,1,1];

    metricTitle = {'A. COM$_{GP}$ inside BOS: $d_{COM}$',...
                   'B. COM$_{GP}$-COP Dist.: $_{COM}d_{COP}$',...
                   'C. COM Speed: $|v_{COM}|_2$',...
                   'D. Ang. Speed: $|\omega_{avg}|_2$',...
                   'E. FPE Asmpt.: $\omega_{avg}\cdot\hat{k} \approx 0$',...
                   'F. FPE Dyn. Stab. Margin: $d_{FPE}$'};

    metricPlotHalfPlaneBox = [1,0,0,...
                              0,0,1];             

    metricAnnotationLines = [1,0,0,...
                             0,0,1];
    
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

xPlotLim = 0;

if(flag_plotSubjectData == 1 && flag_plotGroupData == 1)
  xPlotLim = size(subjectData,1) ...
           + size(groupData,1) ...
           + xDeltaSubjectGroup ...
           + marginWidthForSignificanceStats;
end

if(flag_plotSubjectData == 0 && flag_plotGroupData == 1)
  marginWidthForSignificanceStats = 1.5;
  xPlotLim = size(groupData,1) ...
           + xDeltaSubjectGroup ...
           + marginWidthForSignificanceStats;
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

            figH = plotMetricDistributionEventData(figH, subPlotVec, indexSubject,...
                  subjectData(indexSubject,indexTrial,indexPhase).(metricName),...
                  metricYAxisScale(1,indexMetric),...
                  subjectId, groups(indexGroup).color, axisLimits, ...
                  boxWidth, lineWidth, plotFontName, fontColor);
            hold on;

          end
        end
        
        groupXPosition =   size(subjectData,1) + indexGroup;
      end      
             
      if(flag_plotGroupData == 1)
      
        groups(indexGroup).x = groupXPosition;

        fontColor = groups(indexGroup).color(1,:);

        figH = plotMetricDistributionEventData(figH, subPlotVec, ...
                groupXPosition,...
                groupData(indexGroup,indexTrial,indexPhase).(metricName),...
                metricYAxisScale(1,indexMetric),...
                groups(indexGroup).name, groups(indexGroup).color, ...
                axisLimits,boxWidth,lineWidth,plotFontName,fontColor);
          hold on;

        %Compute the Wilcoxin rank sum test  
        if(indexGroup == size(groupData,1))
          assert(indexGroup == 2)
          [phaseResults, startResults, endResults] = ...
            calcGroupComparison(...
            groupData(indexGroupYoung,indexTrial,indexPhase).(metricName),...
            groupData(indexGroupElderly,indexTrial,indexPhase).(metricName));

          fontColor = [0,0,0];
          figH = plotGroupComparisons(figH, subPlotVec,...
                   groups(1).x, groupData(1,indexTrial,indexPhase).(metricName),...
                   groups(2).x, groupData(2,indexTrial,indexPhase).(metricName),...
                   phaseResults, startResults, [],...
                   'all', 'seat-off', 'standing',...
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
          xTitle = xTitle - 0.3*( axisLimitsScaled(1,2)-axisLimitsScaled(1,1) );
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
end

for indexTrialsToProcess=1:1:length(trialsToProcess)
  formatTag = '';
  if(flag_PaperPlots0Presentation1 == 1)
    formatTag = 'Beamer';
  end
  
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  configPlotExporter;
  print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                  formatTag,...
                  trialsToProcess{indexTrialsToProcess},...
                  plotConfigName,... 
                  nameToeTag,'.pdf']);  
end  
  
        

