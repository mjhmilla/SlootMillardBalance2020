
clc;
close all;
clear all;

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

trialsToProcess = {'Chest','Conv','Leg','Side','Rob'};


metricNameList = {'duration','com2edge','com2cop','comvel',...
                  'fpe2edge','fpewidth','fpelen','fpeasmhkz'};
                
metricYLim = [  0,8;  -6,17; 0,17; 0,65;...
              -6,17; -5,10; -5,10.5; 0,0.015];
offsetYLable = -0.05;            
            
metricYAxisLabel = {'Time (s)','Distance (cm)','Distance (cm)','Velocity (cm/s)',...
                    'Distance (cm)','Distance (cm)','Distance (cm)','Percentage (\%)'};
                  
metricYAxisScale = [1,1,1,1,...
                     1,1,1,100];
                                      
metricTitle = {' Duration',...
               'A. COM to nearest BOS edge',...
               'B. COM to COP',...
               'C. COM Speed',...
               'A. FPE to nearest BOS edge',...
               'B. FPE-COP in $\hat{t}$',...
               'C. FPE-COP in $\hat{s}$',...
               'FPE Assumption: $$\epsilon = 1- H_{GP}\cdot\hat{k}/|H_{GP}| \approx 0 $$'};

metricPlotHalfPlaneBox = [0,1,0,0,1,0,0,0];             
             
numberOfFiguresPerPage = length(metricNameList);
plotConfigFrontiers;


%%
%
% Subject/Grouping/Data
%
%%

gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;

groups(2) = struct('name',[],'color',[],'label',[]);

indexGroupYoung   = 1;
indexGroupElderly = 2;

groups(indexGroupYoung  ).name  = 'Y';
groups(indexGroupYoung  ).label = 'Young (Y)';
groups(indexGroupElderly).name  = 'E';
groups(indexGroupElderly).label = 'Elderly (E)';

colorMod = 0.75;

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
xPlotLim = size(subjectData,1) + size(groupData,1) + xDeltaSubjectGroup;





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
      fill(dataBox(:,1),dataBox(:,2),[1,1,1].*0.5,'EdgeColor','none');
      hold on;      
    else
      x0 = -0.5;
      x1 = xPlotLim;
      plot([x0;x1],[0;0],'Color',[1,1,1].*0.5,'LineWidth',0.5);
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
                subjectId, groups(indexGroup).color, axisLimits,plotFontName);
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
               
      if(indexGroup==2 && indexTrial == 5 && indexPhase == 2 && indexMetric==1)
        here=1;
      end
      figH = plotMetricDistributionEventData(figH, subPlotVec, xPosition,...
              groupData(indexGroup,indexTrial,indexPhase).(metricName),...
              metricYAxisScale(1,indexMetric),...
              groups(indexGroup).name, groups(indexGroup).color, ...
              axisLimits,plotFontName);
        hold on;
        
      %Compute the Wilcoxin rank sum test  
      if(indexGroup == size(groupData,1))
        assert(indexGroup == 2)
        [phaseResults, startResults, endResults] = ...
          calcGroupComparison(...
          groupData(indexGroupYoung,indexTrial,indexPhase).(metricName),...
          groupData(indexGroupElderly,indexTrial,indexPhase).(metricName));
        
        x0 = size(subjectData,1)+1;
        x1 = size(subjectData,1)+2;
        axisLim = axis;
        yA = axisLim(1,3);
        yB = axisLim(1,4);
        dy = yB-yA;
        y1 = yB-2*0.05*dy;
        y0 = yB-3*0.05*dy;
        plot([x0;x0],[y0;y1],'-','Color',[0,0,0],'LineWidth',0.5);
        hold on;
        plot([x1;x1],[y0;y1],'-','Color',[0,0,0],'LineWidth',0.5);
        hold on;
        plot([x0;x1],[y1;y1],'-','Color',[0,0,0],'LineWidth',0.5);
        hold on;
        
        starText = '';
        if(phaseResults.h==1)
          starText = '*';
        else
          starText = '';
        end

        
        
        
        statsText = sprintf('%sp = %1.2e',starText,phaseResults.p);
        text( x0,yB+0.05*dy,statsText,...
            'FontSize',6,'HorizontalAlignment','left',...
            'fontname',plotFontName);
          hold on;        
        
      end
    end 
  end
end

for indexTrialsToProcess=1:1:length(trialsToProcess)
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  configPlotExporter;
  print('-dpdf',['../../plots/Frontiers2020Pub/fig_Results',...
                  trialsToProcess{indexTrialsToProcess},'.pdf']);  
end  
  
        

