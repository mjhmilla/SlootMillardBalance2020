
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

trialsToProcess = {'Side'};%^{'Chest','Conv','Leg','Side','Rob'};


metricNameList = {'duration','com2edge','com2cop','comvel',...
                  'fpe2edge','fpewidth','fpelen','fpeasmhkz'};
                
metricYLim = [  0,7;  -6,17; 0,13; 0,60;...
              -0.5,16; -3,4.5; -3,9; 0,0.1];
offsetYLable = -0.05;            
            
metricYAxisLabel = {'Time (s)','Distance (cm)','Distance (cm)','Velocity (cm/s)',...
                    'Distance (cm)','Distance (cm)','Distance (cm)','Percentage (%)'};
                  
metricYAxisScxale = [1,1,1,1;...
                     1,1,1,100];
                   
metricTitle = {' Duration',...
               'A. COM to nearest BOS edge',...
               'B. COP to nearest BOS edge',...
               'C. COM Speed',...
               'A. FPE to nearest BOS edge',...
               'B. FPE-COP in $\hat{t}$',...
               'C. FPE-COP in $\hat{s}$',...
               '$\epsilon = 1 - (H_{GP}\dot\hat{k})/|H_{GP}|$'};
                   
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

groups(indexGroupYoung  ).color   = [0,0,0];
groups(indexGroupElderly).color   = gorgeousGreen;



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
    end
  end  

  if(indexTrialsToProcess == 1)
    figPlotMatrix(indexTrialsToProcess).h = figure;
    flag_firstCall=1;
  end
  figH = figPlotMatrix(indexTrialsToProcess).h;
  figure(figPlotMatrix(indexTrialsToProcess).h);  
  
  for indexMetric =1:1:length(metricNameList)
    metricName = metricNameList{indexMetric};

    [row,col] = find(subPlotPanelIndex==indexMetric);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);    
    
    for indexGroup=1:1:size(groupData,1)


        
      subjectsInGroup = groupData(indexGroup,indexTrial,indexPhase).subjectIndex;
      for indexGroupMember=1:1:length(subjectsInGroup)
        indexSubject = subjectsInGroup(1,indexGroupMember);
        
        if(isempty(subjectData(indexSubject,indexTrial,indexPhase).(metricName) )==0)

          phaseWidth = 0.33;
          startEndWidth = phaseWidth/2;
          typeBoxAndWhisker=0;
          typeDotAndWhisker=1;

          if( isfield(subjectData(indexSubject,indexTrial,indexPhase).(metricName),'phase'))        
            figH = plotDistributionData(...
              figH,subPlotVec, indexSubject,...
              subjectData(indexSubject,indexTrial,indexPhase).(metricName).phase,...
              groups(indexGroup).color,...
              phaseWidth,...
              groups(indexGroup).color,...
              typeBoxAndWhisker);
          end

          if( isfield(subjectData(indexSubject,indexTrial,indexPhase).(metricName),'start'))        
            startPosition = indexSubject-phaseWidth;
            startMedian = subjectData(indexSubject,indexTrial,indexPhase... 
                                       ).(metricName).start.median;
            plot([startPosition,indexSubject],[startMedian,startMedian],...
                 'Color',groups(indexGroup).color,'LineWidth',0.5);
            hold on;

            figH = plotDistributionData(...
              figH,subPlotVec, startPosition,...
              subjectData(indexSubject,indexTrial,indexPhase).(metricName).start,...
              groups(indexGroup).color,...
              startEndWidth,...
              [1,1,1],...
              typeDotAndWhisker);
          end
          if( isfield(subjectData(indexSubject,indexTrial,indexPhase).(metricName),'end'))
            endPosition = indexSubject+phaseWidth;
            endMedian = subjectData(indexSubject,indexTrial,indexPhase... 
                                       ).(metricName).end.median;
            plot([endPosition,indexSubject],[endMedian,endMedian],...
                 'Color',groups(indexGroup).color,'LineWidth',0.5);
            hold on;

            figH = plotDistributionData(...
              figH,subPlotVec, endPosition,...
              subjectData(indexSubject,indexTrial,indexPhase).(metricName).end,...
              groups(indexGroup).color,...
              startEndWidth,...
              groups(indexGroup).color,...
              typeDotAndWhisker);        
          end
          axis([0, xPlotLim, metricYLim(indexMetric,1),metricYLim(indexMetric,2)]);
          y0 = metricYLim(indexMetric,1);
          dy = metricYLim(indexMetric,2)-metricYLim(indexMetric,1);
          textYOffset = y0-0.05*dy;

          subjectId = subjectData(indexSubject,indexTrial,indexPhase).subjectId(2:3);
          if(contains(subjectId,'0'))
            idx = strfind(subjectId,'0');
            if(idx == 1)
              subjectId = subjectId(2);
            end
          end

          text(indexSubject,textYOffset,subjectId,...
            'FontSize',8,'Interpreter','latex','HorizontalAlignment','center');
          hold on;
          set(gca,'TickLength',[0 0]);
          set(gca,'XTickLabel',[]); 
          box off;
        end
      end
    end 
  end
end


        

