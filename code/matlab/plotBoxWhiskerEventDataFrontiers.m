function [figH, dataSummary, dataRawVec] = plotBoxWhiskerEventDataFrontiers(...
                                 figH, subPlotVec, ...
                                 movementSequence, ...
                                 xPoint, xPointScale,...                   
                                 dataSeries, dataScale, ...
                                 dataIndex,...
                                 lineColor,  ...
                                 boxWidth,...
                                 dataLabel, dataLabelYLocation,...
                                 xLabelText, ...
                                 yLabelText, ...
                                 titleText,...
                                 axisLimits,...
                                 flag_IntervalEventMode, flag_drawBox,...
                                 flag_firstCall)

                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end

if(flag_drawBox==1 && flag_firstCall)
   boxPts = [-10,-10; 2000,-10;2000,0;-10,0;-10,-10];
  fill(boxPts(:,1),boxPts(:,2),[1,1,1].*0.9,'EdgeColor','none');
  hold on;

end
  
data = [];
dataSummary =struct('min',NaN,'p25',NaN,'mean',NaN,'p75',NaN','max',NaN,...
                    'events',zeros(1,2).*NaN);
                  
dataRaw = struct('x',[],'y',[]);
dataRawVec = [];

eventData = [];
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
 
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    idxA = 0;
    idxB = 0;
    switch(flag_IntervalEventMode)
      case 0
        idxA = idx0;
        idxB = idx2;
      case 1
        idxA = idx1;
        idxB = idx2;
      otherwise assert(0);
      end

    yData = dataSeries(idxA:idxB,:).*dataScale;
    
    eventData = [eventData;...
                  dataSeries(idxA,:).*dataScale,...
                  dataSeries(idxB,:).*dataScale];
                
    data      = [data; yData];
    
    dataRaw.x = dataIndex(idxA:idxB,1);
    dataRaw.y = dataSeries(idxA:idxB,1).*dataScale;
    
    dataRawVec = [dataRawVec;dataRaw];
  end

end

if(isempty(data) == 0)

  
  dataSorted = sort(data(:,1));
  n = length(dataSorted);
  
  n25 = round(n*0.25);
  n75 = round(n*0.75);
  
  dataSummary.min=dataSorted(1,1);
  dataSummary.mean=mean(data(:,1));
  dataSummary.max=dataSorted(end,1);

  dataSummary.p25 = dataSorted(n25,1);
  dataSummary.p75 = dataSorted(n75,1);
  
  for k=1:1:size(dataSummary.events,2)
    dataSummary.events(1,k) = mean(eventData(:,k));
  end
  xDataMid = (xPoint*xPointScale);
  
  figH = plotAntsOnALog(figH,subPlotVec, ...
            xDataMid,...
            dataSummary.min, dataSummary.mean, dataSummary.max,...
            dataSummary.p25, dataSummary.p75, eventData,...
            {'o','o'},[1,1,1;lineColor],boxWidth,lineColor);


        


%   minData  = min(data(:,1));
%   maxData  = max(data(:,1));
%   meanData = mean(data(:,1));
%   stdData  = std(data(:,1));
%   
%   dataSorted = sort(data(:,1));
%   n = length(dataSorted);
%   n25 = round(n*0.25);
%   n75 = round(n*0.75);
%   
%   
% 
%   xDataMid = xPoint*xPointScale;
% 
%   boxStd = [ xDataMid+boxStdWidth*0.5, dataSorted(n25,1);...
%              xDataMid+boxStdWidth*0.5, dataSorted(n75,1);... 
%              xDataMid-boxStdWidth*0.5, dataSorted(n75,1);... 
%              xDataMid-boxStdWidth*0.5, dataSorted(n25,1);... 
%              xDataMid+boxStdWidth*0.5, dataSorted(n25,1)];
% 
%   plot([xDataMid;xDataMid],[minData;maxData],'-','Color',lineColor,...
%        'LineWidth',0.5);
%   hold on;
% 
%   plot([xDataMid],[minData],'.','Color',lineColor,'MarkerSize',3);
%   hold on;
%   
%   plot([xDataMid],[maxData],'.','Color',lineColor,'MarkerSize',3);
%   hold on;
%   
% 
% 
%   for j=1:1:size(eventData,1)
%     plot([xDataMid-boxStdWidth, xDataMid],...
%          [eventData(j,1),eventData(j,1)],'-','Color',lineColor,...
%          'LineWidth',0.5);
%     hold on;    
%     plot([xDataMid;xDataMid+boxStdWidth],...
%          [eventData(j,1),eventData(j,1)],'-','Color',lineColor,...
%          'LineWidth',0.5);
%     hold on;    
% 
% 
%     plot([xDataMid-boxStdWidth*1.1],...
%          [eventData(j,1)],'o','Color',lineColor,...
%          'MarkerSize',2,'MarkerFaceColor',[1,1,1],...
%          'LineWidth',0.5);
%     hold on;    
%     plot([xDataMid+boxStdWidth*1.1],...
%          [eventData(j,2)],'o','Color',lineColor,...
%          'MarkerSize',2,'MarkerFaceColor',lineColor,...
%          'LineWidth',0.5);
%     hold on;    
% 
%   end
%   
%   
%   fill(boxStd(:,1),boxStd(:,2),lineColor,'EdgeColor','none');
%   hold on;
% 
%   plot([xDataMid-0.5*boxStdWidth;xDataMid+0.5*boxStdWidth],...
%        [             meanData; meanData            ],...
%        '-','Color',[1,1,1],'LineWidth',1.5);
%   hold on;  
%   
box off;

text( xDataMid, dataLabelYLocation, dataLabel,...
  'FontSize',8,'Interpreter','latex','HorizontalAlignment','center');  
hold on;

if(flag_firstCall ==1)

  if(isempty(xLabelText)==0)
    xlabel(xLabelText);
  end
  if(isempty(yLabelText)==0)
    ylabel(yLabelText);
  end

   
  titleFontSize = get(groot,'defaultAxesFontSize')...
                 *get(groot,'defaultAxesTitleFontSizeMultiplier');
  xTitle = axisLimits(1);% - 0.1*(axisLimits(2)-axisLimits(1));
  yTitle = axisLimits(4) + 0.1*(axisLimits(4)-axisLimits(3));
  text(xTitle,yTitle,titleText,'FontSize',titleFontSize,...
      'Interpreter','latex','HorizontalAlignment','left');
  hold on;
  %title(titleText);     
end  
end