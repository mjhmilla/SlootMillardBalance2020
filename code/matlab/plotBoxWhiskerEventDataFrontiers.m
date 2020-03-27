function [figH, dataSummary, dataRaw] = plotBoxWhiskerEventDataFrontiers(...
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
                                 flag_drawBox,...
                                 flag_firstCall,...
                                 numberOfPhases)

                 
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
if(flag_firstCall ==1)
  plot([-10,2000],[0,0],'Color',[1,1,1].*0.75);
  hold on;
end

  
%Two phases: 
% 1. start-to-seatoff
% 2. seatoff-to-standing
%numberOfPhases = 2;

data(numberOfPhases) = struct('y',[]);

dataSummary(numberOfPhases) =...
  struct('min',NaN,'p25',NaN,'mean',NaN,'p75',NaN','max',NaN,...
         'start',zeros(1,1).*NaN,...
         'end',zeros(1,1).*NaN,...
         'duration',NaN);
  





eventStart(numberOfPhases) = struct('x',[],'y',[]);
eventEnd(numberOfPhases)   = struct('x',[],'y',[]);



dataMin(numberOfPhases) = struct('x',[],'y',[]);
dataMax(numberOfPhases) = struct('x',[],'y',[]);

maxNumberOfRepetitions = 0;
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
    maxNumberOfRepetitions=maxNumberOfRepetitions+1;
  end  
end

                 
dataRaw(numberOfPhases,maxNumberOfRepetitions) = struct('x',[],'y',[]);

indexRepetition=0;
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
    indexRepetition =  indexRepetition+1;
    
    
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;
    
    
    for p=1:1:numberOfPhases

      switch p
        case 1
          idxA = idx0;
          idxB = idx1;
        case 2
          idxA = idx1;
          idxB = idx2;          
        otherwise assert(0);          
      end
       
      eventStart(p).x = [eventStart(p).x; dataIndex(idxA,1)];
      eventStart(p).y = [eventStart(p).y; dataSeries(idxA,1).*dataScale];
                  
      eventEnd(p).x = [eventEnd(p).x; dataIndex(idxB,1)];
      eventEnd(p).y = [eventEnd(p).y; dataSeries(idxB,1).*dataScale];
                    
                   
      data(p).y   = [data(p).y; dataSeries(idxA:idxB,:).*dataScale];
    
      
      [valMin,idxMin] = min(dataSeries(idxA:idxB,1).*dataScale);
      [valMax,idxMax] = max(dataSeries(idxA:idxB,1).*dataScale);
      
      dataMin(p).x = [dataMin(p).x; dataIndex(idxA+idxMin-1,1)];
      dataMin(p).y = [dataMin(p).y; valMin];

      dataMax(p).x = [dataMax(p).x; dataIndex(idxA+idxMax-1,1)];
      dataMax(p).y = [dataMax(p).y; valMax];
      
        
      
      dataRaw(p,indexRepetition).x = ...
        [dataRaw(p,indexRepetition).x; dataIndex(idxA:idxB,1)];
      dataRaw(p,indexRepetition).y = ...
        [dataRaw(p,indexRepetition).y; dataSeries(idxA:idxB,1).*dataScale];
      
    end
    
  end

end



if(isempty(data) == 0)

  for p=1:1:numberOfPhases

    dataSorted = sort(data(p).y(:,1));
    n = length(dataSorted);

    n25 = round(n*0.25);
    n75 = round(n*0.75);

    dataSummary(p).min = mean(dataMin(p).y);
    dataSummary(p).mean= mean(data(p).y(:,1));
    dataSummary(p).max = mean(dataMax(p).y);

    dataSummary(p).p25 = dataSorted(n25,1);
    dataSummary(p).p75 = dataSorted(n75,1);
    
    dataSummary(p).startMin = min(eventStart(p).y);
    dataSummary(p).startMean= mean(eventStart(p).y);    
    dataSummary(p).startMax = max(eventStart(p).y);

    dataSummary(p).start= dataSummary(p).startMean;
    
    dataSummary(p).endMin   = min(eventEnd(p).y);    
    dataSummary(p).endMean  = mean(eventEnd(p).y);
    dataSummary(p).endMax   = max(eventEnd(p).y);
    dataSummary(p).end= dataSummary(p).endMean;
    
    dataSummary(p).durationMin  = min(eventEnd(p).x-eventStart(p).x);    
    dataSummary(p).durationMean = mean(eventEnd(p).x-eventStart(p).x);   
    dataSummary(p).durationMax  = max(eventEnd(p).x-eventStart(p).x);
  
    dataSummary(p).duration = dataSummary(p).durationMean;

    
  end
  
  xDataMid = (xPoint*xPointScale);
  p=2;
  
  figH = plotAntsOnALog(figH,subPlotVec, ...
            xDataMid,...
            dataSummary(p).min, dataSummary(p).mean, dataSummary(p).max,...
            dataSummary(p).p25, dataSummary(p).p75, ...
            [dataSummary(p).start,dataSummary(p).end],...
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

    set(gca,'TickLength',[0 0])
    
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