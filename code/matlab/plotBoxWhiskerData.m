function figH = plotBoxWhiskerData(...
                   figH, subPlotVec, ...
                   movementSequence, ...
                   xPoint, xPointScale,...                   
                   dataSeries, dataScale, ...
                   lineColor,  ...
                   boxStdWidth,...
                   dataLabel, dataLabelYOffset,...
                   xLabelText, ...
                   yLabelText, ...
                   titleText,...
                   flag_IntervalMode)

                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  
data = [];
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
 
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    idxA = 0;
    idxB = 0;
    switch(flag_IntervalMode)
      case 0
        idxA = idx0;
        idxB = idx2;
      case 1
        idxA = idx1;
        idxB = idx2;
      otherwise assert(0);
      end

    yData = dataSeries(idxA:idxB,:).*dataScale;
    
    data = [data; ones(size(yData)).*(xPoint*xPointScale), yData];
  end

end

if(isempty(data) == 0)

  minData  = min(data(:,2));
  maxData  = max(data(:,2));
  meanData = mean(data(:,2));
  stdData  = std(data(:,2));

  xDataMid = xPoint*xPointScale;

  boxStd = [ xDataMid+boxStdWidth*0.5, meanData-stdData;...
             xDataMid+boxStdWidth*0.5, meanData+stdData;... 
             xDataMid-boxStdWidth*0.5, meanData+stdData;... 
             xDataMid-boxStdWidth*0.5, meanData-stdData;... 
             xDataMid+boxStdWidth*0.5, meanData-stdData];

  plot([xDataMid;xDataMid],[minData;maxData],'-','Color',lineColor);
  hold on;

  fill(boxStd(:,1),boxStd(:,2),[1,1,1],'EdgeColor',lineColor);
  hold on;

  plot([xDataMid-0.5*boxStdWidth;xDataMid+0.5*boxStdWidth],...
       [             meanData; meanData            ],...
       '-','Color',lineColor);
  hold on;

  plot([xDataMid-0.5*boxStdWidth; xDataMid+0.5*boxStdWidth],...
       [              maxData; maxData              ],...
       '-','Color',lineColor);
  hold on;

  plot([xDataMid-0.5*boxStdWidth; xDataMid+0.5*boxStdWidth],...
       [              minData; minData              ],...
       '-','Color',lineColor);
  hold on;


  text( xDataMid, maxData+dataLabelYOffset, dataLabel,...
    'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');  
  hold on;
  if(isempty(xLabelText)==0)
    xlabel(xLabelText);
  end
  if(isempty(yLabelText)==0)
    ylabel(yLabelText);
  end
  box off;
  title(titleText);     
end