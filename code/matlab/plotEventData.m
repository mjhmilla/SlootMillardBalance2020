function figH = plotEventData(...
                   figH, subPlotVec, ...
                   movementSequence, ...
                   xPoint, xPointScale,...                   
                   dataSeries, dataScale, ...
                   lineColor,  ...
                   dataLabel, dataLabelYOffset,...
                   xLabelText, ...
                   yLabelText, ...
                   titleText, ...
                   flag_start0Reference1End2)

                 
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
    xData = xPoint*xPointScale;

    yData =0;
    switch flag_start0Reference1End2
      case 0
        yData = dataSeries(movementSequence(z).indexStart, 1).*dataScale;
      case 1
        yData = dataSeries(movementSequence(z).indexReference, 1).*dataScale;
      case 2
        yData = dataSeries(movementSequence(z).indexEnd, 1).*dataScale;
      otherwise assert(0);
    end
      

    
    data = [data; xData, yData];
  end

end

if(isempty(data) == 0)
  plot( data(:,1),data(:,2),'-','Color',lineColor);    
  hold on;
  plot( data(:,1), data(:,2),...
        'o','Color',lineColor,'MarkerSize',3,'MarkerFaceColor',lineColor);
  hold on;

  text( data(1,1), (max(data(:,2))+dataLabelYOffset), dataLabel,...
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