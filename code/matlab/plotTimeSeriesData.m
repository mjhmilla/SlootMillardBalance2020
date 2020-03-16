function figH = plotTimeSeriesData(...
                   figH, subPlotVec, ...
                   movementSequence, ...
                   timeSeries, timeScale, ...
                   dataSeries, dataScale, ...
                   lineColor,  ...
                   xLabelText, ...
                   yLabelText, ...
                   titleText)

                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  
                 
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)

    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    timeReference = 0;
    if(isnan(idx1) == 0)
      timeReference = timeSeries(idx1);
    end
    idx2 = movementSequence(z).indexEnd;    

    plot( (timeSeries(idx0:1:idx2,1)-timeReference).*timeScale,...
          dataSeries(idx0:1:idx2,1).*dataScale,...
          'Color',lineColor);
    hold on;
    xlabel(xLabelText);
    ylabel(yLabelText);
    box off;
    title(titleText);     

  end

end