function figH = plotMetricDistributionEventData(figH, subPlotVec,...
  xPosition, metricData, metricScale, metricLabel, dataColors, axisLimits,...
  boxWidth, lineWidth, plotFontName, plotFontColor) 

 if(isempty(metricData)==0)
   

    startEndWidth = boxWidth/2;
    typeBoxAndWhisker=0;
    typeDotAndWhisker=1;


    if( isfield(metricData,'start'))        
      if(isfield(metricData.start,'median'))
        if(isempty(metricData.start.median)==0)
          startPosition = xPosition-boxWidth;
          startMedian   = metricData.start.median;
          plot( [startPosition,xPosition],...
                [startMedian,startMedian].*metricScale,...
               'Color',dataColors(2,:),'LineWidth',lineWidth);
          hold on;

          figH = plotDistributionData(...
            figH,subPlotVec, startPosition,...
            metricData.start,...
            metricScale,...              
            dataColors,...
            startEndWidth,...
            [1,1,1],...
            lineWidth,...              
            typeDotAndWhisker);
        end
      end
    end
    if( isfield(metricData,'end'))
      if(isfield(metricData.end,'median'))
        if(isempty(metricData.end.median)==0)
          endPosition = xPosition+boxWidth;
          endMedian = metricData.end.median;
          plot( [endPosition,xPosition],...
                [endMedian,endMedian].*metricScale,...
               'Color',dataColors(2,:),'LineWidth',lineWidth);
          hold on;

          figH = plotDistributionData(...
            figH,subPlotVec, endPosition,...
            metricData.end,...
            metricScale,...              
            dataColors,...
            startEndWidth,...
            dataColors(1,:),...
            lineWidth,...                   
            typeDotAndWhisker); 
        end
      end
    end
    
    if( isfield(metricData,'phase'))        
      if(isfield(metricData.phase,'median'))
        if(isempty(metricData.phase.median)==0)      
          figH = plotDistributionData(...
            figH,subPlotVec, xPosition,...
            metricData.phase,...
            metricScale,...
            dataColors,...
            boxWidth,...
            dataColors(1,:),...
            lineWidth,...
            typeBoxAndWhisker);
        end
      end
    end
    
    axisLimitsScaled = [axisLimits(1,1),...
                        axisLimits(1,2),...
                        axisLimits(1,3)*metricScale,...
                        axisLimits(1,4)*metricScale];

    axis(axisLimitsScaled);

    y0 = axisLimitsScaled(1,3);
    dy = axisLimitsScaled(1,4)-axisLimitsScaled(1,3);

    textYOffset = y0-0.075*dy;

    %'Interpreter','latex',
    text(xPosition,textYOffset,metricLabel,...
      'FontSize',8,'HorizontalAlignment','center',...
      'fontname',plotFontName,'Color',plotFontColor);
    hold on;
 end
