function [figH] = plotMetricDistributionEventData2(figH, subPlotVec,...
  xPosition, metricData, metricScale, metricLabel, dataColors, axisLimits,...
  boxWidth, lineWidth, plotFontName, plotFontColor, ...
  flagPlotStart,flagPlotPhase,flagPlotEnd) 

 if(isempty(metricData)==0)
   
    dataColorHalf =[1,1,1].*0.5 + dataColors(1,:).*0.5;    
   
    startEndWidth = boxWidth/2;
    typeBoxAndWhisker  = 0;
    typeDotAndWhisker  = 1;
    typeBoxAndWhisker2 = 2;
    typeBoxAndWhisker3 = 3;

    if( isfield(metricData,'start') && flagPlotStart==1)        
      if(isfield(metricData.start,'median'))
        if(isempty(metricData.start.median)==0)
          startPosition = xPosition-boxWidth;
          startMedian   = metricData.start.median;
          if(flagPlotPhase==1)
            plot( [startPosition,xPosition],...
                  [startMedian,startMedian].*metricScale,...
                 'Color',dataColorHalf,'LineWidth',lineWidth);
            hold on;
          end
          
          
        
          figH = plotDistributionData2(...
              figH,subPlotVec, startPosition,...
              metricData.start,...
              metricScale,...
              dataColors,...
              boxWidth,...
              dataColors(1,:),...
              startEndWidth,...
              typeBoxAndWhisker2);          
          
%           figH = plotDistributionData(...
%             figH,subPlotVec, startPosition,...
%             metricData.start,...
%             metricScale,...              
%             dataColors,...
%             startEndWidth,...
%             [1,1,1],...
%             lineWidth,...              
%             typeDotAndWhisker);
        end
      end
    end
    if( isfield(metricData,'end') && flagPlotEnd==1)
      if(isfield(metricData.end,'median'))
        if(isempty(metricData.end.median)==0)
          
      
          
          endPosition = xPosition+boxWidth;
          endMedian = metricData.end.median;
          if(flagPlotPhase==1)
            plot( [endPosition,xPosition],...
                  [endMedian,endMedian].*metricScale,...
                 'Color',dataColorHalf,'LineWidth',lineWidth);
            hold on;
            plot( [endPosition],...
                  [endMedian].*metricScale,...
                 '.','Color',dataColorHalf,...
                 'MarkerSize',8,'MarkerFaceColor',dataColorHalf);
            hold on;
            plot( [endPosition],...
                  [endMedian].*metricScale,...
                 '.','Color',[1,1,1],...
                 'MarkerSize',4,'MarkerFaceColor',[1,1,1]);
            hold on;
            
          end
%           figH = plotDistributionData2(...
%             figH,subPlotVec, endPosition,...
%             metricData.end,...
%             metricScale,...              
%             dataColorHalf,...
%             startEndWidth,...
%             dataColorHalf,...
%             lineWidth,...                   
%             typeDotAndWhisker); 

%         figH = plotDistributionData2(...
%             figH,subPlotVec, endPosition,...
%             metricData.end,...
%             metricScale,...
%             dataColorHalf,...
%             boxWidth*0.5,...
%             dataColorHalf,...
%             startEndWidth*0.5,...
%             typeBoxAndWhisker3);   
        end
      end
    end
    
    if( isfield(metricData,'phase') && flagPlotPhase==1)        
      if(isfield(metricData.phase,'median'))
        if(isempty(metricData.phase.median)==0)
          
        figH = plotDistributionData2(...
            figH,subPlotVec, xPosition,...
            metricData.phase,...
            metricScale,...
            dataColorHalf,...
            boxWidth*0.5,...
            dataColorHalf,...
            startEndWidth,...
            typeBoxAndWhisker3);   
          
          
%           figH = plotDistributionData2(...
%             figH,subPlotVec, xPosition,...
%             metricData.phase,...
%             metricScale,...
%             dataColorHalf,...
%             startEndWidth,...
%             dataColorHalf,...
%             lineWidth,...
%             typeDotAndWhisker);
          
%           figH = plotDistributionData(...
%             figH,subPlotVec, startPosition,...
%             metricData.start,...
%             metricScale,...              
%             dataColors,...
%             startEndWidth,...
%             [1,1,1],...
%             lineWidth,...              
%             typeDotAndWhisker);          
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

    if(flagPlotStart==1 && flagPlotPhase==0 && flagPlotEnd==0)
      xPosition=startPosition;
    end
    if(flagPlotStart==0 && flagPlotPhase==0 && flagPlotEnd==1)
      xPosition=endPosition;
    end
    
    
    %'Interpreter','latex',
    text(xPosition,textYOffset,metricLabel,...
      'FontSize',8,'HorizontalAlignment','center',...
      'fontname',plotFontName,'Color',plotFontColor,...
      'fontweight','bold');
    hold on;
 end
