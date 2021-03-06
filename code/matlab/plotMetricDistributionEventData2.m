function [figH] = plotMetricDistributionEventData2(figH, subPlotVec,...
  xPosition, metricData, metricScale, metricLabel, dataColors, axisLimits,...
  boxWidth, lineWidth, plotFontName, plotFontColor, ...
  flagPlotStart,flagPlotPhase,flagPlotEnd,...
  flag_motionSegmentEmphasis, flagEnableMinMaxLabels) 

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
          
          
          if(flag_motionSegmentEmphasis==0)
            figH = plotDistributionData2(...
                figH,subPlotVec, startPosition,...
                metricData.start,...
                metricScale,...
                dataColors,...
                boxWidth,...
                dataColors(1,:),...
                startEndWidth,...
                typeBoxAndWhisker2);          
          else
            plot( [startPosition,xPosition],...
                  [startMedian,startMedian].*metricScale,...
                 'Color',dataColorHalf,'LineWidth',lineWidth);
            hold on;
            plot( [startPosition],...
                  [startMedian].*metricScale,...
                 '.','Color',dataColorHalf,...
                 'MarkerSize',8,'MarkerFaceColor',dataColorHalf);
            hold on;
            plot( [startPosition],...
                  [startMedian].*metricScale,...
                 '.','Color',[1,1,1],...
                 'MarkerSize',4,'MarkerFaceColor',[1,1,1]);
            hold on;                     
          end

            
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
          
          if(flag_motionSegmentEmphasis == 2)
            figH = plotDistributionData2(...
                      figH,subPlotVec, endPosition,...
                      metricData.end,...
                      metricScale,...
                      dataColors,...
                      boxWidth,...
                      dataColors(1,:),...
                      startEndWidth,...
                      typeBoxAndWhisker2);                               
          else                      
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
          

          if(flag_motionSegmentEmphasis == 1)
            figH = plotDistributionData2(...
                figH,subPlotVec, xPosition,...
                metricData.phase,...
                metricScale,...
                dataColorHalf,...
                boxWidth*0.5,...
                dataColors(1,:),...
                startEndWidth,...
                typeBoxAndWhisker2);   
            
          else
          
            figH = plotDistributionData2(...
                figH,subPlotVec, xPosition,...
                metricData.phase,...
                metricScale,...
                dataColorHalf,...
                boxWidth*0.5,...
                dataColorHalf,...
                startEndWidth,...
                typeBoxAndWhisker3);   
          end
          
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
    
    dataList = {};
    
    switch(flag_motionSegmentEmphasis)
      case 0
        dataList = {'start','phase'};        
      case 1
        dataList = {'phase'};                
      case 2
        dataList = {'phase','end'};                
      otherwise 
        assert(0);
    end
      
    
    minData = Inf;
    maxData = -Inf;
    for z=1:1:length(dataList)
      if(metricData.(dataList{z}).min*metricScale < minData)
        minData = metricData.(dataList{z}).min*metricScale;
      end
      if(metricData.(dataList{z}).max*metricScale > maxData)
        maxData = metricData.(dataList{z}).max*metricScale;
      end
    end
    
    
    if(minData < axisLimits(1,3) && flagEnableMinMaxLabels==1)
      text(xPosition+0.25*boxWidth, axisLimits(1,3), ...
        {'$\downarrow$',sprintf('%1.1f',minData)},...
        'FontSize',6,...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color',dataColorHalf);
      hold on;
    end

    if(maxData > axisLimits(1,4) && flagEnableMinMaxLabels==1)
      text(xPosition+0.25*boxWidth, axisLimits(1,4), ...
        {sprintf('%1.1f',maxData),'$\uparrow$'},...
        'FontSize',6,...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color',dataColorHalf);
      hold on;              
    end    
    
    axisLimitsScaled = [axisLimits(1,1),...
                        axisLimits(1,2),...
                        axisLimits(1,3)*metricScale,...
                        axisLimits(1,4)*metricScale];

    axis(axisLimitsScaled);

    y0 = axisLimitsScaled(1,3);
    dy = axisLimitsScaled(1,4)-axisLimitsScaled(1,3);

    textYOffset = y0-0.075*dy;

    if(flag_motionSegmentEmphasis==0)
      xPosition=startPosition;
    end
    if(flag_motionSegmentEmphasis==2)
      xPosition=endPosition;
    end
    
    
    %'Interpreter','latex',
    text(xPosition,textYOffset,metricLabel,...
      'FontSize',8,'HorizontalAlignment','center',...
      'fontname',plotFontName,'Color',plotFontColor,...
      'fontweight','bold');
    hold on;
 end
