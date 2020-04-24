function figH = plotGroupComparisons( figH, subPlotVec,...
                              xPositionGroupA, statsStructA,...
                              xPositionGroupB, statsStructB,...
                              statsABPhase, statsABStart, statsABEnd,...
                              phaseLabel,startLabel,endLabel, ...
                              flag_drawAnnotationLines,boxWidth,...
                              plotFontName, plotFontColor)

figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end

axisLim = axis;

xMin = axisLim(1,1);
xMax = axisLim(1,2);
yMin = axisLim(1,3);
yMax = axisLim(1,4);

yDelta  = 0.05*(yMax-yMin);
xDelta  = yDelta;

yBottom = yMin + 0.05*(yMax-yMin);
yTop    = yMax - 0.05*(yMax-yMin);
yMiddle = 0.5*(yTop+yBottom);

startEndWidth = boxWidth;

xP = max(xPositionGroupA,xPositionGroupB) + boxWidth*2 ;
printOrder = 3;

%%
% Add the start p values
%%
if(isempty(statsABStart)==0)
  if( isfield(statsStructA,'start') && isfield(statsStructB,'start'))        
    if(isfield(statsStructA.start,'median') && isfield(statsStructB.start,'median'))
      if(isempty(statsStructA.start.median)==0 && isempty(statsStructB.start.median)==0  )

        xA = xPositionGroupA - startEndWidth;
        xB = xPositionGroupB - startEndWidth;
        yA = statsStructA.start.p75 + yDelta;
        yB = statsStructB.start.p75 + yDelta;

        xLeft = min(xA,xB);

        if(flag_drawAnnotationLines==1)
          plot([xA;xA],[yTop;yA],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
          plot([xB;xB],[yTop;yB],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
          plot([xLeft;xP],[yTop;yTop],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
        end

        starText = '';
        if(statsABStart.h==1)
          starText = '*';
        end

        statsText  = '';
        if(statsABStart.p*(10^(printOrder)) >= 1)
          statsText = sprintf(['%s:\n%sp = %1.',num2str(printOrder),'f'],...
                                startLabel,starText,statsABStart.p);
        else
          statsText = sprintf('%s:\n%sp = %1.2e',startLabel,starText,statsABStart.p);
        end

        text( xP,yTop,statsText,...
            'FontSize',6,'HorizontalAlignment','left',...
            'VerticalAlignment','cap',...
            'fontname',plotFontName);
          hold on;   

      end
    end
  end
end
%%
% Add the end p values
%%
if(isempty(statsABEnd)==0)
  if( isfield(statsStructA,'end') && isfield(statsStructB,'end'))        
    if(isfield(statsStructA.end,'median') && isfield(statsStructB.end,'median'))
      if(isempty(statsStructA.end.median)==0 && isempty(statsStructB.end.median)==0  )

        xA = xPositionGroupA + startEndWidth;
        xB = xPositionGroupB + startEndWidth;
        yA = statsStructA.end.p25 - yDelta;
        yB = statsStructB.end.p25 - yDelta;

        xLeft = min(xA,xB);

        if(flag_drawAnnotationLines==1)
          plot([xA;xA],[yA;yBottom],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
          plot([xB;xB],[yB;yBottom],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
          plot([xLeft;xP],[yBottom;yBottom],'-','Color',[0,0,0],'LineWidth',0.5);
          hold on;
        end

        starText = '';
        if(statsABEnd.h==1)
          starText = '*';
        end

        statsText  = '';
        if(statsABEnd.p*(10^printOrder) >= 1)
          statsText = sprintf(['%s:\n%sp = %1.',num2str(printOrder),'f'],...
                                        endLabel,starText,statsABEnd.p);
        else
          statsText = sprintf('%s:\n%sp = %1.2e',endLabel,starText,statsABEnd.p);
        end

        text( xP,yBottom,statsText,...
            'FontSize',6,'HorizontalAlignment','left',...
            'VerticalAlignment','bottom',...
            'fontname',plotFontName);
          hold on;   
      end
    end
  end
end

%%
% Add the whole movement p values
%%
if(isempty(statsABPhase)==0)
  if( isfield(statsStructA,'phase') && isfield(statsStructB,'phase'))        
    if(isfield(statsStructA.phase,'median') && isfield(statsStructB.phase,'median'))
      if(isempty(statsStructA.phase.median)==0 && isempty(statsStructB.phase.median)==0  )


        xA = xPositionGroupA + boxWidth*1.5;
        xB = xPositionGroupB + boxWidth*1.5;

        xData     = xA;
        yMedian   = statsStructA.phase.median; 
        if(xB > xA)
          xData     = xB;
          yMedian   = statsStructB.phase.median;         
        end
        
        %plot([xData;xP],[yMedian;yMiddle],'-','Color',[0,0,0],'LineWidth',0.5);
        %hold on;

        starText = '';
        if(statsABPhase.h==1)
          starText = '*';
        end

        statsText  = '';
        if(statsABPhase.p*(10^printOrder) >= 1)
          statsText = sprintf(['%s:\n%sp = %1.',num2str(printOrder),'f'],...
                                        phaseLabel,starText,statsABPhase.p);
        else
          statsText = sprintf('%s:\n%sp = %1.1e',phaseLabel,starText,statsABPhase.p);
        end

        vertAlign = 'cap';
        if(yMiddle > yMedian)
          vertAlign = 'baseline';
        end

        text( xP,yMiddle,statsText,...
            'FontSize',6,'HorizontalAlignment','left',...
            'VerticalAlignment', vertAlign,...
            'fontname',plotFontName,'Color',plotFontColor);
          hold on;   

      end
    end
  end
end
