function figH = plotGrfZTimeSeries(figH, subPlotVec, ...
                   c3dTime, movementSequence, ...
                   c3dGrfChair,...
                   lineColor,...
                   figureTitle)

                 
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
    timeBias = c3dTime(idx1);
    idx2 = movementSequence(z).indexEnd;    
    

    plot( c3dTime(idx0:1:idx2,1)-timeBias,...
          c3dGrfChair.force(idx0:1:idx2,3),...
          'Color',lineColor);
    hold on;
    xlabel('Time (s)');
    ylabel('Force (N)');
    box off;
    title(figureTitle);     

  end

end