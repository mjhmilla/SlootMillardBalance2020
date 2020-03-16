function figH = plotFpeCapErrorTimeSeries(figH, subPlotVec, ...
                   c3dTime, movementSequence, ...
                   fpeData, capData,...
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
    idx2 = movementSequence(z).indexEnd;

    timeBias = c3dTime(idx1);

    
    rErr = fpeData.r0F0(idx0:1:idx2,:) ...
          -capData.r0F0(idx0:1:idx2,:);
    dFpeCap = sum(rErr.^2,2).^0.5;

    plot( c3dTime(idx0:1:idx2,1)-timeBias, dFpeCap.*100,...
          'Color',lineColor);
    hold on;
    xlabel('Time (s)');
    ylabel('Distance (cm)');
    box off;
    title(figureTitle);
   
  end

end