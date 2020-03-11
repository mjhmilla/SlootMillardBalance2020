function figH = plotGrfZTimeSeries(figH, subPlotVec, ...
                   c3dTime, movementSequence, c3dGrfChair,...
                   idSittingDynamic, idCrouchingStable,...
                   lineColor,...
                   figureTitle,...
                   flag_zeroAtSeatOff)

if(flag_zeroAtSeatOff==1)
  assert(idSittingDynamic==1);
  assert(idCrouchingStable==2);
end
                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  
                 
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).phaseTransitions))==0)

    idx0 = movementSequence(z).phaseTransitionIndices(1);
    timeBias = c3dTime(idx0);
    
    if(flag_zeroAtSeatOff==1)
      flag_alertMultipleSeatOffs=1;
      idx1 = getLastSeatOffIndex( movementSequence(z).phaseTransitions,...
                          movementSequence(z).phaseTransitionIndices,...
                          idSittingDynamic,idCrouchingStable,...
                          flag_alertMultipleSeatOffs);
        
      timeBias = c3dTime(idx1);
    end


    idx2 = movementSequence(z).phaseTransitionIndices(end);

    


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