function figH = plotGrfZAtPointInTime(figH, subPlotVec, ...
                   subjectNumber, subjectLabel,...
                   c3dTime, movementSequence, c3dGrfChair,...
                   idSittingDynamic, idCrouchingStable,...
                   lineColor,...
                   figureTitle,...
                   flag_pointInTimeSelector)

if(flag_pointInTimeSelector==1)
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

data = zeros(length(movementSequence),2);


for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).phaseTransitions))==0)

    idx0 = movementSequence(z).phaseTransitionIndices(1);
    
    idxSample = 0;
    
    if(flag_pointInTimeSelector==1)
      flag_alertMultipleSeatOffs=1;
      idxSample = getLastSeatOffIndex( movementSequence(z).phaseTransitions,...
                          movementSequence(z).phaseTransitionIndices,...
                          idSittingDynamic,idCrouchingStable,...
                          flag_alertMultipleSeatOffs, c3dGrfChair);
        
    end


    idx2 = movementSequence(z).phaseTransitionIndices(end);

    
    data(z,1)=subjectNumber;
    data(z,2) = c3dGrfChair.force(idxSample,3);
   
  
  end

end

plot( data(:,1),...
      data(:,2),...
      '-','Color',lineColor);    
hold on;
plot( data(:,1), data(:,2),...
      'o','Color',lineColor,'MarkerSize',3,'MarkerFaceColor',lineColor);
hold on;

text( data(1,1), (max(data(:,2))+2.5), subjectLabel,...
  'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');

ylabel('Force (N)');
box off;
title(figureTitle);  