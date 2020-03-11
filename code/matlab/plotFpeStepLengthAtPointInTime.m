function figH = plotFpeStepLengthAtPointInTime(figH, subPlotVec, ...
                   subjectNumber, subjectLabel, movementSequence, fpeData,...
                   idSittingDynamic, idCrouchingStable,...
                   lineColor,...
                   figureTitle,...
                   flag_pointInTimeSelector)

assert(flag_pointInTimeSelector==1); %For now: this is the only option                 
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
                          flag_alertMultipleSeatOffs);
        
    end



    

    rGF0 = fpeData.r0F0(idxSample,:) ...
          -fpeData.r0G0(idxSample,:);
    data(z,1) = subjectNumber;
    data(z,2) = sum(rGF0.^2,2).^0.5;
  end

end

plot( data(:,1), data(:,2).*100,...
      '-','Color',lineColor);
hold on;
plot( data(:,1), data(:,2).*100,...
      'o','Color',lineColor,'MarkerSize',3,'MarkerFaceColor',lineColor);
hold on;
%axisLim = axis;
%axisLim(1) = 0;
%axisLim(2) = data(1,1);
%axisLim(3) = 0;
%axis(axisLim);

text( data(1,1), (max(data(:,2).*100)+0.5), subjectLabel,...
  'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');

hold on;

%xlabel('Subject No');
ylabel('Distance (cm)');
box off;
title(figureTitle);





