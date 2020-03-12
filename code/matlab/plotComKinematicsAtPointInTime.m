function figH = plotComKinematicsAtPointInTime(figH, subPlotVec, ...
                   subjectNumber, subjectLabel,...
                   c3dTime, movementSequence, c3dGrfChair,c3dGrfFeet,  ...
                   comPosition,comVelocity,gravityVec,...
                   idSittingDynamic, idCrouchingStable,...
                   lineColor,...
                   figureTitle,...
                   flag_pointInTimeSelector,...
                   flag_ComVel0ComGpVsCop1)

assert(flag_pointInTimeSelector==1);                 
assert(idSittingDynamic==1);
assert(idCrouchingStable==2);
                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  
eV = -gravityVec./norm(gravityVec);  


data = zeros(length(movementSequence),2);

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).phaseTransitions))==0)

    idx0 = movementSequence(z).phaseTransitionIndices(1);
    timeBias = c3dTime(idx0);
    
    idxSample = 0;
    
    if(flag_pointInTimeSelector==1)
      flag_alertMultipleSeatOffs=1;
      idxSample = getLastSeatOffIndex( ...
                          movementSequence(z).phaseTransitions,...
                          movementSequence(z).phaseTransitionIndices,...
                          idSittingDynamic,idCrouchingStable,...
                          flag_alertMultipleSeatOffs,c3dGrfChair);
        
    end


    idx2 = movementSequence(z).phaseTransitionIndices(end);
    data(z,1)=subjectNumber;
    
    switch flag_ComVel0ComGpVsCop1
      case 0        
        v0C0 = comVelocity(idxSample,:);
        data(z,2) = norm(v0C0); 
      case 1
        r0C0 = comPosition(idxSample,:);          
        rCP0 = c3dGrfFeet.cop(idxSample,:)-r0C0;
        rCP0 = rCP0 - (rCP0*eV).*(eV'); %Ground projection
        data(z,2) = norm(rCP0);
      otherwise
        assert(0);
    end
    

  end

end

plot( data(:,1),...
      data(:,2).*100,...
      '-','Color',lineColor);    
hold on;
plot( data(:,1), data(:,2).*100,...
      'o','Color',lineColor,'MarkerSize',3,'MarkerFaceColor',lineColor);
hold on;
text( data(1,1), (max(data(:,2).*100)+1), subjectLabel,...
  'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');

switch flag_ComVel0ComGpVsCop1
  case 0            
    ylabel('Velocity (cm/s)'); 
  case 1
    ylabel('Distance (cm)'); 
    
  otherwise
    assert(0);
end


box off;
title(figureTitle);  
