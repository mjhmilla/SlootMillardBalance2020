function figH = plotCapStepLengthAtPointInTime(figH, subPlotVec, ...
                   subjectNumber, subjectLabel, movementSequence, ...
                   c3dGrfChair,c3dGrfFeet,capData,...
                   idSittingDynamic, idCrouchingStable,...
                   lineColor,...
                   figureTitle,...
                   flag_pointInTimeSelector,...
                   flag_BalancePointsVsCom0VsCop1,...
                   flag_AnalyzeBalanceAlong0Across1)

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
                          flag_alertMultipleSeatOffs,c3dGrfChair);
        
    end



    distance =[];
    width = [];

    
    switch flag_BalancePointsVsCom0VsCop1
      case 0 
        distance = zeros( 1,1);
        width    = zeros( 1,1);
        %for k = idx0:1:idx2
          eU = capData.u(idxSample,:)';
          eK = capData.k(idxSample,:)';
          eN = getCrossProductMatrix(eK)*eU;
          rCF0  = capData.r0F0(idxSample,:)-capData.r0G0(idxSample,:);
          lCF0  = rCF0*eU;          
          wCF0  = rCF0*eN;
          distance(1,1) = lCF0;
          width(1,1)    =  wCF0;
        %end         
        
      case 1
        distance = zeros( 1,1);
        width    = zeros( 1,1);
        %for k = idx0:1:idx2
          eU = capData.u(idxSample,:)';
          eK = capData.k(idxSample,:)';
          eN = getCrossProductMatrix(eK)*eU;
          rCF0  = capData.r0F0(idxSample,:)-c3dGrfFeet.cop(idxSample,:);
          wCF0  = rCF0*eN;
          lCF0  = rCF0*eU;
          distance(1,1) = lCF0;
          width(1,1)    =  wCF0;
        %end 
      otherwise
        assert(0);
    end    
    
    data(z,1) = subjectNumber;
    switch flag_AnalyzeBalanceAlong0Across1
      case 0
        data(z,2) = distance;
      case 1
        data(z,2) = width;        
      otherwise
        assert(0)        
    end  

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

switch flag_AnalyzeBalanceAlong0Across1
  case 0
    ylabel('Length (cm)');
  case 1
    ylabel('Width (cm)'); 
  otherwise
    assert(0)        
end

box off;
title(figureTitle);





