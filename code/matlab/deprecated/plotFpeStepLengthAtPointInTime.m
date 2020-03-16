function figH = plotFpeStepLengthAtPointInTime(figH, subPlotVec, ...
                  subjectNumber, subjectLabel,...
                  movementSequence, ...
                  c3dGrfFeet,...
                  fpeData, gravityVec,...
                  lineColor,...
                  figureTitle,...
                  flag_BalancePointsVsCom0VsCop1,...
                  flag_AnalyzeBalanceAlong0Across1)


                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  

data = zeros(length(movementSequence),2);
eV = -gravityVec./norm(gravityVec);  
flag_data=0;

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)
    flag_data=1;
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    idxSample = idx1;

    distance =[];
    width = [];    
    switch flag_BalancePointsVsCom0VsCop1
      case 0 
        distance = zeros( 1, 1);
        width    = zeros( 1, 1);
        %for k = idx0:1:idx2
          eN = fpeData.n(idxSample,:)';
          eU = getCrossProductMatrix(eN)*eV;
          rCF0  = fpeData.r0F0(idxSample,:)-fpeData.r0G0(idxSample,:);
          wCF0  = rCF0*eN;
          lCF0  = rCF0*eU;
          distance = lCF0;
          width    =  wCF0;
        %end         
        
      case 1
        distance = zeros(1,1);
        width    = zeros(1,1);
        %for k = idx0:1:idx2
          eN = fpeData.n(idxSample,:)';
          eU = getCrossProductMatrix(eN)*eV;
          rCF0  = fpeData.r0F0(idxSample,:)-c3dGrfFeet.cop(idxSample,:);
          wCF0  = rCF0*eN;
          lCF0  = rCF0*eU;
          distance = lCF0;
          width    =  wCF0;
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

if(flag_data==1)

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

  %xlabel('Subject No');

  box off;
  title(figureTitle);

end


