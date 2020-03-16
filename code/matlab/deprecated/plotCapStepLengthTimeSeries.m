function figH = plotCapStepLengthTimeSeries(figH, subPlotVec, ...
                  c3dTime, movementSequence,...
                  c3dGrfFeet,...
                  capData,...                   
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
  
           

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)

    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    timeBias = c3dTime(idx1);
    idx2 = movementSequence(z).indexEnd;

    distance =[];
    width = [];
    
    switch flag_BalancePointsVsCom0VsCop1
      case 0 
        distance = zeros( (idx2-idx0+1),1);
        width    = zeros( (idx2-idx0+1),1);
        for k = idx0:1:idx2
          eU = capData.u(k,:)';
          eK = capData.k(k,:)';
          eN = getCrossProductMatrix(eK)*eU;
          rCF0  = capData.r0F0(k,:)-capData.r0G0(k,:);
          lCF0  = rCF0*eU;          
          wCF0  = rCF0*eN;
          distance(k-idx0+1,1) = lCF0;
          width(k-idx0+1,1)    =  wCF0;
        end         
        
      case 1
        distance = zeros( (idx2-idx0+1),1);
        width    = zeros( (idx2-idx0+1),1);
        for k = idx0:1:idx2
          eU = capData.u(k,:)';
          eK = capData.k(k,:)';
          eN = getCrossProductMatrix(eK)*eU;
          rCF0  = capData.r0F0(k,:)-c3dGrfFeet.cop(k,:);
          wCF0  = rCF0*eN;
          lCF0  = rCF0*eU;
          distance(k-idx0+1,1) = lCF0;
          width(k-idx0+1,1)    =  wCF0;
        end 
      otherwise
        assert(0);
    end    
    
    switch flag_AnalyzeBalanceAlong0Across1
      
      case 0
        plot( c3dTime(idx0:1:idx2,1)-timeBias, distance.*100,...
              'Color',lineColor);
        hold on;
        xlabel('Time (s)');
        ylabel('Length (cm)');        
      case 1
        plot( c3dTime(idx0:1:idx2,1)-timeBias, width.*100,...
              'Color',lineColor);
        hold on;
        xlabel('Time (s)');
        ylabel('Width (cm)');            
      otherwise
        assert(0)
        
    end
    box off;
    title(figureTitle); 

  end

end