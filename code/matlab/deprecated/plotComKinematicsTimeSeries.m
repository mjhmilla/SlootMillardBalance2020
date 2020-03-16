function figH = plotComKinematicsTimeSeries(figH, subPlotVec, ...
                  c3dTime, movementSequence,...
                  c3dGrfFeet,  ...
                  comPosition,comVelocity, ...
                  gravityVec,...
                  lineColor,...
                  figureTitle,...
                  flag_ComVel0ComGpVsCop1)


                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
  
eV = -gravityVec./norm(gravityVec);  

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)

    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    timeBias = c3dTime(idx1);
  
    data = zeros(idx2-idx0+1,2);
    
    switch flag_ComVel0ComGpVsCop1
      case 0        
        for k = idx0:1:idx2
          v0C0 = comVelocity(k,:);
          data(k-idx0+1,1) = c3dTime(k,1)-timeBias;
          data(k-idx0+1,2) = norm(v0C0);
        end   
        plot( data(:,1), data(:,2).*100,...
              'Color',lineColor);
        hold on;
        xlabel('Time (s)');
        ylabel('Velocity (cm/s)'); 
      case 1
        for k = idx0:1:idx2
          r0C0 = comPosition(k,:);          
          rCP0 = c3dGrfFeet.cop(k,:)-r0C0;
          rCP0 = rCP0 - (rCP0*eV).*(eV'); %Ground projection
          data(k-idx0+1,1) = c3dTime(k,1)-timeBias;
          data(k-idx0+1,2) = norm(rCP0);
        end
        plot( data(:,1), data(:,2).*100,...
              'Color',lineColor);
        hold on;
        xlabel('Time (s)');
        ylabel('Distance (cm)'); 
      otherwise
        assert(0);
    end
    
    box off;
    title(figureTitle);       

  end

end