function figH = plotComKinematicsAtPointInTime(...
                  figH, subPlotVec, ...
                  subjectNumber, subjectLabel,...
                  c3dTime, movementSequence,...
                  c3dGrfFeet,  ...
                  comPosition,comVelocity,...
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


data = zeros(length(movementSequence),2);
flag_data= 0;
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)
    flag_data=1;
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    data(z,1)=subjectNumber;
    
    switch flag_ComVel0ComGpVsCop1
      case 0        
        v0C0 = comVelocity(idx1,:);
        data(z,2) = norm(v0C0); 
      case 1
        r0C0 = comPosition(idx1,:);          
        rCP0 = c3dGrfFeet.cop(idx1,:)-r0C0;
        rCP0 = rCP0 - (rCP0*eV).*(eV'); %Ground projection
        data(z,2) = norm(rCP0);
      otherwise
        assert(0);
    end
    

  end

end

if(flag_data==1)
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
end