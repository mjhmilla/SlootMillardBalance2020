function figH = plotDistributionData(figH, subPlotVec,...
                  xPosition, statsStruct, ...
                  dataColor, dataBodyWidth, dataBodyColor, ...                  
                  flag_plotType)
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end
                



%Phase
plot([xPosition;xPosition],...
     [statsStruct.min;statsStruct.max],...
     '-','Color',dataColor,'LineWidth',0.5);
hold on;   
plot([xPosition],...
     [statsStruct.min],...
     '.','Color',dataColor,'MarkerSize',3);
hold on;   
plot([xPosition],...
     [statsStruct.max],...
     '.','Color',dataColor,'MarkerSize',3);
hold on;   

switch flag_plotType
  case 0
    dataBoxMin=statsStruct.p25;
    dataBoxMax=statsStruct.p75;
    dataBox = [ xPosition+dataBodyWidth*0.5, dataBoxMin;...
               xPosition+dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMin;... 
               xPosition+dataBodyWidth*0.5,  dataBoxMin];

    normWhite = norm([1,1,1]-dataBodyColor);
    medianLineColor = [1,1,1];
    
    if(normWhite < 1e-6)
      fill(dataBox(:,1),dataBox(:,2),dataBodyColor,'EdgeColor',dataColor);
      hold on;           
      
      medianLineColor = dataColor;
    else
      fill(dataBox(:,1),dataBox(:,2),dataBodyColor,'EdgeColor','none');
      hold on;           
      
      medianLineColor = [1,1,1];
    end
             
    
    plot([    xPosition-0.5*dataBodyWidth;xPosition+0.5*dataBodyWidth],...
         [statsStruct.median; statsStruct.median],...
         '-','Color',[1,1,1],'LineWidth',0.5);
    hold on;  
  case 1
    plot([xPosition],...
     [statsStruct.median],...
     'o','Color',dataColor,'MarkerSize',5,'MarkerFaceColor',dataBodyColor);
    hold on;    
  otherwise assert(0)
end






