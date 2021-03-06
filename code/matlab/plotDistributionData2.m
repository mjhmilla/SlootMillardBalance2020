function figH = plotDistributionData2(figH, subPlotVec,...
                  xPosition, statsStruct, dataScale, ...
                  dataColor, dataBodyWidth, dataBodyColor, ...
                  lineWidth,...
                  flag_plotType)
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end


normWhite = norm([1,1,1]-dataBodyColor);
boxEdgeColor = [1,1,1];

if(normWhite < 1e-6)
  boxEdgeColor = dataColor(1,:);
else
  boxEdgeColor = [1,1,1];
end


switch flag_plotType
  case 0

    plot([xPosition;xPosition],...
        [statsStruct.min;statsStruct.max].*dataScale,...
        '-','Color',dataColor(2,:),'LineWidth',lineWidth);
    hold on;   
    plot([xPosition],...
         [statsStruct.min].*dataScale,...
         '.','Color',dataColor(2,:),'MarkerSize',3);
    hold on;   
    plot([xPosition],...
         [statsStruct.max].*dataScale,...
         '.','Color',dataColor(2,:),'MarkerSize',3);
    hold on;  
    
    dataBoxMin=statsStruct.p05;
    dataBoxMax=statsStruct.p95;
    dataBox = [ xPosition+dataBodyWidth*0.5, dataBoxMin;...
               xPosition+dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMin;... 
               xPosition+dataBodyWidth*0.5,  dataBoxMin];
    fill(dataBox(:,1),dataBox(:,2).*dataScale,dataColor(2,:),...
        'EdgeColor',boxEdgeColor);
    hold on;           
  case 1
    plot([xPosition;xPosition],...
        [statsStruct.p05;statsStruct.p95].*dataScale,...
        '-','Color',dataColor(1,:),'LineWidth',lineWidth);
    hold on;   
    plot([xPosition],...
         [statsStruct.p05].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;   
    plot([xPosition],...
         [statsStruct.p95].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;    
    
  case 2
    plot([xPosition;xPosition],...
        [statsStruct.min;statsStruct.max].*dataScale,...
        '-','Color',dataColor(1,:),'LineWidth',lineWidth);
    hold on;   
    plot([xPosition],...
         [statsStruct.min].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;   
    plot([xPosition],...
         [statsStruct.max].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;    
    
  case 3
    plot([xPosition;xPosition],...
        [statsStruct.min;statsStruct.max].*dataScale,...
        '-','Color',dataColor(1,:),'LineWidth',lineWidth);
    hold on;   
    plot([xPosition],...
         [statsStruct.min].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;   
    plot([xPosition],...
         [statsStruct.max].*dataScale,...
         '.','Color',dataColor(1,:),'MarkerSize',3);
    hold on;    
    
                        
    
  otherwise assert(0);
end
         
%Phase




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
    boxEdgeColor = [1,1,1];
    
    if(normWhite < 1e-6)
      fill(dataBox(:,1),dataBox(:,2).*dataScale,dataBodyColor,...
          'EdgeColor',dataColor(1,:),'LineWidth',lineWidth);
      hold on;           
      
      boxEdgeColor = dataColor(1,:);
    else
      fill(dataBox(:,1),dataBox(:,2).*dataScale,dataBodyColor,'EdgeColor','none');
      hold on;           
      
      boxEdgeColor = [1,1,1];
    end
             
    
    plot([    xPosition-0.5*dataBodyWidth;xPosition+0.5*dataBodyWidth],...
         [statsStruct.median; statsStruct.median].*dataScale,...
         '-','Color',boxEdgeColor,'LineWidth',lineWidth);
    hold on;  
  case 1
    colorDifference = norm(dataColor(1,:)-dataBodyColor);
    
    if(colorDifference == 0)
      plot([xPosition],...
       [statsStruct.median].*dataScale,...
       '.','Color',dataColor(1,:),'LineWidth',lineWidth,...
       'MarkerSize',10,'MarkerFaceColor',dataBodyColor);
      hold on;    
      
    else    
      plot([xPosition],...
       [statsStruct.median].*dataScale,...
       'o','Color',dataColor(1,:),'LineWidth',0.25*lineWidth,...
       'MarkerSize',3,'MarkerFaceColor',dataBodyColor);
      hold on;    
    end

  case 2
    
    dataBoxMin=statsStruct.p25;
    dataBoxMax=statsStruct.p75;
    dataBox = [ xPosition+dataBodyWidth*0.5, dataBoxMin;...
               xPosition+dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMin;... 
               xPosition+dataBodyWidth*0.5,  dataBoxMin];

    normWhite = norm([1,1,1]-dataBodyColor);
    boxEdgeColor = dataColor(1,:);
    
    boxColorHalfTone = [1,1,1].*0.5 + dataColor(1,:).*0.5;
    
    fill(dataBox(:,1),dataBox(:,2).*dataScale,boxColorHalfTone,...
        'EdgeColor',dataColor(1,:));
    hold on;           

    %plot([    xPosition-0.5*dataBodyWidth;xPosition+0.5*dataBodyWidth],...
    %     [statsStruct.median; statsStruct.median].*dataScale,...
    %     '-','Color',[1,1,1],'LineWidth',2*lineWidth);
    %hold on;     
    
    plot([    xPosition-0.7*dataBodyWidth;xPosition+0.5*dataBodyWidth],...
         [statsStruct.median; statsStruct.median].*dataScale,...
         '-','Color',boxEdgeColor,'LineWidth',lineWidth);
    hold on;      
    
  case 3
    
    dataBoxMin=statsStruct.p25;
    dataBoxMax=statsStruct.p75;
    dataBox = [ xPosition+dataBodyWidth*0.5, dataBoxMin;...
               xPosition+dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMax;... 
               xPosition-dataBodyWidth*0.5,  dataBoxMin;... 
               xPosition+dataBodyWidth*0.5,  dataBoxMin];

    normWhite = norm([1,1,1]-dataBodyColor);
    boxEdgeColor = dataColor(1,:);
    
    boxColorHalfTone = [1,1,1].*0.5 + dataColor(1,:).*0.5;
    
    fill(dataBox(:,1),dataBox(:,2).*dataScale,[1,1,1],...
        'EdgeColor',dataColor(1,:));
    hold on;           

    plot([    xPosition-0.5*dataBodyWidth;xPosition+0.5*dataBodyWidth],...
         [statsStruct.median; statsStruct.median].*dataScale,...
         '-','Color',dataColor(1,:),'LineWidth',lineWidth);
    hold on;      
    
  otherwise assert(0)
end






