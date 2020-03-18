function figH = plotAntsOnALog(figH, subPlotVec, ...
                          xPosition, dataMin,dataMean,dataMax, ...
                          dataBoxMin,dataBoxMax, dataEvents, ...
                          eventMarker, eventMarkerFaceColor,...
                          width, lineColor)

                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end



plot([xPosition;xPosition],[dataMin;dataMax],'-','Color',lineColor,...
     'LineWidth',0.5);
hold on;

plot([xPosition],[dataMin],'.','Color',lineColor,'MarkerSize',3);
hold on;

plot([xPosition],[dataMax],'.','Color',lineColor,'MarkerSize',3);
hold on;


if(isempty(dataEvents)==0)
  signOffset = -1;
  for k=1:1:size(dataEvents,2)  
    for j=1:1:size(dataEvents,1)


      plot([xPosition+width*signOffset, xPosition],...
           [dataEvents(j,k),dataEvents(j,k)],'-','Color',lineColor,...
           'LineWidth',0.5);
      hold on;    


      plot([xPosition+width*signOffset],...
           [dataEvents(j,k)],eventMarker{k},'Color',lineColor,...
           'MarkerSize',2,'MarkerFaceColor',eventMarkerFaceColor(k,:),...
           'LineWidth',0.5);
      hold on;    

    end
    signOffset = -1*signOffset;
  end
end


dataBox = [ xPosition+width*0.5, dataBoxMin;...
           xPosition+width*0.5,  dataBoxMax;... 
           xPosition-width*0.5,  dataBoxMax;... 
           xPosition-width*0.5,  dataBoxMin;... 
           xPosition+width*0.5,  dataBoxMin];


fill(dataBox(:,1),dataBox(:,2),lineColor,'EdgeColor','none');
hold on;

plot([xPosition-0.5*width;xPosition+0.5*width],...
     [             dataMean; dataMean            ],...
     '-','Color',[1,1,1],'LineWidth',1.5);
hold on;  