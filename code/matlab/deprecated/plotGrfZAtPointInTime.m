function figH = plotGrfZAtPointInTime(figH, subPlotVec, ...
                   subjectNumber, subjectLabel,...
                   c3dTime, movementSequence,...
                   c3dGrfChair,...
                   lineColor,...
                   figureTitle)


                 
figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end

data = zeros(length(movementSequence),2);

flag_data=0;
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)
    flag_data=1;
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;

    idxSample = idx1;
    
    data(z,1)=subjectNumber;
    data(z,2) = c3dGrfChair.force(idxSample,3);
   
  
  end

end

if(flag_data==1)
  plot( data(:,1),...
        data(:,2),...
        '-','Color',lineColor);    
  hold on;
  plot( data(:,1), data(:,2),...
        'o','Color',lineColor,'MarkerSize',3,'MarkerFaceColor',lineColor);
  hold on;

  text( data(1,1), (max(data(:,2))+2.5), subjectLabel,...
    'FontSize',6,'Interpreter','latex','HorizontalAlignment','center');

  ylabel('Force (N)');
  box off;
  title(figureTitle);  
end