function figH =plotFootBosModel(footData,figH,subPlotPanel,subPlotPanelIndex,...
                                normFootAxis,normXTicks,normYTicks, subplotIndexOffset)

figure(figH);                              
footwareType = fields(footData);

for indexFootware=1:1:length(footwareType)
  
  [row,col] = find(subPlotPanelIndex==(indexFootware+subplotIndexOffset));          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);  
  
  footware = footwareType{indexFootware};
  markerSummaryFields = fields(footData.(footware).markersNorm);
  
  %%
  %Plot the markers
  %%
  for z=1:1:length(markerSummaryFields)
    markerName =markerSummaryFields{z};
    
    x0    = footData.(footware).markersNorm.(markerName)(:,1);
    y0    = footData.(footware).markersNorm.(markerName)(:,2);    
    xStd  = footData.(footware).markersNorm.(markerName)(:,3);
    yStd  = footData.(footware).markersNorm.(markerName)(:,4); 
        
    mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);
                
    fill( mkrEllipse(:,1),...
          mkrEllipse(:,2),...
          [1,1,1].*0.75,'EdgeColor',[1,1,1].*0.5);
    hold on;
    
    vAlign = 'bottom';
    hAlign = 'center';
    dy = 0.05;   
    if(y0<0)
      dy = dy*-1.;
      vAlign = 'top';
    end
    
    text(x0,y0+dy,markerSummaryFields{z},...
          'VerticalAlignment',vAlign,...
          'HorizontalAlignment',hAlign);
    hold on;
    
  end  
  
  %%
  %Plot the convex hull
  %%
  plot(footData.(footware).convhullNorm(:,1),...
       footData.(footware).convhullNorm(:,2),...
       '-','Color',[0,0,0],'LineWidth',2);
  hold on;
    
  %%
  %Plot the decoration
  %%
  pointMetrics = {'heel','toe','med','lat'};
  direction = [2;2;1;1];
  hAlign='';
  vAlign='';
  for w=1:1:length(pointMetrics)
  %Plot the standard deviations at the heel, toe, med, la
    xA = mean(footData.(footware).([pointMetrics{w},'XY'])(1,:));
    xB = mean(footData.(footware).([pointMetrics{w},'XY'])(1,:));  

    yA = mean(footData.(footware).([pointMetrics{w},'XY'])(2,:));
    yB = mean(footData.(footware).([pointMetrics{w},'XY'])(2,:));      
    
    delta = 0;
    xLbl = 0;
    yLbl = 0;
    dirLabel='';
    if(direction(w,1)==2)
      delta = mean(footData.(footware).([pointMetrics{w},'XYStd'])(2,:));
      yA = yA -delta;
      yB = yB +delta;
      
      xLbl = xA;
      if(w ==1)
        yLbl = yA;
        vAlign = 'top';
      else
        yLbl = yB;
        vAlign = 'bottom';
      end
      hAlign = 'center';
      dirLabel = 'Y';
    else
      delta = mean(footData.(footware).([pointMetrics{w},'XYStd'])(1,:));
      xA = xA -delta;
      xB = xB +delta;
      
      yLbl = yA;
      
      vAlign = 'bottom';
      if(w ==3)
        hAlign = 'right';
        xLbl = xB;
      else
        hAlign = 'left';        
        xLbl = xA;
      end
      dirLabel = 'X';
    end
    plot( [xA;xB], [yA;yB],'-','Color',[0,0,0],'LineWidth',1); 
    hold on;
    text( xLbl, yLbl, sprintf('%s%1.3f%s',['$\sigma_',dirLabel,'='],delta,'$'),...
         'HorizontalAlignment',hAlign,...
         'VerticalAlignment',vAlign);
    hold on;
  end
  
  
  plot(0,0,'ok','MarkerFaceColor','k');
  hold on;
  plot([0;0.2],[0;0],'-k','MarkerFaceColor','k');
  hold on;
  plot([0;0],[0;0.2],'-k','MarkerFaceColor','k');
  hold on;
  
  xMinCH = min(footData.(footwareType{indexFootware}).convhullNorm(:,1));
  xMaxCH = max(footData.(footwareType{indexFootware}).convhullNorm(:,1));
  xZeroCH = 0;
  xTicksCH = sort([xMinCH,xZeroCH,xMaxCH]);
  yMinCH = min(footData.(footwareType{indexFootware}).convhullNorm(:,2));
  yMaxCH = max(footData.(footwareType{indexFootware}).convhullNorm(:,2));
  yZeroCH = 0;
  yTicksCH = sort([yMinCH,yZeroCH,yMaxCH]);
  
  xlabel('Norm. Width');
  ylabel('Norm. Length');
  title(['Functional BOS : ',footwareType{indexFootware}]);
    
  box off;
  axis(normFootAxis);
  xticks(round(xTicksCH,2));
  yticks(round(yTicksCH,2));
  grid on;

  
  
  
  
end