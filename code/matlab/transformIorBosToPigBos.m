function footDataPig = transformIorBosToPigBos(footData,footwareType)

footDataPig=footData;

flag_debug=1;

for indexFootware=1:1:length(footwareType)

  if(flag_debug==1)
    if(exist('figDebug')==1)
      clf(figDebug);
    else
      figDebug=figure;
    end
  end
  
  footware = footwareType{indexFootware};
  
  footLength  = footData.(footware).length;
  footWidth   = footData.(footware).width;
  
  %Put the average marker positions into units of cm
  FM1 = footData.(footware).markersNorm.('FM1')(1,1:2).*[footWidth,footLength];
  FCC = footData.(footware).markersNorm.('FCC')(1,1:2).*[footWidth,footLength];
  FAL = footData.(footware).markersNorm.('FAL')(1,1:2).*[footWidth,footLength];
  
  %Evaluate the location and orientation of the new foot frame 
  %root using the PiG foot model: FCC, FM1, and FAL  

  %PiG foot y axis
  ey = FM1-FCC;
  ey = ey./norm(ey);
  
  rHA = FAL-FCC;
  eHA = rHA./norm(rHA);
  
  

  %PiG foot x axis
  ex = -(eHA - (eHA*ey').*ey);
  ex = ex./norm(ex);
  EP0 = [ex', ey'];


  %PiG foot origin
  r0P0 = FCC + (rHA*ey').*ey;

  if(flag_debug==1)
    
    axisLength = 0.25;
    subplot(1,2,1);
    plot(0,0,'ok');
    hold on;
    plot([0,1].*axisLength,...
         [0,0].*axisLength,'-k');
    hold on
    plot([0,0].*axisLength,...
         [0,1].*axisLength,'-k');
    hold on
    grid on;
    
    
    subplot(1,2,1);    
    xAxisPts = [r0P0;(r0P0 + ex.*(footWidth*0.5))];
    yAxisPts = [r0P0;(r0P0 + ey.*(footWidth*0.5))];
    plot(r0P0(1,1).*(1/footLength),...
         r0P0(1,2).*(1/footLength),'ob');
    hold on;
    plot(xAxisPts(:,1).*(1/footLength),...
         xAxisPts(:,2).*(1/footLength),'b');
    hold on;
    plot(yAxisPts(:,1).*(1/footLength),...
         yAxisPts(:,2).*(1/footLength),'b');
    hold on;

    subplot(1,2,2);
    plot(0,0,'ok');
    hold on;
    plot([0,1].*axisLength,...
         [0,0].*axisLength,'-b');
    hold on
    plot([0,0].*axisLength,...
         [0,1].*axisLength,'-b');
    hold on
    grid on;
    
    
    grid on;
    
  end
  
  assert( abs(ex*ey') < 1e-6 );

  %Transform all of the fields in footData s.t. all quantities
  %look like they came from the PiG

  %'convhullNorm',[],
  convhull = zeros(size(footData.(footware).convhullNorm));

  r0X0 = footData.(footware).convhullNorm;
  r0X0(:,1) = r0X0(:,1).*footWidth;
  r0X0(:,2) = r0X0(:,2).*footLength;
  convhull = transformXY(r0X0,r0P0,EP0');

  footDataPig.(footware).convhullNorm = convhull.*(1/footLength);

  if(flag_debug==1)
    subplot(1,2,1);        
    plot( footData.(footware).convhullNorm(:,1).*(footWidth/footLength),...
          footData.(footware).convhullNorm(:,2),...
          '-k');
    hold on;
    
    subplot(1,2,2);    
    plot( footDataPig.(footware).convhullNorm(:,1),...
          footDataPig.(footware).convhullNorm(:,2),...
          '-b');
    hold on;
    here=1;
    
  end
  
  %'markersNorm',[],...
  markers = footData.(footware).markersNorm;
  markersLabels = fields(markers);
  for indexMarker=1:1:length(markersLabels)
    markerName = markersLabels{indexMarker};
    
    r0X0 = footData.(footware).markersNorm.(markerName);
    r0X0(:,1) = r0X0(:,1).*footWidth;
    r0X0(:,2) = r0X0(:,2).*footLength;
    r0X0(:,3) = r0X0(:,3).*footWidth;
    r0X0(:,4) = r0X0(:,4).*footLength;

    markers.(markerName)(1,1:2) = transformXY(r0X0(1,1:2),r0P0,EP0');
    markers.(markerName)(1,3:4) = transformXY(r0X0(1,3:4),[0,0],EP0');
    markers.(markerName)(1,:) = markers.(markerName)(1,:) .* (1/footLength);
    here=1;
  end

  footDataPig.(footware).markersNorm = markers;

  %%
  %Plot the markers
  %%
  if(flag_debug==1)
      
    %IOR model
    subplot(1,2,1);
    markerSummaryFields = fields(footData.(footware).markersNorm);  
    for z=1:1:length(markerSummaryFields)
      markerName =markerSummaryFields{z};

      x0    = footData.(footware).markersNorm.(markerName)(:,1);
      y0    = footData.(footware).markersNorm.(markerName)(:,2);    
      xStd  = footData.(footware).markersNorm.(markerName)(:,3);
      yStd  = footData.(footware).markersNorm.(markerName)(:,4); 

      mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);

      fill( mkrEllipse(:,1).*(footWidth/footLength),...
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
    
    %PiG model
    subplot(1,2,2);    
    markerSummaryFields = fields(footDataPig.(footware).markersNorm);  
    for z=1:1:length(markerSummaryFields)
      markerName =markerSummaryFields{z};

      x0    = footDataPig.(footware).markersNorm.(markerName)(:,1);
      y0    = footDataPig.(footware).markersNorm.(markerName)(:,2);    
      xStd  = footDataPig.(footware).markersNorm.(markerName)(:,3);
      yStd  = footDataPig.(footware).markersNorm.(markerName)(:,4); 

      mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);

      fill( mkrEllipse(:,1),...
            mkrEllipse(:,2),...
            [0,0,1].*0.75,'EdgeColor',[0,0,1].*0.5);
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
    
  end
  %'x',[],
  %'y',[],
  
  %processes the left and right foot
  for indexFoot=1:1:2
  
    r0X0x = footData.(footware).x(:,indexFoot) .* footWidth;
    r0X0y = footData.(footware).y(:,indexFoot) .* footLength;

    rPXP = transformXY([r0X0x,r0X0y],r0P0,EP0');
    x = rPXP(:,1);
    y = rPXP(:,2);
      
    footDataPig.(footware).x(:,indexFoot) = x .* (1/footLength);
    footDataPig.(footware).y(:,indexFoot) = y .* (1/footLength);

    r0X0x = footData.(footware).xStd(:,indexFoot) .* footWidth;
    r0X0y = footData.(footware).yStd(:,indexFoot) .* footLength;
    rPXP = transformXY([r0X0x,r0X0y],[0,0],EP0');
    xStd = rPXP(:,1);
    yStd = rPXP(:,2);

    footDataPig.(footware).xStd(:,indexFoot) = xStd .* (1/footLength);
    footDataPig.(footware).yStd(:,indexFoot) = yStd .* (1/footLength);

    xyFields = {'heelXY','toeXY','latXY','medXY',...
                'heelXYStd','toeXYStd','latXYStd','medXYStd'};

    for indexFields = 1:1:length(xyFields)
      r = r0P0;
      if(contains( xyFields{indexFields},'Std')==1)
        r = [0,0];
      end

      r0X0 = footData.(footware).(xyFields{indexFields})(:,indexFoot)';
      r0X0 = r0X0 .* [footWidth,footLength];
      
      footDataPig.(footware).(xyFields{indexFields})(:,indexFoot) = ...
        transformXY( r0X0,...
                     r,...
                     EP0')';
      footDataPig.(footware).(xyFields{indexFields})(:,indexFoot) = ...
        footDataPig.(footware).(xyFields{indexFields})(:,indexFoot).*(1/footLength);
    end    
  end

  %'length',[],
  %'width',[],
  %'lengthStd',[],
  %'widthStd',[]

end