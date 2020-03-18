function figH = plotBalancePortraitAtPointInTime(figH, subPlotVec,...
                      subjectNumber, subjectId,...
                      c3dFootMarkerRightNames,c3dFootMarkerLeftNames,...
                      movementSequence,...
                      scaleTime,scaleDistance,scaleAngle,scaleForce,...
                      c3dTime, ...
                      c3dMarkers, ...                      
                      comPosition,...
                      balancePoint, balanceTangentDir, balanceNormalDir,...
                      c3dGrfFeet,...
                      lineColorCom,...
                      lineColorBalancePoint,lineColorBalanceDir,...
                      lineColorCop,...
                      footPatchColor,...
                      xLabelText,yLabelText,...
                      titleText)


figure(figH);
if(length(subPlotVec) == 3)
  subplot(subPlotVec(1,1),subPlotVec(1,2),subPlotVec(1,3));
end  
if(length(subPlotVec) == 4)
  subplot('Position',subPlotVec);
end


patchColor = footPatchColor;
footColor  = footPatchColor.*0.5 +[1,1,1].*0.5;

c3dFootMarkerNames = {c3dFootMarkerRightNames{:},c3dFootMarkerLeftNames{:}};

footMarkers      = zeros(length(c3dFootMarkerNames),3);
footMarkersLeft  = zeros(length(c3dFootMarkerLeftNames),3);
footMarkersRight = zeros(length(c3dFootMarkerRightNames),3);

xExtents = [Inf,-Inf];
yExtents = [Inf,-Inf];

r1C1V = zeros(length(movementSequence),2).*NaN;
r1P1V = zeros(length(movementSequence),2).*NaN;
r1F1V = zeros(length(movementSequence),2).*NaN;
r1H1V = zeros(length(movementSequence),2).*NaN;
r1T1V = zeros(length(movementSequence),2).*NaN;

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)


    idx1 = movementSequence(z).indexReference;
    
    %%
    %
    % Get the convex hull of the feet and project it onto the 
    % tangent direction of the balance line
    %
    %%
    for j=1:1:size(footMarkers,1)
      footMarkers(j,:) = c3dMarkers.(c3dFootMarkerNames{j})(idx1,:);
      
      if(footMarkers(j,1)==0 && footMarkers(j,2)==0 && footMarkers(j,3)==0)
        footMarkers(j,:) = NaN;
      end
      
    end
    for j=1:1:size(footMarkersLeft,1)
      footMarkersLeft(j,:) = c3dMarkers.(c3dFootMarkerLeftNames{j})(idx1,:);
      if(footMarkersLeft(j,1)==0 && footMarkersLeft(j,2)==0 && footMarkersLeft(j,3)==0)
        footMarkersRight(j,:) = NaN;
      end
      
    end
    for j=1:1:size(footMarkersRight,1)
      footMarkersRight(j,:) = c3dMarkers.(c3dFootMarkerRightNames{j})(idx1,:);
      if(footMarkersRight(j,1)==0 && footMarkersRight(j,2)==0 && footMarkersRight(j,3)==0)
        footMarkersRight(j,:) = NaN;
      end      
    end

    footMarkers      = footMarkers(isnan(footMarkers(:,1))==0,:);
    footMarkersLeft  = footMarkersLeft(isnan(footMarkersLeft(:,1))==0,:);
    footMarkersRight = footMarkersRight(isnan(footMarkersRight(:,1))==0,:);
    
    idxFeetCH       = convhull(footMarkers(:,1:2));
    idxLeftFootCH   = convhull(footMarkersLeft(:,1:2));
    idxRightFootCH  = convhull(footMarkersRight(:,1:2));

    [dHC, nH0, r0H0t] = ...
      calcDistanceToConvexHull(comPosition(idx1,1:2), ...
                               footMarkers(idxFeetCH,1:2));
                             
    %%
    %
    % Form the reference coordinate system
    %
    %%

    eX10 = balanceTangentDir(idx1,1:2)';
    eY10 = balanceNormalDir(idx1,1:2)';    
    rm10 = [eX10,eY10];
    
    
    
    r0H0 = r0H0t';
    r0H1 = rm10'*r0H0;    
    r0C0 = comPosition(idx1,1:2)';
    r0C1 = rm10'*r0C0;
    
    r010 = rm10*[r0H1(1,1);r0C1(2,1)];

    
    %Find the location of the far side of the foot
    r020 = r010 + eX10.*0.325;
    [dT2t,nT0t,r0T0t]=calcDistanceToConvexHull(r020',footMarkers(idxFeetCH,1:2));
    
    r0T0 = r0T0t';
    
    %Get the points related to balance
    r0F0 = balancePoint(idx1,1:2)';    
    r0P0 = c3dGrfFeet.cop(idx1,1:2)';
    
    r1H1 = rm10'*(r0H0-r010);
    r1F1 = rm10'*(r0F0-r010);
    r1C1 = rm10'*(r0C0-r010);
    r1P1 = rm10'*(r0P0-r010);
  
    r1T1 = rm10'*(r0T0-r010);
    
    r1C1V(z,:) = r1C1';
    r1P1V(z,:) = r1P1';
    r1F1V(z,:) = r1F1';
    r1H1V(z,:) = r1H1';
    r1T1V(z,:) = r1T1';    
    

  
    
%     
%     xyOffset = 0.5.*c3dMarkers.(c3dFootMarkerRightNames{end})(idx1,:)...
%              + 0.5.*c3dMarkers.(c3dFootMarkerLeftNames{end})(idx1,:); 
%     
%     
%     for j=1:1:size(footMarkers,1)
%       if(footMarkers(j,1)-xyOffset(1,1) < xExtents(1,1))
%         xExtents(1,1) = footMarkers(j,1)-xyOffset(1,1);
%       end
%       if(footMarkers(j,1)-xyOffset(1,1) > xExtents(1,2))
%         xExtents(1,2) = footMarkers(j,1)-xyOffset(1,1);
%       end
%       if(footMarkers(j,2)-xyOffset(1,2) < yExtents(1,1))
%         yExtents(1,1) = footMarkers(j,2)-xyOffset(1,2);
%       end
%       if(footMarkers(j,2)-xyOffset(1,2) > yExtents(1,2))
%         yExtents(1,2) = footMarkers(j,2)-xyOffset(1,2);
%       end
%     end
% 
% 
% 
%     plot((footMarkers(idxFeetCH,1)-xyOffset(1,1)).*100,...
%          (footMarkers(idxFeetCH,2)-xyOffset(1,2)).*100,...
%          '-','Color',patchColor);
%     hold on;
%     plot((footMarkersLeft(idxLeftFootCH,1)-xyOffset(1,1)).*100,...
%          (footMarkersLeft(idxLeftFootCH,2)-xyOffset(1,2)).*100,...
%          '-','Color',footColor);
%     hold on;
%     plot((footMarkersRight(idxRightFootCH,1)-xyOffset(1,1)).*100,...
%          (footMarkersRight(idxRightFootCH,2)-xyOffset(1,2)).*100,...
%          '-','Color',footColor);
%     hold on;
% 
% 
%     %Balance plane
% 
%     if(flag_sitToStand0SeatOff1 == 0)
%       plot((comPosition(idx1:idx2,1)-xyOffset(1,1)).*100,...
%            (comPosition(idx1:idx2,2)-xyOffset(1,2)).*100,...
%            '-','Color',lineColorCom);
%       hold on;
%     else
%       plot((comPosition(idx1,1)-xyOffset(1,1)).*100,...
%            (comPosition(idx1,2)-xyOffset(1,2)).*100,...
%            '-','Color',lineColorCom);
%       hold on;
%     end
% 
%     minXVal = min(comPosition(idx1:idx2,1))-xyOffset(1,1);
%     maxXVal = max(comPosition(idx1:idx2,1))-xyOffset(1,1);
%     minYVal = min(comPosition(idx1:idx2,2))-xyOffset(1,2);
%     maxYVal = max(comPosition(idx1:idx2,2))-xyOffset(1,2);
%     if(minXVal < xExtents(1,1))
%       xExtents(1,1) = minXVal;
%     end
%     if(maxXVal > xExtents(1,2))
%       xExtents(1,2) = maxXVal;
%     end
%     if(minYVal < yExtents(1,1))
%       yExtents(1,1) = minYVal;
%     end
%     if(maxYVal > yExtents(1,2))
%       yExtents(1,2) = maxYVal;
%     end
% 
%     %Cop-Balance point trajectory
%     xyPlane = comPosition(idx1,:)-xyOffset;
%     %for k = idx1:1:idx2
%     k = idx1;
%     pt0 = comPosition(k,:)-xyOffset;
%     
%     pt1 = pt0 + balanceStepDir(k,:).*0.3; %balancePoint(k,:)-xyOffset;
%     plot([pt0(1,1);pt1(1,1)].*100,...
%          [pt0(1,2);pt1(1,2)].*100,'-','Color',lineColorBalanceDir);
%     hold on;
%     plot([pt0(1,1)].*100,...
%          [pt0(1,2)].*100,'.','Color',lineColorCom,...
%          'MarkerSize',10);
%     hold on;
%     %end
% 
% 
%     %Balance point trajectory
%     %plot((balancePoint(idx1:idx2,1)-xyOffset(1,1)).*100,...
%     %     (balancePoint(idx1:idx2,2)-xyOffset(1,2)).*100,...
%     %     '-','Color',lineColorBalancePoint);
%     %hold on;
% 
%     %Cop trajectory
%     if(flag_sitToStand0SeatOff1 == 0)
%       plot((c3dGrfFeet.cop(idx1:idx2,1)-xyOffset(1,1)).*100,...
%            (c3dGrfFeet.cop(idx1:idx2,2)-xyOffset(1,2)).*100,...
%            '-','Color',lineColorCop);
%       hold on;
%     else
%       plot((c3dGrfFeet.cop(idx1,1)-xyOffset(1,1)).*100,...
%            (c3dGrfFeet.cop(idx1,2)-xyOffset(1,2)).*100,...
%            '-','Color',lineColorCop);
%       hold on;      
%     end
%     
%     %Com end point
%     if(flag_sitToStand0SeatOff1==0)    
%       plot([comPosition(idx2,1)-xyOffset(1,1)].*100,...
%            [comPosition(idx2,2)-xyOffset(1,2)].*100,'o','Color',[0,0,0],...
%            'MarkerSize',5,'MarkerFaceColor',lineColorCom);
%       hold on;
%     else
%       plot([comPosition(idx1,1)-xyOffset(1,1)].*100,...
%            [comPosition(idx1,2)-xyOffset(1,2)].*100,'o','Color',[0,0,0],...
%            'MarkerSize',5,'MarkerFaceColor',lineColorCom);
%       hold on;    
%     end
%     
%     %Balance end point
%     if(flag_sitToStand0SeatOff1==0)
%       plot((balancePoint(idx2,1)-xyOffset(1,1)).*100,...
%            (balancePoint(idx2,2)-xyOffset(1,2)).*100,...
%            'o','Color',lineColorBalancePoint,'MarkerSize',5,...
%            'MarkerFaceColor',[1,1,1]);
%       hold on;
%     else
%       plot((balancePoint(idx1,1)-xyOffset(1,1)).*100,...
%            (balancePoint(idx1,2)-xyOffset(1,2)).*100,...
%            'o','Color',lineColorBalancePoint,'MarkerSize',5,...
%            'MarkerFaceColor',[1,1,1]);
%       hold on;
%     end
%     minXVal = min(balancePoint(idx1:idx2,1))-xyOffset(1,1);
%     maxXVal = max(balancePoint(idx1:idx2,1))-xyOffset(1,1);
%     minYVal = min(balancePoint(idx1:idx2,2))-xyOffset(1,2);
%     maxYVal = max(balancePoint(idx1:idx2,2))-xyOffset(1,2);
%     if(minXVal < xExtents(1,1))
%       xExtents(1,1) = minXVal;
%     end
%     if(maxXVal > xExtents(1,2))
%       xExtents(1,2) = maxXVal;
%     end
%     if(minYVal < yExtents(1,1))
%       yExtents(1,1) = minYVal;
%     end
%     if(maxYVal > yExtents(1,2))
%       yExtents(1,2) = maxYVal;
%     end
%     
%     
% 
%     %Cop end point
%     if(flag_sitToStand0SeatOff1 == 0)
%       plot((c3dGrfFeet.cop(idx2,1)-xyOffset(1,1)).*100,...
%            (c3dGrfFeet.cop(idx2,2)-xyOffset(1,2)).*100,...
%            'o','Color',[0,0,0],'MarkerSize',5,...
%            'MarkerFaceColor',lineColorCop);
%       hold on;
%     else
%       plot((c3dGrfFeet.cop(idx1,1)-xyOffset(1,1)).*100,...
%            (c3dGrfFeet.cop(idx1,2)-xyOffset(1,2)).*100,...
%            'o','Color',[0,0,0],'MarkerSize',5,...
%            'MarkerFaceColor',lineColorCop);
%       hold on;
%     end
%     minXVal = min(c3dGrfFeet.cop(idx1:idx2,1))-xyOffset(1,1);
%     maxXVal = max(c3dGrfFeet.cop(idx1:idx2,1))-xyOffset(1,1);
%     minYVal = min(c3dGrfFeet.cop(idx1:idx2,2))-xyOffset(1,2);
%     maxYVal = max(c3dGrfFeet.cop(idx1:idx2,2))-xyOffset(1,2);
%     if(minXVal < xExtents(1,1))
%       xExtents(1,1) = minXVal;
%     end
%     if(maxXVal > xExtents(1,2))
%       xExtents(1,2) = maxXVal;
%     end
%     if(minYVal < yExtents(1,1))
%       yExtents(1,1) = minYVal;
%     end
%     if(maxYVal > yExtents(1,2))
%       yExtents(1,2) = maxYVal;
%     end    
% 
%     xDelta = diff(xExtents);
%     yDelta = diff(yExtents);
%     delta = max(xDelta,yDelta);
% 
%     if(xExtents(1,1) < 0)
%       here=1;
%     end
%     
%     axis([xExtents(1,1),xExtents(1,1)+delta,...
%           yExtents(1,1),yExtents(1,1)+delta].*100);
%     axis square    
  end
end

%The XY coordinate system of the lab is X forward, Y left
%Remap this to X right, Y Up
r1T1V = [-r1T1V(:,2),r1T1V(:,1)];
r1H1V = [-r1H1V(:,2),r1H1V(:,1)];
r1C1V = [-r1C1V(:,2),r1C1V(:,1)];
r1F1V = [-r1F1V(:,2),r1F1V(:,1)];
r1P1V = [-r1P1V(:,2),r1P1V(:,1)];

deltaDx = 0.025;
dx = subjectNumber*2*deltaDx;

behindHeel = [(dx-deltaDx*0.5),(-.075);...
              (dx+deltaDx*0.5),(-.075);...
              (dx+deltaDx*0.5),(0. );...
              (dx-deltaDx*0.5),(0. );...
              (dx-deltaDx*0.5),(-.075)].*scaleDistance;                
fill(behindHeel(:,1),behindHeel(:,2),footPatchColor,'EdgeColor','none');
hold on;

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexStart))==0)
    

    
%     aheadToe   = [(dx-deltaDx),(r1T1V(z,2));...
%                  (dx+deltaDx),(r1T1V(z,2));...
%                  (dx+deltaDx),(r1T1V(z,2)+0.1);...
%                  (dx-deltaDx),(r1T1V(z,2)+0.1);...
%                  (dx-deltaDx),(r1T1V(z,2))].*scaleDistance;                
%     fill(aheadToe(:,1),aheadToe(:,2),footPatchColor,'EdgeColor','none');
%     hold on;
    plot([dx-deltaDx*0.5;dx+deltaDx*0.5].*scaleDistance,...
         [r1T1V(z,2),r1T1V(z,2)].*scaleDistance,...
         '-','Color', footPatchColor);
    hold on;
    
    plot( [r1C1V(z,1)+dx;r1F1V(z,1)+dx].*scaleDistance, ...
          [   r1C1V(z,2);   r1F1V(z,2)].*scaleDistance,...
         '-','Color',lineColorBalanceDir);
    hold on;
    plot([r1C1V(z,1)+dx].*scaleDistance,...
         [r1C1V(z,2)   ].*scaleDistance,...
         'o','Color',lineColorCom,...
          'MarkerSize',3, 'MarkerFaceColor',[1,1,1]);
    hold on;
    plot([r1F1V(z,1)+dx].*scaleDistance,...
         [r1F1V(z,2)   ].*scaleDistance,...
         '+','Color',lineColorBalancePoint,...
         'MarkerSize',3, 'MarkerFaceColor',lineColorBalancePoint);
    hold on;
    plot([r1P1V(z,1)+dx].*scaleDistance,...
         [r1P1V(z,2)   ].*scaleDistance,...
         '.','Color',lineColorCop,...
         'MarkerSize',3, 'MarkerFaceColor',lineColorCop);
    hold on;  
  end
end
    
  dy = -0.08;  
  
  text( [                dx].*scaleDistance, ...
        [                dy].*scaleDistance, subjectId,...
        'FontSize',6,'Interpreter','latex',...
        'HorizontalAlignment','center');  
  hold on;   
    
  if(isempty(xLabelText)==0)
    xlabel(xLabelText);
  end
  if(isempty(yLabelText)==0)
    ylabel(yLabelText);
  end
  box off;
  title(titleText);  
end