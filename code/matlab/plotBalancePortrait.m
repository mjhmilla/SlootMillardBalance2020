function figH = plotBalancePortrait(figH, subPlotVec,...
                      subjectNumber, subjectId,...
                      c3dTime, c3dMarkers, ...
                      c3dFootMarkerRightNames,c3dFootMarkerLeftNames,...
                      movementSequence, segmentedData, ...
                      c3dGrfChair, c3dGrfFeet,...
                      comPosition, balanceStepDir, balancePoint,...
                      idSittingDynamic, idCrouchingStable,...
                      lineColorBalanceDir,lineColorBalancePoint,...
                      lineColorCop, lineColorCom,...
                      footPatchColor,...
                      figTitle,...
                      flag_zeroAtSeatOff,...
                      allowedFootMarkerMovement)

assert(flag_zeroAtSeatOff==1); %For now: this is the only option                 
if(flag_zeroAtSeatOff==1)
  assert(idSittingDynamic==1);
  assert(idCrouchingStable==2);
end                      

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

for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).phaseTransitions))==0)

    idx0 = movementSequence(z).phaseTransitionIndices(1);
    idx1 = movementSequence(z).phaseTransitionIndices(2);

    if(flag_zeroAtSeatOff==1)
      flag_alertMultipleSeatOffs=1;
      idx1 = getLastSeatOffIndex( movementSequence(z).phaseTransitions,...
                          movementSequence(z).phaseTransitionIndices,...
                          idSittingDynamic,idCrouchingStable,...
                          flag_alertMultipleSeatOffs,c3dGrfChair);        
    end

    %Refine the standing position
    idx2 = movementSequence(z).phaseTransitionIndices(end);
    idx3 = idx2;
    phase2 = movementSequence(z).phaseTransitions(end,2);
    phase3 = phase2;
    while    idx3 < size(segmentedData.phase,1) ...
          && segmentedData.phase(idx3,1) == phase2
      idx3 = idx3+1;
    end
    
    maxComHeight = max(comPosition(idx2:idx3,3));
    idx2ComHeight = comPosition(idx2,3);
    targetComHeight = idx2ComHeight + (maxComHeight-idx2ComHeight)*0.95;
    
    while idx2 < idx3 && comPosition(idx2,3) < targetComHeight
      idx2 = idx2+1;
    end
    
    
   
    %Plot the foot outline, and the convex hull outline
    for j=1:1:size(footMarkers,1)
      footMarkers(j,:) = c3dMarkers.(c3dFootMarkerNames{j})(idx1,:);
      
      if(footMarkers(j,1)==0 && footMarkers(j,2)==0 && footMarkers(j,3)==0)
        footMarkers(j,:) = NaN;
      end
      
      maxMarkerMovement = ...
        max(norm(diff(c3dMarkers.(c3dFootMarkerNames{j})(idx1:idx2,:))));
      if(maxMarkerMovement > allowedFootMarkerMovement)
        disp(['    :',c3dFootMarkerNames{j},' moved by ', ...
          num2str(maxMarkerMovement), ' > ', num2str(allowedFootMarkerMovement) ]);
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

    footMarkers = footMarkers(isnan(footMarkers(:,1))==0,:);
    footMarkersLeft = footMarkersLeft(isnan(footMarkersLeft(:,1))==0,:);
    footMarkersRight = footMarkersRight(isnan(footMarkersRight(:,1))==0,:);
    
    idxFeetCH       = convhull(footMarkers(:,1:2));
    idxLeftFootCH   = convhull(footMarkersLeft(:,1:2));
    idxRightFootCH = convhull(footMarkersRight(:,1:2));


    xyOffset = 0.5.*c3dMarkers.(c3dFootMarkerRightNames{end})(idx1,:)...
             + 0.5.*c3dMarkers.(c3dFootMarkerLeftNames{end})(idx1,:); 
    
    
    for j=1:1:size(footMarkers,1)
      if(footMarkers(j,1)-xyOffset(1,1) < xExtents(1,1))
        xExtents(1,1) = footMarkers(j,1)-xyOffset(1,1);
      end
      if(footMarkers(j,1)-xyOffset(1,1) > xExtents(1,2))
        xExtents(1,2) = footMarkers(j,1)-xyOffset(1,1);
      end
      if(footMarkers(j,2)-xyOffset(1,2) < yExtents(1,1))
        yExtents(1,1) = footMarkers(j,2)-xyOffset(1,2);
      end
      if(footMarkers(j,2)-xyOffset(1,2) > yExtents(1,2))
        yExtents(1,2) = footMarkers(j,2)-xyOffset(1,2);
      end
    end

%     fill(footMarkers(idxFeetCH,1)-xyOffset(1,1),...
%          footMarkers(idxFeetCH,2)-xyOffset(1,2),...
%          patchColor, 'EdgeColor','none');
%     hold on;
%     fill(footMarkersLeft(:,1)-xyOffset(1,1),...
%          footMarkersLeft(:,2)-xyOffset(1,2),...
%          footColor, 'EdgeColor','none');
%     hold on;
%     fill(footMarkersRight(:,1)-xyOffset(1,1),...
%          footMarkersRight(:,2)-xyOffset(1,2),...
%          footColor, 'EdgeColor','none');
%     hold on;

    plot((footMarkers(idxFeetCH,1)-xyOffset(1,1)).*100,...
         (footMarkers(idxFeetCH,2)-xyOffset(1,2)).*100,...
         '-','Color',patchColor);
    hold on;
    plot((footMarkersLeft(idxLeftFootCH,1)-xyOffset(1,1)).*100,...
         (footMarkersLeft(idxLeftFootCH,2)-xyOffset(1,2)).*100,...
         '-','Color',footColor);
    hold on;
    plot((footMarkersRight(idxRightFootCH,1)-xyOffset(1,1)).*100,...
         (footMarkersRight(idxRightFootCH,2)-xyOffset(1,2)).*100,...
         '-','Color',footColor);
    hold on;


    %Balance plane

    plot((comPosition(idx1:idx2,1)-xyOffset(1,1)).*100,...
         (comPosition(idx1:idx2,2)-xyOffset(1,2)).*100,...
         '-','Color',lineColorCom);
    hold on;


    minXVal = min(comPosition(idx1:idx2,1))-xyOffset(1,1);
    maxXVal = max(comPosition(idx1:idx2,1))-xyOffset(1,1);
    minYVal = min(comPosition(idx1:idx2,2))-xyOffset(1,2);
    maxYVal = max(comPosition(idx1:idx2,2))-xyOffset(1,2);
    if(minXVal < xExtents(1,1))
      xExtents(1,1) = minXVal;
    end
    if(maxXVal > xExtents(1,2))
      xExtents(1,2) = maxXVal;
    end
    if(minYVal < yExtents(1,1))
      yExtents(1,1) = minYVal;
    end
    if(maxYVal > yExtents(1,2))
      yExtents(1,2) = maxYVal;
    end

    %Cop-Balance point trajectory
    xyPlane = comPosition(idx1,:)-xyOffset;
    %for k = idx1:1:idx2
    k = idx1;
      pt0 = comPosition(k,:)-xyOffset;
      
      pt1 = pt0 + balanceStepDir(k,:).*0.3; %balancePoint(k,:)-xyOffset;
      plot([pt0(1,1);pt1(1,1)].*100,...
           [pt0(1,2);pt1(1,2)].*100,'-','Color',lineColorBalanceDir);
      hold on;
      plot([pt0(1,1)].*100,...
           [pt0(1,2)].*100,'.','Color',lineColorCom,...
           'MarkerSize',10);
      hold on;
    %end


    %Balance point trajectory
    %plot((balancePoint(idx1:idx2,1)-xyOffset(1,1)).*100,...
    %     (balancePoint(idx1:idx2,2)-xyOffset(1,2)).*100,...
    %     '-','Color',lineColorBalancePoint);
    %hold on;

    %Cop trajectory
    plot((c3dGrfFeet.cop(idx1:idx2,1)-xyOffset(1,1)).*100,...
         (c3dGrfFeet.cop(idx1:idx2,2)-xyOffset(1,2)).*100,...
         '-','Color',lineColorCop);
    hold on;

    
    %Com end point
    plot([comPosition(idx2,1)-xyOffset(1,1)].*100,...
         [comPosition(idx2,2)-xyOffset(1,2)].*100,'o','Color',[0,0,0],...
         'MarkerSize',5,'MarkerFaceColor',lineColorCom);
    hold on;
    
    %Balance end point
    hold on;
    plot((balancePoint(idx2,1)-xyOffset(1,1)).*100,...
         (balancePoint(idx2,2)-xyOffset(1,2)).*100,...
         'o','Color',lineColorBalancePoint,'MarkerSize',5,...
         'MarkerFaceColor',[1,1,1]);
    hold on;

    minXVal = min(balancePoint(idx1:idx2,1))-xyOffset(1,1);
    maxXVal = max(balancePoint(idx1:idx2,1))-xyOffset(1,1);
    minYVal = min(balancePoint(idx1:idx2,2))-xyOffset(1,2);
    maxYVal = max(balancePoint(idx1:idx2,2))-xyOffset(1,2);
    if(minXVal < xExtents(1,1))
      xExtents(1,1) = minXVal;
    end
    if(maxXVal > xExtents(1,2))
      xExtents(1,2) = maxXVal;
    end
    if(minYVal < yExtents(1,1))
      yExtents(1,1) = minYVal;
    end
    if(maxYVal > yExtents(1,2))
      yExtents(1,2) = maxYVal;
    end
    
    

    %Cop end point
    plot((c3dGrfFeet.cop(idx2,1)-xyOffset(1,1)).*100,...
         (c3dGrfFeet.cop(idx2,2)-xyOffset(1,2)).*100,...
         'o','Color',[0,0,0],'MarkerSize',5,...
         'MarkerFaceColor',lineColorCop);
    hold on;

    minXVal = min(c3dGrfFeet.cop(idx1:idx2,1))-xyOffset(1,1);
    maxXVal = max(c3dGrfFeet.cop(idx1:idx2,1))-xyOffset(1,1);
    minYVal = min(c3dGrfFeet.cop(idx1:idx2,2))-xyOffset(1,2);
    maxYVal = max(c3dGrfFeet.cop(idx1:idx2,2))-xyOffset(1,2);
    if(minXVal < xExtents(1,1))
      xExtents(1,1) = minXVal;
    end
    if(maxXVal > xExtents(1,2))
      xExtents(1,2) = maxXVal;
    end
    if(minYVal < yExtents(1,1))
      yExtents(1,1) = minYVal;
    end
    if(maxYVal > yExtents(1,2))
      yExtents(1,2) = maxYVal;
    end    

    xDelta = diff(xExtents);
    yDelta = diff(yExtents);
    delta = max(xDelta,yDelta);

    if(xExtents(1,1) < 0)
      here=1;
    end
    
    axis([xExtents(1,1),xExtents(1,1)+delta,...
          yExtents(1,1),yExtents(1,1)+delta].*100);
    axis square    
  end
  box off;

  xlabel('Distance (cm)');
  ylabel('Distance (cm)');
  title(figTitle);
end