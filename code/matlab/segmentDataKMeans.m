function [segmentedData, segmentInfo] = segmentDataKMeans(...
                      c3DTime,...
                      c3DMarkers,...
                      c3DFootMarkerRightNames,...
                      c3DFootMarkerLeftNames,...                      
                      comPosition,...                      
                      c3DGrfChair,...
                      c3DGrfFeet,...
                      comgp2FootConvexHullDist,...
                      fpeDataGroundPlane,...                      
                      standingHeightLowerBound, ...
                      footContactZMovementTolerance,...
                      numberOfLowFootMarkersForStance,...
                      segmentedDataFileName,...
                      flag_loadTimingDataFromFile,...
                      c3dGrfDataAvailable)
%highNormalForce,...                  
%segmentedData = [];  
%segmentInfo = [];
if(flag_loadTimingDataFromFile == 0)
  
  c3DFootMarkerNames = {c3DFootMarkerRightNames{:},...
                        c3DFootMarkerLeftNames{:}};

  assert(size(c3DMarkers.(c3DFootMarkerNames{1}),1)==size(c3DTime,1));
  
  if(c3dGrfDataAvailable==1)

    assert(size(c3DTime,1) == size(c3DGrfChair.force,1));
    assert(size(c3DGrfChair.force,1)==size(c3DGrfFeet.force,1));
    assert(size(c3DGrfChair.force,1)==size(fpeDataGroundPlane.r0G0,1));

  end
  
  indexSittingStatic      = 0;
  indexSittingDynamic     = 1;

  indexCrouchingStable    = 2;
  indexCrouchingUnstable  = 3;
  
  indexStandingStable     = 4;
  indexStandingUnstable   = 5;  

  indexFootGroundContactBroken = -1;
  
  indexUnclassified = NaN;
  
  segmentInfo(7) = struct('phaseId',0,...
                          'phaseName', {''});
  
  segmentInfo(1).phaseId    = indexSittingStatic;
  segmentInfo(1).phaseName = 'sittingStatic';
  segmentInfo(2).phaseId    = indexSittingDynamic;
  segmentInfo(2).phaseName = 'sittingDynamic';
  segmentInfo(3).phaseId    = indexCrouchingStable;
  segmentInfo(3).phaseName = 'crouchingStable';
  segmentInfo(4).phaseId    = indexCrouchingUnstable;
  segmentInfo(4).phaseName = 'crouchingUnstable';
  segmentInfo(5).phaseId    = indexStandingStable;
  segmentInfo(5).phaseName = 'standingStable';
  segmentInfo(6).phaseId    = indexStandingUnstable;
  segmentInfo(6).phaseName = 'standingUnstable';

  numberOfPoints = 0;
  if(c3dGrfDataAvailable==1)
    numberOfPoints = size(c3DGrfChair.force,1);
  else
    numberOfPoints = size(comPosition,1);
  end
  
  
  segmentedData = struct('phase', zeros(numberOfPoints,1),...                      
                      'phaseTransitions', zeros(1,2).*NaN,... 
                      'phaseTransitionTimes', zeros(1,1).*NaN,...
                      'phaseTransitionIndices',zeros(1,1).*NaN,...
                      'footContactLeft',zeros(numberOfPoints,1),...
                      'footContactRight',zeros(numberOfPoints,1));                    

  phasePrev = indexUnclassified;
  
  footMarkers = zeros(length(c3DFootMarkerNames),3);
  
  optionsKMeans = statset('Display','off');
  
  %Here the letter 'X' is used in the name because these thresholds are
  %identified using different methods depending on what data is available:
  % Ground forces, if available
  % com-gp distance to BOS if available.
  idXStanding   = 0;
  idXTransition = 0;
  idXSitting    = 0;  
  
  clusterX = [];
  
  if(c3dGrfDataAvailable == 1)
    %Solve for 3 clusters on the chair side force plate
    clusterX = kmeans(c3DGrfChair.force(:,3),3,...
                      'Options',optionsKMeans);  

    id1 = find(clusterX ==1);
    id2 = find(clusterX ==2);
    id3 = find(clusterX ==3);

    fz1 = mean(c3DGrfChair.force(id1,3));
    fz2 = mean(c3DGrfChair.force(id2,3));
    fz3 = mean(c3DGrfChair.force(id3,3));

    [val,idx] = sort([fz1;fz2;fz3]);

    idXStanding   = idx(1);
    idXTransition = idx(2);
    idXSitting    = idx(3);
  else
    
    %Negative distances are within the convex hull in the data
    %structure: this sign convention is changed when the plots 
    %are generated.
    clusterX = kmeans(comgp2FootConvexHullDist.distance,3,...
                      'Options',optionsKMeans);

    id1 = find(clusterX ==1);
    id2 = find(clusterX ==2);
    id3 = find(clusterX ==3);
 
    d1 = mean(comgp2FootConvexHullDist.distance(id1,1));
    d2 = mean(comgp2FootConvexHullDist.distance(id2,1));
    d3 = mean(comgp2FootConvexHullDist.distance(id3,1));

    [val,idx] = sort([d1;d2;d3]);

    idXStanding   = idx(1); %Com cluster furthest inside BOS    
    idXTransition = idx(2); %Transition
    idXSitting    = idx(3); %Com cluster farthest outside BOS
            
  end
  %Solve for 3 clusters on com height
  clusterComHeight = kmeans(comPosition(:,3),3,...
                            'Options',optionsKMeans);
  id1 = find(clusterComHeight ==1);
  id2 = find(clusterComHeight ==2);
  id3 = find(clusterComHeight ==3);
  
  h1 = mean(comPosition(id1,3));
  h2 = mean(comPosition(id2,3));
  h3 = mean(comPosition(id3,3));

  [val,idx] = sort([h1;h2;h3]);

  idComZSitting   = idx(1);
  idComZCrouching = idx(2);
  idComZStanding  = idx(3);
  
  %Solve for 2 foot placement stability clusters
  rGF0 = fpeDataGroundPlane.r0F0 - fpeDataGroundPlane.r0G0;
  lfpe = sum(rGF0.^2,2).^0.5;
    
  
  lfpe( isnan(lfpe)==1 ) =  0.;
  
  clusterFpeStepLength = kmeans(lfpe,3, 'Options',optionsKMeans);
  id1 = find(clusterFpeStepLength ==1);
  id2 = find(clusterFpeStepLength ==2);
  id3 = find(clusterFpeStepLength ==3);
  
  lfpe1 = mean(lfpe(id1,1));
  lfpe2 = mean(lfpe(id2,1));
  lfpe3 = mean(lfpe(id3,1));

  [val,idx] = sort([lfpe1;lfpe2;lfpe3]);

  idFpeStatic   = idx(1);
  idFpeTransition= idx(2);
  idFpeDynamic  = idx(3);
  
  c3DFootMarkerRightMinimumHeight =...
    zeros(length(c3DFootMarkerRightNames),1);
  c3DFootMarkerLeftMinimumHeight =...
    zeros(length(c3DFootMarkerLeftNames),1);
  
  for j=1:1:length(c3DFootMarkerRightNames)
    c3DFootMarkerRightMinimumHeight(j,1) = ... 
      min(c3DMarkers.(c3DFootMarkerRightNames{j})(:,3)) ...
      + footContactZMovementTolerance;
  end
  
  for j=1:1:length(c3DFootMarkerLeftNames)
    c3DFootMarkerLeftMinimumHeight(j,1) = ... 
      min(c3DMarkers.(c3DFootMarkerLeftNames{j})(:,3)) ...
      + footContactZMovementTolerance;
  end
  
  footMarkerRightZHeight = zeros(length(c3DFootMarkerRightNames),1);
  footMarkerLeftZHeight  = zeros(length(c3DFootMarkerLeftNames),1);

  nLeftFootBrokeContact = 0;
  nRightFootBrokeContact = 0;
  for i=1:1:numberOfPoints
    
    comHeight = comPosition(i,3);

    
    
    for j=1:1:size(footMarkers,1)
      footMarkers(j,:) = c3DMarkers.(c3DFootMarkerNames{j})(i,:);
    end
    
    %If all of the foot markers break the height threshold the foot
    %is off the ground then flag_bothFeetOnGround will remain zero. 
    %flag_bothFeetOnGround   = 0; 
    flag_leftFootOnGround   = 0;
    flag_rightFootOnGround  = 0;

    for j=1:1:length(c3DFootMarkerRightNames)
      footMarkerRightZHeight(j,1) = ...
        c3DMarkers.(c3DFootMarkerRightNames{j})(i,3) ...
        -c3DFootMarkerRightMinimumHeight(j,1);
    end
    for j=1:1:length(c3DFootMarkerLeftNames)
      footMarkerLeftZHeight(j,1) = ...
        c3DMarkers.(c3DFootMarkerLeftNames{j})(i,3) ...
        -c3DFootMarkerLeftMinimumHeight(j,1);
    end

    if( sum( footMarkerRightZHeight < 0 ) >= numberOfLowFootMarkersForStance)
      flag_rightFootOnGround=1;
    end       
    if( sum( footMarkerLeftZHeight < 0 ) >= numberOfLowFootMarkersForStance)
      flag_leftFootOnGround=1;
    end
    
    %if(flag_leftFootOnGround == 1 && flag_rightFootOnGround == 1)
    %  flag_bothFeetOnGround=1;
    %end
    
    
    
    r0F0 = fpeDataGroundPlane.r0F0(i,:);
    
    %Check to see if r0F0 is within the convex hull of the foot markers.
    indexCH = convhull(footMarkers(:,1:2));
    %convexHullCenter = [mean(footMarkers(indexCH,1:2))];
    %rCF0 = r0F0(1,1:2) - convexHullCenter;
    
    %footMarkersCentered = footMarkers(indexCH,1:2) - convexHullCenter;
        
    distanceToConvexHull = calcDistanceToConvexHull(r0F0(1,1:2), ...
                            footMarkers(indexCH,1:2));
    
    
    if(i==832)
      here=1;
    end
                
    
    phase = indexUnclassified;
    
 
    
    if(    isnan(comHeight)==0 && isnan(lfpe(i,1))==0 ...
        && sum(isnan(r0F0))==0 && isnan(distanceToConvexHull) == 0 )
      
      if(    (clusterX(i,1)==idXSitting) ...
          || (clusterX(i,1)==idXTransition) )   
        
        %c3DGrfChair.force(i,3) >= smallNormalForce ...
        %  && comHeight <= sittingHeightUpperBound)
      %On the seat
        %if(distanceBalance <= staticSittingBalanceThreshold)
        if(    (clusterFpeStepLength(i,1) == idFpeStatic) ...
            || (clusterFpeStepLength(i,1) == idFpeTransition))
        %Static sitting
          phase = indexSittingStatic;
        else
        %Dynamic sitting
          phase = indexSittingDynamic;
        end
      elseif( comHeight <= standingHeightLowerBound)     

          %(clusterComHeight(i,1) == idComZCrouching) ...
          %|| (clusterComHeight(i,1) == idComZSitting))
        %Crouching
          if(distanceToConvexHull > 0)
          %Unstable
            phase = indexCrouchingUnstable;
          else
          %Stable
            phase = indexCrouchingStable;
          end
      elseif(comHeight > standingHeightLowerBound)
        %Standing
          if(distanceToConvexHull > 0)
          %Unstable
            phase = indexStandingUnstable;
          else
          %Stable
            phase = indexStandingStable;
          end
      end

      if(isnan(phase)==1)
        here=1;
      end
      assert(isnan(phase)==0);

    end
    
    
    
%     if(phase == indexFootGroundContactBroken)
%       if(flag_leftFootOnGround == 0)
%         fprintf('    : %i. %s markers low enough %i < %i \n',i, ...
%           'left-foot',...
%            sum( footMarkerLeftZHeight < 0 ),...
%            numberOfLowFootMarkersForStance);        
%       end
%       if(flag_rightFootOnGround ==0)
%         fprintf('    : %i. %s markers low enough %i < %i \n',i, ...
%           'right-foot',...
%            sum( footMarkerRightZHeight < 0 ),...
%            numberOfLowFootMarkersForStance);         
%       end
%     end
    
    %Check if there has been a phase transition
    if(   isnan(phasePrev) == 0  && isnan(phase) == 0 ...
       && phasePrev ~= phase)
     

     
      if(sum(isnan(segmentedData.phaseTransitions))>0)
        segmentedData.phaseTransitions = [phasePrev,phase];
        segmentedData.phaseTransitionTimes = [c3DTime(i,1)];
        segmentedData.phaseTransitionIndices = [i];
      else
        segmentedData.phaseTransitions = ...
            [segmentedData.phaseTransitions;...
             phasePrev,phase];
        segmentedData.phaseTransitionTimes = ...
            [segmentedData.phaseTransitionTimes;...
             c3DTime(i,1)];  
        segmentedData.phaseTransitionIndices = ...
          [segmentedData.phaseTransitionIndices;i];
        here=1;
      end     
    end
    
    
    phasePrev = phase;
    segmentedData.phase(i,1) = phase;
    segmentedData.footContactLeft(i,1) = flag_leftFootOnGround;
    segmentedData.footContactRight(i,1)= flag_rightFootOnGround;
    
    nLeftFootBrokeContact = nLeftFootBrokeContact+(1-flag_leftFootOnGround);
    nRightFootBrokeContact = nRightFootBrokeContact+(1-flag_rightFootOnGround);

    
  end
  
  if( nLeftFootBrokeContact > 0)
    fprintf('      %i Number samples left foot broke contact\n',...
      nLeftFootBrokeContact);    
  end
  if( nRightFootBrokeContact > 0)
    fprintf('      %i Number samples right foot broke contact\n',...
      nRightFootBrokeContact);    
  end
  
  
  save(segmentedDataFileName,'segmentedData','segmentInfo');
else
%  if(flag_verbose)
%    disp('Timing Data: reading in saved mat structures');
%  end  
  load(segmentedDataFileName);  
end
                    