function [segmentedData, segmentInfo] = segmentDataKMeans(...
                      c3DTime,...
                      c3DMarkers,...
                      c3DFootMarkerNames,...
                      comPosition,...                      
                      c3DGrfChair,...
                      c3DGrfFeet,...
                      fpeDataGroundPlane,...                      
                      standingHeightLowerBound,...
                      segmentedDataFileName,...
                      flag_loadTimingDataFromFile)
%highNormalForce,...                  
%segmentedData = [];  
%segmentInfo = [];
if(flag_loadTimingDataFromFile == 0)

  
  assert(size(c3DMarkers.(c3DFootMarkerNames{1}),1)==size(c3DTime,1));
  assert(size(c3DTime,1) == size(c3DGrfChair.force,1));
  assert(size(c3DGrfChair.force,1)==size(c3DGrfFeet.force,1));
  assert(size(c3DGrfChair.force,1)==size(fpeDataGroundPlane.r0G0,1));
  
  indexSittingStatic      = 0;
  indexSittingDynamic     = 1;

  indexCrouchingStable    = 2;
  indexCrouchingUnstable  = 3;
  
  indexStandingStable     = 4;
  indexStandingUnstable   = 5;  
    
  indexUnclassified = NaN;

  segmentInfo(6) = struct('phaseId',0,...
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

  
  segmentedData = struct('phase', zeros(size(c3DGrfChair.force,1),1),...                      
                      'phaseTransitions', zeros(1,2).*NaN,... 
                      'phaseTransitionTimes', zeros(1,1).*NaN,...
                      'phaseTransitionIndices',zeros(1,1).*NaN);                    

  phasePrev = indexUnclassified;
  
  footMarkers = zeros(length(c3DFootMarkerNames),3);
  
  optionsKMeans = statset('Display','off');
  
  %Solve for 3 clusters on the chair side force plate
  clusterChairForcePlate = kmeans(c3DGrfChair.force(:,3),3,...
                                  'Options',optionsKMeans);  
  
  id1 = find(clusterChairForcePlate ==1);
  id2 = find(clusterChairForcePlate ==2);
  id3 = find(clusterChairForcePlate ==3);
  
  fz1 = mean(c3DGrfChair.force(id1,3));
  fz2 = mean(c3DGrfChair.force(id2,3));
  fz3 = mean(c3DGrfChair.force(id3,3));

  [val,idx] = sort([fz1;fz2;fz3]);

  idFzStanding   = idx(1);
  idFzTransition = idx(2);
  idFzSitting    = idx(3);
  
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
  
  for i=1:1:size(c3DGrfChair.force,1)
    
    comHeight = comPosition(i,3);


    
    for j=1:1:size(footMarkers,1)
      footMarkers(j,:) = c3DMarkers.(c3DFootMarkerNames{j})(i,:);
    end
             
    r0F0 = fpeDataGroundPlane.r0F0(i,:);
    
    %Check to see if r0F0 is within the convex hull of the foot markers.
    indexCH = convhull(footMarkers(:,1:2));
    convexHullCenter = [mean(footMarkers(indexCH,1:2))];
    rCF0 = r0F0(1,1:2) - convexHullCenter;
    
    footMarkersCentered = footMarkers(indexCH,1:2) - convexHullCenter;
    
    distanceToConvexHull = calcDistanceToConvexHull(rCF0, ...
                            footMarkersCentered);
    
    
    if(i==832)
      here=1;
    end
                
    
    phase = indexUnclassified;
    
    if(    isnan(comHeight)==0 && isnan(lfpe(i,1))==0 ...
        && sum(isnan(r0F0))==0 && isnan(distanceToConvexHull) == 0)
      
      if(    (clusterChairForcePlate(i,1)==idFzSitting) ...
          || (clusterChairForcePlate(i,1)==idFzTransition) )   
        
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
      else
      %Off the seat
        if( comHeight <= standingHeightLowerBound)  
          
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
        else
        %Standing
          if(distanceToConvexHull > 0)
          %Unstable
            phase = indexStandingUnstable;
          else
          %Stable
            phase = indexStandingStable;
          end
        end
      end

      assert(isnan(phase)==0);

    end
    
    %Check if there has been a phase transition
    if(   isnan(phasePrev) == 0 ...
       && isnan(phase) == 0 ...
       && phasePrev ~= phase)
     
     if(phasePrev == indexSittingStatic && phase == indexStandingStable)
       here=1;
     end
     
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
        
  end
  save(segmentedDataFileName,'segmentedData','segmentInfo');
else
%  if(flag_verbose)
%    disp('Timing Data: reading in saved mat structures');
%  end  
  load(segmentedDataFileName);  
end
                    