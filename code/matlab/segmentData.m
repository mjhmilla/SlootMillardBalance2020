function [segmentedData, segmentInfo] = segmentData(...
                      c3DTime,...
                      c3DMarkers,...
                      c3DFootMarkerNames,...
                      comPosition,...                      
                      c3DGrfChair,...
                      c3DGrfFeet,...
                      fpeDataGroundPlane,...
                      smallNormalForce,...
                      staticSittingBalanceThreshold, ...
                      standingHeightLowerBound,...
                      sittingHeightUpperBound,...
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
  assert(smallNormalForce > 0);% && smallNormalForce < highNormalForce)
  assert(staticSittingBalanceThreshold > 0);
  assert(sittingHeightUpperBound < standingHeightLowerBound);
  
  indexSittingStatic      = 0;
  indexSittingDynamic     = 1;

  indexCrouchingStable    = 2;
  indexCrouchingUnstable  = 3;
  
  indexStandingStable     = 4;
  indexStandingUnstable   = 5;  
    
  indexUnclassified = NaN;
%   
%   [indexSittingStatic,indexSittingDynamic,...
%                                         indexCrouchingStable,indexCrouchingUnstable,...
%                                         indexStandingStable,indexStandingUnstable,...
%                                         indexUnclassified]
%   
% ...
%                               {'sittingStatic',  'sittingDynamic',...
%                                'crouchingStable','crouchingUnstable',...
%                                'standingStable','standingUnstable',...
%                                'unclassified'}
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
                    
  for i=1:1:size(c3DGrfChair.force,1)
    
    comHeight = comPosition(i,3);

    distanceBalance = norm(  fpeDataGroundPlane.r0F0(i,:) ...
                           - fpeDataGroundPlane.r0G0(i,:));
    
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
    
    
    if(i==877)
      here=1;
    end
                
    
    phase = indexUnclassified;
    
    if(    isnan(comHeight)==0 && isnan(distanceBalance)==0 ...
        && sum(isnan(r0F0))==0 && isnan(distanceToConvexHull) == 0)
      
      if(    c3DGrfChair.force(i,3) >= smallNormalForce ...
          && comHeight <= sittingHeightUpperBound)
      %On the seat
        if(distanceBalance <= staticSittingBalanceThreshold)
        %Static sitting
          phase = indexSittingStatic;
        else
        %Dynamic sitting
          phase = indexSittingDynamic;
        end
      else
      %Off the seat
        if(comHeight < standingHeightLowerBound)
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
                    