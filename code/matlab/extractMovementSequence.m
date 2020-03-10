function [sitToStandSequence] = ...
    extractMovementSequence(c3DTime,segmentedData, segmentInfo, ...
               thresholdQuietTime,...
               movementSequenceDataFileName,...
               flag_loadSequenceDataFromFile)
                



indexSittingStatic      = 0;
indexSittingDynamic     = 1;

indexCrouchingStable    = 2;
indexCrouchingUnstable  = 3;

indexStandingStable     = 4;
indexStandingUnstable   = 5;  


assert(strcmp(segmentInfo(1).phaseName,'sittingStatic')     ==1);
assert(strcmp(segmentInfo(2).phaseName,'sittingDynamic')     ==1);
 
assert(strcmp(segmentInfo(3).phaseName,'crouchingStable')   ==1);
assert(strcmp(segmentInfo(4).phaseName,'crouchingUnstable') ==1);

assert(strcmp(segmentInfo(5).phaseName,'standingStable')    ==1);
assert(strcmp(segmentInfo(6).phaseName,'standingUnstable')  ==1);



if(flag_loadSequenceDataFromFile == 0)

  %sequenceData = struct('sit2stand',[]);
                  
  sitToStandSequence(1) = struct('phaseTransitions', zeros(1,2).*NaN,... 
                    'phaseTransitionTimes', zeros(1,2).*NaN,...
                    'phaseTransitionIndices',zeros(1,2).*NaN);
  countS2S = 1;
  %Identify quiet-sitting to quiet-standing movement sequences
  indexListStableStanding = find(segmentedData.phaseTransitions(:,2) == indexStandingStable);
      
  
  %1. Ensure that the quiet standing period is longer than the required
  %   threshold
  %2. Go backwards from the standing trial until the first static sitting 
  %   phase is reached.
  %3. If this static sitting phase lasts longer than the required threshold
  %   then copy the sequence over and append it sit-to-stand

  if(isempty(indexListStableStanding)==0)
    for i=1:1:length(indexListStableStanding)
      idx = indexListStableStanding(i);
      t0 = segmentedData.phaseTransitionTimes(idx);
      t1 = 0;
      idxStart=0;
      idxEnd = 0;
      if(idx < length(segmentedData.phaseTransitionTimes))
        t1 = segmentedData.phaseTransitionTimes(idx+1,1);
        idxEnd = segmentedData.phaseTransitionIndices(idx+1,1);
      else
        %At the very end of the sequence
        t1 = c3DTime(end);
        idxEnd = length(c3DTime);
      end
      
      %Make sure the standing phase is sufficiently long
      if( (t1-t0) >= thresholdQuietTime)
        %Go back to the very first transition from static-sitting to 
        %something else
        k=idx;        
        while( k > 1 && segmentedData.phaseTransitions(k,1) ~= indexSittingStatic)          
          k=k-1;
        end
        t1 = segmentedData.phaseTransitionTimes(k,1);
        t0 = 0;
        if(k > 1)
          t0 = segmentedData.phaseTransitionTimes(k-1,1);
          idxStart = segmentedData.phaseTransitionIndices(idx-1,1);
        else
          t0 = c3DTime(1,1);
          idxStart=1;
        end
        
        %Ensure that the length of this static-sitting period exceeds the
        %required threshold
        if((t1-t0) >= thresholdQuietTime)
       
          sitToStandSequence(countS2S).phaseTransitions = ...
            segmentedData.phaseTransitions(k:1:idx,:);
          sitToStandSequence(countS2S).phaseTransitionTimes = ...
            segmentedData.phaseTransitionTimes(k:1:idx,1);
          sitToStandSequence(countS2S).phaseTransitionIndices = ...
            segmentedData.phaseTransitionIndices(k:1:idx,1);
          
          countS2S= countS2S + 1;
        end
        
      end
      
      
    end
  end
                  
  save(movementSequenceDataFileName, 'sitToStandSequence');
else
  load(movementSequenceDataFileName);
end