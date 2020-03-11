function [sitToStandQuietSequence, sitToStandDynamicSequence] = ...
    extractMovementSequence(c3DTime,segmentedData, segmentInfo, ...
               thresholdQuietTime,...
               movementSequenceDataFileName,...
               flag_loadSequenceDataFromFile, flag_verbose)
                



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
                  
  sitToStandQuietSequence(1) = struct('phaseTransitions', zeros(1,2).*NaN,... 
                    'phaseTransitionTimes', zeros(1,2).*NaN,...
                    'phaseTransitionIndices',zeros(1,2).*NaN);
                 
  sitToStandDynamicSequence(1) = struct('phaseTransitions', zeros(1,2).*NaN,... 
                    'phaseTransitionTimes', zeros(1,2).*NaN,...
                    'phaseTransitionIndices',zeros(1,2).*NaN);
                  
  countS2Sq = 1;
  countS2Sd = 1;
  %Identify quiet-sitting to quiet-standing movement sequences
  indexListStableStanding = find(segmentedData.phaseTransitions(:,2) == indexStandingStable);
   
  if(flag_verbose==1)
    fprintf('\t%i : Standing Stable Events\n',...
              length(indexListStableStanding));
  end
  
  %1. Ensure that the quiet standing period is longer than the required
  %   threshold
  %2. Go backwards from the standing trial until the first static sitting 
  %   phase is reached.
  %3. If this static sitting phase lasts longer than the required threshold
  %   then copy the sequence over and append it sit-to-stand

  if(isempty(indexListStableStanding)==0)
    for i=1:1:length(indexListStableStanding)
      idx = indexListStableStanding(i,1);
      
      prevPhase = segmentedData.phaseTransitions(idx,1);
      flag_validTransition = 1;
      if(prevPhase == indexSittingDynamic ...
           || prevPhase == indexSittingStatic)
        flag_validTransition = 0; 
        %Happens when the subjects feet touch the chair force plate.
      end
      
      %Make sure this sequence goes back to a static sitting sequence
      %without switching between standing sequences ...
      if(i > 1)
        flag_sittingFound = 0;
        for z = idx:-1:indexListStableStanding(i-1,1) 
          if(segmentedData.phaseTransitions(z,1) == indexSittingStatic ...
              || segmentedData.phaseTransitions(z,1) == indexSittingDynamic)
            flag_sittingFound=1;
          end
        end
        
        if(flag_sittingFound==0)
          flag_validTransition = -1;
        end
        
      end
      
      
      if(flag_validTransition ~= 1)
        if(flag_verbose==1)
          if(flag_validTransition == 0)
            fprintf('\t#%i : false sequence (feet touching chair force plate)\n',i);
          end
          if(flag_validTransition == -1)
            fprintf('\t#%i : incomplete sequence (no sitting between standing events)\n',i);
          end
          
        end
      else
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
        timeQuietStanding = t1-t0;
        %Make sure the standing phase is sufficiently long
        if( (timeQuietStanding) >= thresholdQuietTime)
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
          timeQuietSitting = (t1-t0);
          %Ensure that the length of this static-sitting period exceeds the
          %required threshold
          if((timeQuietSitting) >= thresholdQuietTime)

            sitToStandQuietSequence(countS2Sq).phaseTransitions = ...
              segmentedData.phaseTransitions(k:1:idx,:);
            sitToStandQuietSequence(countS2Sq).phaseTransitionTimes = ...
              segmentedData.phaseTransitionTimes(k:1:idx,1);
            sitToStandQuietSequence(countS2Sq).phaseTransitionIndices = ...
              segmentedData.phaseTransitionIndices(k:1:idx,1);
            countS2Sq= countS2Sq + 1;

            if(flag_verbose==1)
              fprintf('\t#%i : quiet-sit to quiet-stand\n',i);
            end
          else
            sitToStandDynamicSequence(countS2Sd).phaseTransitions = ...
              segmentedData.phaseTransitions(k:1:idx,:);
            sitToStandDynamicSequence(countS2Sd).phaseTransitionTimes = ...
              segmentedData.phaseTransitionTimes(k:1:idx,1);
            sitToStandDynamicSequence(countS2Sd).phaseTransitionIndices = ...
              segmentedData.phaseTransitionIndices(k:1:idx,1);              
            countS2Sd = countS2Sd+1;

            if(flag_verbose==1)
              fprintf(['\t#%i : dynamic-sit to quiet-stand ',...
                      '(quiet sitting time: %f)\n'],i, timeQuietSitting);
            end          
          end

        else
          if(flag_verbose==1)
            fprintf('\t#%i Rejected: Quiet standing too short (%f)\n',i,...
                    (timeQuietStanding));
          end
        end      
      end
    end
  end
      
  countS2Sq = countS2Sq-1;
  countS2Sd = countS2Sd-1;
  
  fprintf('\t%i: Successful quiet-sit-to-quiet-stand\n',countS2Sq);
  fprintf('\t%i: Successful dynamic-sit-to-quiet-stand\n',countS2Sd);
  save(movementSequenceDataFileName,...
        'sitToStandQuietSequence','sitToStandDynamicSequence');
else
  load(movementSequenceDataFileName);
end