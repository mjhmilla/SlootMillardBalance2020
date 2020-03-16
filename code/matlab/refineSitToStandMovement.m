function [indexQuietSitting, indexSeatOff, indexStanding] = ...
    refineSitToStandMovement(phaseTransitions,...
                            phaseTransitionIndices,...
                            indexSittingDynamic,...
                            indexCrouchingStable,...
                            comVelocity,...
                            toleranceQuietSittingComVelocity,...
                            c3dGrfChair,...
                            toleranceSeatOffChairForce,...
                            comPosition,...
                            toleranceStandingHeight,...
                            flag_printAlertForMultipleSeatOffs)


indexQuietSittig = NaN;
indexSetOff = NaN;
indexStanding = NaN;

numberOfTransitions =  size(phaseTransitions,1);
indexLastSeatOff = -1;
index = 0;
for indexPhase=numberOfTransitions:-1:1
  phase1 = phaseTransitions(indexPhase,1);
  phase2 = phaseTransitions(indexPhase,2);

  %Mark the index where sit-to-stand first occurs
  if(phase1 <= indexSittingDynamic ...
      && phase2 >= indexCrouchingStable)
    if(indexLastSeatOff ==-1)
      indexLastSeatOff = phaseTransitionIndices(indexPhase);
      index = indexPhase;
    else
      if(flag_printAlertForMultipleSeatOffs==1)
        disp('    :Multiple sit-to-crouch transitions');
      end
      %assert(0,'Error: more than 1 sit-to-stand transition');
    end
  end
end 

%Polish the estimate of seatoff
indexStand   = phaseTransitionIndices(end);
grfZ = c3dGrfChair.force(indexLastSeatOff:1:indexStand,3);
fzMean = mean(grfZ);
grfZ( grfZ < 0) = fzMean;



fzMin  = min(grfZ);
fzMean = mean(grfZ);

fzDelta = fzMin + 0.125*(fzMean - fzMin);

flag_minFound = 0;
i = 1;
while i < (indexStand-indexLastSeatOff) && flag_minFound == 0
  if(isnan(grfZ(i))==0)
    if(grfZ(i) <= (fzMin + fzDelta))
      flag_minFound = 1;
    end
  end
  i=i+1;
end
assert(flag_minFound ==1);

i = i-2;

indexLastSeatOff = indexLastSeatOff+i;

