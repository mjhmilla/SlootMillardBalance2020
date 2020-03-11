function indexLastSeatOff = getLastSeatOffIndex(phaseTransitions,...
                            phaseTransitionIndices,...
                            idSittingDynamic,idCrouchingStable,...
                            flag_printAlertForMultipleSeatOffs)

numberOfTransitions =  size(phaseTransitions,1);
indexLastSeatOff = -1;
for indexPhase=numberOfTransitions:-1:1
  phase1 = phaseTransitions(indexPhase,1);
  phase2 = phaseTransitions(indexPhase,2);

  %Mark the index where sit-to-stand first occurs
  if(phase1 <= idSittingDynamic ...
      && phase2 >= idCrouchingStable)
    if(indexLastSeatOff ==-1)
      indexLastSeatOff = phaseTransitionIndices(indexPhase);
    else
      if(flag_printAlertForMultipleSeatOffs==1)
        disp('    :Multiple sit-to-crouch transitions');
      end
      %assert(0,'Error: more than 1 sit-to-stand transition');
    end
  end
end 