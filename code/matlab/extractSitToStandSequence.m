function [sitToStandSequence, sitToStandSequenceInfo] = ...
    extractSitToStandSequence(c3dTime,...
                segmentedData, segmentInfo, ...
                comPosition,...
                comVelocity,...
                c3dGrfChair,...
                fpeStepLength,...
                quietDwellTime, ...
                startS2SComVelocityTolerance,...
                seatOffChairForceZTolerance,...
                endS2SComHeightTolerance,...
                endS2SComVelocityTolerance,...
                flag_rejectTrialsFootGroundContactBroken,...
                flag_verbose,...
                c3dGrfDataAvailable)
                
                %movementSequenceDataFileName,...
                %flag_writeCsvSequenceFile,...


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

figVerbose = [];
flag_verbosePlots = 0;
if(flag_verbose==1 && flag_verbosePlots)
  figVerbose = figure;
end


comSpeed = sum(comVelocity.^2,2).^0.5;
%angSpeed = sum(wholeBodyAngularVelocity.^2,2).^0.5;
%sequenceData = struct('sit2stand',[]);

sitToStandSequence(1) = ...
  struct('indexStart'                    , zeros(1,1).*NaN,... 
         'valueSwitchConditionStart'     , zeros(1,1).*NaN,...
         'indexEnd'                      , zeros(1,1).*NaN,...                                 
         'valueSwitchConditionEnd'       , zeros(1,1).*NaN,...
         'indexReference'                , zeros(1,1).*NaN,...
         'valueSwitchConditionReference' , zeros(1,1).*NaN);

sitToStandSequenceInfo = ...
  struct('switchConditionStart','comSpeed',...
        'switchConditionEnd', 'comHeight',...
        'switchConditionReference','chairForceZ');

%sitToStandDynamicSequence(1) = struct('indexStart', zeros(1,1).*NaN,... 
%                                    'indexEnd', zeros(1,1).*NaN,...
%                                    'indexSeatOff', zeros(1,1).*NaN);
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
    indexPhase = indexListStableStanding(i,1);
    indexStandingStart = ...
       segmentedData.phaseTransitionIndices(indexPhase,1);
    indexSeatOff = 0;

    prevPhase = segmentedData.phaseTransitions(indexPhase,1);
    flag_validTransition = 1;
    if(prevPhase == indexSittingDynamic ...
         || prevPhase == indexSittingStatic)
      flag_validTransition = 0; 
      %Happens when the subjects feet touch the chair force plate.
    end

    %Make sure this sequence goes back to a static sitting sequence
    %without switching between standing sequences ...


    %if(i > 1)
    indexPrevPhase = 1;
    %if(i>1) 
    %  indexListStableStanding(i-1,1);        
    %end

    flag_sittingFound = 0;
    for z = indexPhase:-1:indexPrevPhase 
      if(segmentedData.phaseTransitions(z,1) == indexSittingStatic ...
          || segmentedData.phaseTransitions(z,1) == indexSittingDynamic)
        %Save the last seat off index prior to the standing phase
        if(flag_sittingFound == 0)
          indexSeatOff = segmentedData.phaseTransitionIndices(z,1);
        end          
        flag_sittingFound=1;
      end
    end

    if(flag_sittingFound==0)
      flag_validTransition = -1;
    end

    %end


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




      t0 = segmentedData.phaseTransitionTimes(indexPhase,1);
      t1 = 0;
      idxStart=0;
      indexStandingEnd = 0;
      if(indexPhase < length(segmentedData.phaseTransitionTimes))
        t1                = segmentedData.phaseTransitionTimes(indexPhase+1,1);
        indexStandingEnd  = segmentedData.phaseTransitionIndices(indexPhase+1,1);
      else
        %At the very end of the sequence
        t1 = c3dTime(end);
        indexStandingEnd = length(c3dTime);
      end
      timeQuietStanding = t1-t0;
      %Make sure the standing phase is sufficiently long
      if( (timeQuietStanding) >= quietDwellTime)
        %Go back to the very first transition from static-sitting to 
        %something else
        k=indexPhase;        
        while( k > 1 && segmentedData.phaseTransitions(k,1) ~= indexSittingStatic)          
          k=k-1;
        end
        t1 = segmentedData.phaseTransitionTimes(k,1);
        t0 = 0;
        indexSittingStaticEnd = segmentedData.phaseTransitionIndices(k,1);
        indexSittingStaticStart=0;
        if(k > 1)
          t0 = segmentedData.phaseTransitionTimes(k-1,1);
          indexSittingStaticStart = segmentedData.phaseTransitionIndices(k-1,1);
        else
          t0 = c3dTime(1,1);
          indexSittingStaticStart = 1;
        end
        timeQuietSitting = (t1-t0);
        %Ensure that the length of this static-sitting period exceeds the
        %required threshold
        if((timeQuietSitting) >= quietDwellTime)

          %Refine
          % indexSittingStaticEnd : within tol of the min velocity within the interval
          % indexSeatOff          : within tol of the min chair force within the interval
          % indexStandingStart    : within tol of the max in the interval
          % 
          indexStart      = indexSittingStaticEnd;
          indexEnd        = indexStandingStart;
          indexReference  = indexSeatOff;

          %%
          %
          % Refine the start of the movement: it will begin when the
          %   com velocity is within some tolerance of the minimum com
          %   velocity of the sitting phase
          %%
          minSittingComVelocity = Inf;
          indexMinSittingComVelocity = 0;
          maxSittingComVelocity = -Inf;
          indexMaxSittingComVelocity = 0;
          for k=indexSittingStaticStart:indexSeatOff
            if(sum(isnan(comVelocity(k,:)))==0)
              comSpeedVal = comSpeed(k,1);%norm(comVelocity(k,:));
              if(comSpeedVal < minSittingComVelocity)
                minSittingComVelocity = comSpeedVal;
                indexMinSittingComVelocity = k;
              end
              if(comSpeedVal > maxSittingComVelocity)
                maxSittingComVelocity = comSpeedVal;
                indexMaxSittingComVelocity = k;
              end

            end
          end

          sittingVelocityChange = ...
            (maxSittingComVelocity - minSittingComVelocity);
          startingVelocity = 0.;
          if( sittingVelocityChange < ...
                  startS2SComVelocityTolerance)
            startingVelocity = minSittingComVelocity + 0.5*sittingVelocityChange;
          else
            startingVelocity = minSittingComVelocity + startS2SComVelocityTolerance;  
          end

          flag_threshold = 0;
          comSpeedVal = 0;
          while (flag_threshold == 0 && indexStart > indexMinSittingComVelocity)              
            indexStart = indexStart-1;
            if(sum(isnan(comVelocity(indexStart,:)))==0)
              comSpeedVal = comSpeed(indexStart,1);%norm(comVelocity(indexStart,:));
              if(comSpeedVal <= startingVelocity)
                flag_threshold=1;
              end
            end
          end
          assert(flag_threshold==1);
          valueStart = comSpeedVal;

          %%
          %
          %Refine seat off index: 
          %
          %%
          valueSeatOff = 0;
          if(c3dGrfDataAvailable == 1)
            %%
            % If ground forces are available: seat off will
            % occur when the chair reaches some tolerance value above 
            % the minimum recorded velue
            %%

            fzChairSeated = c3dGrfChair.force(indexSittingStaticEnd:indexStandingEnd,3);

            fzChairStanding = c3dGrfChair.force(indexStandingStart:indexStandingEnd,3);

            %Signal sometimes bounces below zero from the Kistler plates:
            % these values are ignored
            fzChairSeated(   find(fzChairSeated   < 0)) = NaN;
            fzChairStanding( find(fzChairStanding < 0)) = NaN;

            [minFzVal, indexMinFz] = min(fzChairStanding, [], 'omitnan');
            [maxFzVal, indexMaxFz] = max(fzChairStanding, [], 'omitnan');

            %seatOffChangeInForce = maxFzVal-minFzVal;
            %assert(seatOffChangeInForce > seatOffChairForceZTolerance);

            seatOffForce = maxFzVal + seatOffChairForceZTolerance;


            idxSeg = 1;
            idxSegMax = (indexStandingEnd-indexSittingStaticEnd+1);
            while(   idxSeg < idxSegMax ...
                  && ( fzChairSeated(idxSeg) > seatOffForce || isnan(fzChairSeated(idxSeg))==1 ))
              idxSeg=idxSeg+1;
            end
            if(idxSeg >= idxSegMax)
              here=1;
            end
            assert(idxSeg < idxSegMax);
            indexSeatOff = indexSittingStaticEnd+idxSeg-1;
            valueSeatOff = fzChairSeated(idxSeg)-seatOffForce;
          else

            %Seatoff occurs very close to the maximum whole body
            %angular velocity.
            [wnMax, idxSeg]= max(abs(fpeStepLength(...
                                  indexStart:indexStandingStart,1)));

            indexSeatOff = indexStart+idxSeg-1;
            valueSeatOff = ...
              fpeStepLength(indexSeatOff,1);

          end

          indexStandingStart = indexSeatOff;
          %%
          %
          %Refine standing index : standing will occur when the com
          % height is within some tolerance of the maximum value that 
          % occurs during this standing phase.
          %%
          [maxComHeight, indexMaxComHeight] = ...
            max(comPosition(indexStandingStart:indexStandingEnd,3));
          [minComHeight, indexMinComHeight] = ...
            min(comPosition(indexStandingStart:indexStandingEnd,3));

          medComHeight = ...
            median(comPosition(indexStandingStart:indexStandingEnd,3));

          %endS2SComVelocityTolerance

          [maxComSpeed, indexMaxComSpeed] = ...
            max(comSpeed(indexStandingStart:indexStandingEnd,1));
          [minComSpeed, indexMinComSpeed] = ...
            min(comSpeed(indexStandingStart:indexStandingEnd,1));


          %[maxAngSpeed, indexMaxAngSpeed] = ...
          %  max(angSpeed(indexStandingStart:indexStandingEnd,1));
          %[minAngSpeed, indexMinAngpeed] = ...
          %  min(angSpeed(indexStandingStart:indexStandingEnd,1));


          standingHeight = 0;
          standingChangeInHeight = maxComHeight-minComHeight;
          if(standingChangeInHeight < endS2SComHeightTolerance)
            standingHeight = maxComHeight - 0.5*standingChangeInHeight;
          else
            standingHeight = maxComHeight - endS2SComHeightTolerance;
          end

          standingSpeed = 0;
          standingChangeInSpeed = maxComSpeed-minComSpeed;
          if(standingChangeInSpeed < endS2SComVelocityTolerance)
            standingSpeed = minComSpeed + 0.5*standingChangeInSpeed;
          else
            standingSpeed = minComSpeed + endS2SComVelocityTolerance;              
          end


          %standingAngSpeed = 0;
          %standingChangeInAngSpeed = maxAngSpeed-minAngSpeed;
          %if(standingChangeInAngSpeed < endS2SAngularVelocityTolerance)
          %  standingAngSpeed = minAngSpeed + 0.5*standingChangeInAngSpeed;
          %else
          %  standingAngSpeed = minAngSpeed + endS2SAngularVelocityTolerance;                 
          %end

          %Picks off the first time point in which
          % :Vertical velocity is within tolerance
          % :Vertical position is within tolerance of the local max of
          % this sit to stand trial.
          indexStanding=indexStandingStart;

          errRadius = sqrt(abs(comPosition(indexStanding,3)-medComHeight)^2 ...
                         + comVelocity(indexStanding,3)^2);
          errRadiusTolerance = sqrt(endS2SComHeightTolerance^2 ...
                                  + endS2SComVelocityTolerance^2);


          while(indexStanding < indexStandingEnd ...
                && ( comPosition(indexStanding,3) < (medComHeight-endS2SComHeightTolerance) ...
                   || abs(comVelocity(indexStanding,3)) > endS2SComVelocityTolerance ) )
            indexStanding = indexStanding+1;        
            %errRadius = sqrt(abs(comPosition(indexStanding,3)-medComHeight)^2 ...
            %             + comVelocity(indexStanding,3)^2);
          end 

          %if(indexStanding == indexStandingEnd)
          %indexStanding = indexMaxComHeight + indexStandingStart;
          %end

          if(flag_verbose==1 && flag_verbosePlots ==1)
            figure(figVerbose);
            plot(c3dTime(indexStart:indexStandingEnd,1)-c3dTime(indexSeatOff,1),...
                 comPosition(indexStart:indexStandingEnd,3),'b');
            hold on;
            plot(c3dTime(indexStart:indexStandingEnd,1)-c3dTime(indexSeatOff,1),...
                 comVelocity(indexStart:indexStandingEnd,3),'r');
            hold on;              
            plot([1,1]*(c3dTime(indexStanding,1)-c3dTime(indexSeatOff,1)),...
                 [0,1].*comPosition(indexStanding,3),...
                 '-k');
            hold on;

            plot([1,1]*(c3dTime(indexSeatOff,1)-c3dTime(indexSeatOff,1)),...
                 [0,1].*comPosition(indexSeatOff,3),...
                 '-','Color',[1,1,1].*0.5);
            hold on;

            xlabel('Time (s)');
            ylabel('Com Z Pos/Vel');
          end



          valueEnd = comPosition(indexStanding,3);



          %%
          %
          % Check to ensure that both feet remain in contact with the
          % floor before accepting this sit-to-stand
          %
          %%
          flag_bothFeetOnFloor = 1;
          index=indexStart;
          while(flag_bothFeetOnFloor == 1 &&  index  <= indexEnd)
            if(   segmentedData.footContactLeft(index,1) == 0 ...
               || segmentedData.footContactRight(index,1) == 0)
              flag_bothFeetOnFloor = 0;
            end
            index=index+1;
          end

          if(flag_bothFeetOnFloor==1 ...
              || flag_rejectTrialsFootGroundContactBroken == 0)

            if(indexSeatOff <= indexStart)
              here=1;
            end
            if(indexStanding <= indexSeatOff)
              here=1;
            end
            assert(indexSeatOff > indexStart);              
            assert(indexStanding > indexSeatOff);

            sitToStandSequence(countS2Sq).indexStart = indexStart;
            sitToStandSequence(countS2Sq).valueSwitchConditionStart = ... 
              valueStart;

            sitToStandSequence(countS2Sq).indexReference = indexSeatOff;
            sitToStandSequence(countS2Sq).valueSwitchConditionReference = ... 
              valueSeatOff;

            sitToStandSequence(countS2Sq).indexEnd = indexStanding;
            sitToStandSequence(countS2Sq).valueSwitchConditionEnd = ... 
              valueEnd;

            countS2Sq= countS2Sq + 1;
            if(flag_verbose==1)
              phase1duration = c3dTime(indexSeatOff)-c3dTime(indexStart,1);
              phase2duration = c3dTime(indexStanding,1)-c3dTime(indexSeatOff);
              fprintf('\t#%i : sit to stand timing \t[%1.2f, %1.2f]\t(%1.1f, %1.1f, %1.1f), %1.3f: com_max_h-com_med_h)\n',i,...
                       phase1duration,phase2duration,...
                       c3dTime(indexStart,1), c3dTime(indexSeatOff), c3dTime(indexStanding),...
                       maxComHeight-medComHeight);
            end
          else
            if(flag_verbose==1)
              fprintf('\t#%i : Rejected: feet broke ground contact\n',i);
            end
          end

        else


          %sitToStandDynamicSequence(countS2Sd).phaseTransitions = ...
          %  segmentedData.phaseTransitions(k:1:idx,:);
          %sitToStandDynamicSequence(countS2Sd).phaseTransitionTimes = ...
          %  segmentedData.phaseTransitionTimes(k:1:idx,1);
          %sitToStandDynamicSequence(countS2Sd).phaseTransitionIndices = ...
          %  segmentedData.phaseTransitionIndices(k:1:idx,1);              
          %countS2Sd = countS2Sd+1;

          if(flag_verbose==1)
            fprintf(['\t#%i : Rejected: sitting too short (%f)\n'],i,...
                     timeQuietSitting);
          end          
        end

      else
        if(flag_verbose==1)
          fprintf('\t#%i Rejected: standing too short (%f)\n',i,...
                  (timeQuietStanding));
        end
      end      
    end
  end
end

countS2Sq = countS2Sq-1;
%countS2Sd = countS2Sd-1;


fprintf('\t%i: Successful sit-to-stand\n',countS2Sq);
%fprintf('\t%i: Successful dynamic-sit-to-quiet-stand\n',countS2Sd);


% if(flag_writeCsvSequenceFile == 1)
%   idx = strfind(movementSequenceDataFileName,'.mat');
%   movementSequenceDataFileNameCsv = ...
%     [movementSequenceDataFileName(1,1:1:(idx-1)),'.csv'];
%   eventTable = zeros(length(sitToStandSequence),3);
%   for i=1:1:length(sitToStandSequence)
%     eventTable(i,1) = sitToStandSequence(i).indexStart;
%     eventTable(i,2) = sitToStandSequence(i).indexReference;
%     eventTable(i,3) = sitToStandSequence(i).indexEnd;  
%   end
%   csvwrite(movementSequenceDataFileNameCsv, eventTable);
% end

if(flag_verbose==1 && flag_verbosePlots == 1)
  clf(figVerbose);
  close(figVerbose);
end
