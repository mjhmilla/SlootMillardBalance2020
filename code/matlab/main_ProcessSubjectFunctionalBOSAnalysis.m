clc;
close all;
clear all;

flag_plotRawData               = 0;
flag_plotFootPrintCopPortrait  = 1;
flag_plotConvexHullOfAllTrials = 1;

forceLowerBound = 0.5; %As in half body weight
barefootPadCompressionMax = (2.0*max(25*0.5/1000,15*0.5/1000)); %A high limit to detect clear lift off
shodPadCompressionMax     = barefootPadCompressionMax+5/1000;  

flag_limitFootAngle = 1;

lowPassFilterFrequencyHz = -1; %as in 10 Hz

minimumNumberOfPoints = 150;

%Accept 95% of the points, reject the 5% that are outliers
fractionOfPointsToAccept = 0.95; 

                       
inputDirRelative = '../../inputData/TrueBOS_SubAnalysis';
outputDirRelative = '../../outputData/TrueBOS_SubAnalysis';

pathToBTK='/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk/';
addpath(pathToBTK);


subjectsToProcess = { 'configP07',...
                      'configP08',...
                      'configP09',...
                      'configP10'};
% subjectsToProcess = { 'configPilot1',...
%                       'configPilot2',...
%                       'configP01',...
%                       'configP02',...
%                       'configP03',...
%                       'configP04',...
%                       'configP05',...
%                       'configP06',...
%                       'configP07',...
%                       'configP08',...
%                       'configP09',...
%                       'configP10'};
 

%From inputData/COP_SubAnalysis/main_lframe.m : avg
tmp = load(['../../inputData/COP_SubAnalysis/',...
                    'errorForcePlateAtIndex1.mat']);
fpAtIndex1ErrorStruct = tmp.errorStruct;                  

tmp = load(['../../inputData/COP_SubAnalysis/',...
                    'errorForcePlateAtIndex2.mat']);
fpAtIndex2ErrorStruct = tmp.errorStruct;                  
%%
%Plot settings
%%
lineWidth  = 0.75;
boxWidth   = 0.33;
panelHeight = 20;
panelWidth  = 10;

numberOfFiguresPerPage        = 2;
numberOfVerticalPlotRows      = 1;
numberOfHorizontalPlotColumns = numberOfFiguresPerPage;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;           
pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;           
plotHeight  = panelHeight;
plotWidth   = panelWidth ;

plotConfigGeneric;


%%
% Set up the input/output directory structure
%%
codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);

%%
%Output figures
%%
figFootPrint = figure;



 

for indexSubject = 1:1:length(subjectsToProcess)
  
  copRightBarefoot        = [];%withShoes == 0
  copLeftBarefoot         = [];
  
  copLeftWithShoes        = [];%withShoes == 1
  copRightWithShoes       = [];
  
  copLeftWithHikingShoes  = [];%withShoes == 2
  copRightWithHikingShoes = [];
  
  

  
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
  
  numberOfTrials = length(inputC3DFiles);

  footData = struct('isValid',0,'r0F0',[],'ea321',[],'rFCF',[],'f0',[],...
      'bos',[],'markers',[]);
  
  trialData(numberOfTrials) = ...
    struct('isValid',0,'withShoes',0,...
      'midFootLength',0,'footWidth',0,...
      'footLength',0,'footCompressionMax',0,...
      'angleXUB',0,'angleYUB',0,...
      'leftFoot',footData,'rightFoot',footData);
  
  disp(['Processing: ', subjectId]);  
  
  frameLeftOffset     = [];
  frameRightOffset    = [];
  frameLeft           = [];
  frameRight          = [];
  functionalBosLeft   = [];
  functionalBosRight  = [];

  midFootLength = 0;
  footWidth = 0;
  footLength = 0;
  footCompressionMax = 0;
  angleXUB = 0;
  angleYUB = 0;
  
  
  for indexTrial = 1:1:length(inputC3DFiles)

    trialData(indexTrial).isValid = trialIsValid(indexTrial,1);
    trialData(indexTrial).withShoes= withShoes(indexTrial,1);
    
    if(trialIsValid(indexTrial,1)==1)
      
      %%
      % Update the foot reference frame
      %%
      if(updateFootFrames(indexTrial,1)==1 )

        c3dFileNameAndPathRef = [inputPath,'/',inputFolder,...
          '/',inputC3DFiles{updateFootFrames(indexTrial,2)}];
        
        flag_createC3DFilesForRBDLRef = 0;
        c3dRbdlPlanarSettingsRef=[]; 
        c3dRbdlPathRef=[];
        flag_MetersRadiansRef=1; 
        flag_grfDataRecordedRef=1;
        flag_verboseRef=1;
                
        [c3dTimeRef, c3dMarkersRef,c3dMarkerNamesRef, c3dMarkerUnitsRef, c3dForcePlatesRef,...
         c3dForcePlateInfoRef, c3dGrfRef, c3dGrfDataAvailableRef] = ...
            preprocessC3DData(  c3dFileNameAndPathRef, ...
                            fpAtIndex1ErrorStruct,...
                            fpAtIndex2ErrorStruct,...
                            lowPassFilterFrequencyHz,...
                            flag_createC3DFilesForRBDLRef,...
                            c3dRbdlPlanarSettingsRef, c3dRbdlPathRef,...                       
                            flag_MetersRadiansRef, flag_grfDataRecordedRef,...
                            flag_verboseRef );          
          
        [frameLeftOffset, frameRightOffset]=...
          getFootOffsetFrames(updateFootFrames(indexTrial,3), ...
                                c3dMarkersRef,c3dMarkerNamesRef); 

        [frameLeft,frameRight] = ...
          getFootFrames(updateFootFrames(indexTrial,3),c3dMarkersRef,...
            c3dMarkerNamesRef,frameLeftOffset,frameRightOffset);  

         [functionalBosLeft, functionalBosRight] = ...
          getFunctionalFootBaseOfSupport2(updateFootFrames(indexTrial,3),...
            c3dMarkersRef, c3dMarkerNamesRef, frameLeft, frameRight);             

      end
      
      %%
      % Update the size of the foot
      %%
      if(updateFootScale(indexTrial,1)==1 )

        c3dFileNameAndPathRef = [inputPath,'/',inputFolder,...
          '/',inputC3DFiles{updateFootScale(indexTrial,2)}];
        
        flag_createC3DFilesForRBDLRef = 0;
        c3dRbdlPlanarSettingsRef=[]; 
        c3dRbdlPathRef=[];
        flag_MetersRadiansRef=1; 
        flag_grfDataRecordedRef=1;
        flag_verboseRef=1;
                
        [c3dTimeRef, c3dMarkersRef,c3dMarkerNamesRef, c3dMarkerUnitsRef, c3dForcePlatesRef,...
         c3dForcePlateInfoRef, c3dGrfRef, c3dGrfDataAvailableRef] = ...
            preprocessC3DData(  c3dFileNameAndPathRef, ...
                            fpAtIndex1ErrorStruct,...
                            fpAtIndex2ErrorStruct,...
                            lowPassFilterFrequencyHz,...
                            flag_createC3DFilesForRBDLRef,...
                            c3dRbdlPlanarSettingsRef, c3dRbdlPathRef,...                       
                            flag_MetersRadiansRef, flag_grfDataRecordedRef,...
                            flag_verboseRef );          

        [frameLeftOffsetTemp, frameRightOffsetTemp]=...
          getFootOffsetFrames(updateFootScale(indexTrial,3), ...
                                c3dMarkersRef,c3dMarkerNamesRef); 

        [frameLeftTemp,frameRightTemp] = ...
          getFootFrames(updateFootScale(indexTrial,3),c3dMarkersRef,...
            c3dMarkerNamesRef,frameLeftOffsetTemp,frameRightOffsetTemp);  

        idxRefTemp = updateFootScale(indexTrial,3);
        [footLength, midFootLength, footWidth] = ...
         getFootSize(idxRefTemp, frameLeftTemp, frameRightTemp, c3dMarkersRef);
         
         
        footCompressionMax = 0;
        switch(withShoes(indexTrial,1))
          case 0
            footCompressionMax = barefootPadCompressionMax;
          case 1
            footCompressionMax = shodPadCompressionMax;
          case 2
            footCompressionMax = shodPadCompressionMax;
          otherwise
            assert(0)            
        end
         
        %Doubled to allow for foot flexibility
        angleXUB = atan2(footCompressionMax,midFootLength);
        angleYUB = atan2(footCompressionMax,footWidth);

        if(flag_limitFootAngle == 0)
          footCompressionMax = Inf;
          angleXUB = Inf;
          angleYUB = Inf;
        end
         
      end      


      
      
      
      fprintf('  Trial %i/%i\n',indexTrial, length(inputC3DFiles));

      c3dFileNameAndPath = [inputPath,'/',inputFolder,'/',inputC3DFiles{indexTrial}];
      fname = inputC3DFiles{indexTrial};
      idx = strfind(fname,'.c3d')-1;
      pname = fname;
      pname = ['fig_',fname(1:1:idx)];  
      plotNameAndPath = [outputPath,'/',outputFolder,'/',pname];

      flag_createC3DFilesForRBDL = 0;
      c3dRbdlPlanarSettings=[]; 
      c3dRbdlPath=[];
      flag_MetersRadians=1; 
      flag_grfDataRecorded=1;
      flag_verbose=1;

      [c3dTime, c3dMarkers,c3dMarkerNames, c3dMarkerUnits, c3dForcePlates,...
       c3dForcePlateInfo, c3dGrf, c3dGrfDataAvailable] = ...
          preprocessC3DData(  c3dFileNameAndPath, ...
                          fpAtIndex1ErrorStruct,...
                          fpAtIndex2ErrorStruct,...
                          lowPassFilterFrequencyHz,...
                          flag_createC3DFilesForRBDL,...
                          c3dRbdlPlanarSettings, c3dRbdlPath,...                       
                          flag_MetersRadians, flag_grfDataRecorded,...
                          flag_verbose );


      idxFpLeft   = indexLeftForcePlate(indexTrial,1);
      idxFpRight  = indexRightForcePlate(indexTrial,1);

      rFootL = zeros(1,3);
      rFootR = zeros(1,3);

      countFootL=0;
      countFootR=0;

      %%
      %Identify the index of the force plate that the left and right foot
      %is placed on.
      %%
      for m=1:1:length(c3dMarkerNames)
        if(contains(c3dMarkerNames{m},'L_')==1)
          rFootL = rFootL + c3dMarkers.(c3dMarkerNames{m})(1,:);
          countFootL = countFootL+1;
        else
          rFootR = rFootR + c3dMarkers.(c3dMarkerNames{m})(1,:);      
          countFootR = countFootR+1;
        end
      end
      rFootL = rFootL./countFootL;
      rFootR = rFootR./countFootR;

      rFpCenter = mean(c3dForcePlates(idxFpLeft).corners,2);

%       if(    norm(rFpCenter'-rFootL) ...
%           >  norm(rFpCenter'-rFootR))
%         tmp = idxFpLeft;
%         idxFpLeft = idxFpRight;
%         idxFpRight=tmp;
%       end          


      if(flag_plotFootPrintCopPortrait == 1)
        figure(figFootPrint);
        clf(figFootPrint,'reset');
      end

      %Project COP onto the respective foot print frames
      rCopL = zeros(length(c3dTime),3);
      rCopR = zeros(length(c3dTime),3);

      c3dMarkersInFootFrame = c3dMarkers;

      markerExtents = [NaN,NaN,NaN,NaN];

      r0L0 = zeros(length(c3dTime),3);
      ea321L = zeros(length(c3dTime),3);
      fL = zeros(length(c3dTime),3);
      r0R0 = zeros(length(c3dTime),3);
      ea321R = zeros(length(c3dTime),3);
      fR = zeros(length(c3dTime),3);

      for k=1:1:length(c3dTime)
        [frameLeft,frameRight] = getFootFrames(k,c3dMarkers,c3dMarkerNames,...
                                          frameLeftOffset,frameRightOffset);



        rCopL(k,:) = ( frameLeft.E'*(c3dGrf(idxFpLeft).cop(k,:)'...
                                     -frameLeft.r))';
        rCopR(k,:) = ( frameRight.E'*(c3dGrf(idxFpRight).cop(k,:)'...
                                      -frameRight.r))';

        for m=1:1:length(c3dMarkerNames)
          r = zeros(3,1);
          E = zeros(3,3);
          if(contains(c3dMarkerNames{m},'L_'))
            r = frameLeft.r;
            E = frameLeft.E;

            r0L0(k,:) = r;
            ea321L(k,:) = calcEA321(E);
            fL(k,:) = ( E*(c3dGrf(idxFpLeft).force(k,:)') )';

          elseif(contains(c3dMarkerNames{m},'R_'))            
            r = frameRight.r;
            E = frameRight.E;

            r0R0(k,:) = r;
            ea321R(k,:) = calcEA321(E);
            fR(k,:) = ( E*(c3dGrf(idxFpRight).force(k,:)') )';

          else
            if(k==1)
              disp(['Warning, extra marker: ', c3dMarkerNames{m}]);             
            end
            %assert('Marker naming convention broken');
          end

          c3dMarkersInFootFrame.(c3dMarkerNames{m})(k,:) = ...
            (E'*(c3dMarkers.(c3dMarkerNames{m})(k,:)'-r))';
        end
      end

      if(flag_plotFootPrintCopPortrait == 1)
        [row,col] = find(subPlotPanelIndex==1);          
        subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
        subplot('Position',subPlotVec);
      end
      
      trialData(indexTrial).midFootLength = midFootLength;
      trialData(indexTrial).footLength = footLength;      
      trialData(indexTrial).footWidth  = footWidth;
      trialData(indexTrial).footCompressionMax = footCompressionMax;
      trialData(indexTrial).angleXUB = angleXUB;
      trialData(indexTrial).angleYUB = angleYUB;
      
      if( sum(c3dGrf(idxFpLeft).force(:,3)>forceLowerBound*bodyWeight) > 0 )
        trialData(indexTrial).leftFoot.isValid = 1;
        trialData(indexTrial).leftFoot.r0F0   = r0L0;
        trialData(indexTrial).leftFoot.ea321 = ea321L;
        trialData(indexTrial).leftFoot.rFCF   = rCopL;
        trialData(indexTrial).leftFoot.f0    = c3dGrf(idxFpLeft).force;
        trialData(indexTrial).leftFoot.bos   = functionalBosLeft;
        
        for m=1:1:length(c3dMarkerNames)
          if(contains(c3dMarkerNames{m},'L_'))  
            trialData(indexTrial).leftFoot.markers.(c3dMarkerNames{m}) = ...
              c3dMarkersInFootFrame.(c3dMarkerNames{m});
          end
        end            

      end
      
      idxL = find(c3dGrf(idxFpLeft).force(:,3)>forceLowerBound*bodyWeight ...
                  & abs(ea321L(:,3)) < angleXUB ...
                  & abs(ea321L(:,2)) < angleYUB);      
      
      if(length(idxL)>minimumNumberOfPoints)
            
        %Remove the 5% least dense points
        [idxLCH,idxLEx] = removeCenterOfPressureOutliers( rCopL(idxL,1:2)./[footWidth footLength],...
                                                 fractionOfPointsToAccept);
        
                                               
        idxLAll = idxL;
        
        idxLMain = idxL(idxLCH);
        idxLEx   = idxL(idxLEx);
        idxL     = idxLMain;
                                                     
        chL = convhull(rCopL(idxL,1), rCopL(idxL,2));  
        if(flag_plotFootPrintCopPortrait == 1)
          fill(functionalBosLeft(:,1),functionalBosLeft(:,2),[1,1,1].*0.75,...
               'EdgeColor','none');
          hold on;          

          plot(rCopL(idxLAll,1), rCopL(idxLAll,2), 'b');
          hold on;      

          plot(rCopL(idxLEx,1),rCopL(idxLEx,2),'oc');
          hold on;           
                  
          plot(rCopL(idxL(chL),1),rCopL(idxL(chL),2),'k','LineWidth',2);
          hold on;       
               
          
        end

        if(    min(rCopL(idxL(chL),1)) > -0.10 && max(rCopL(idxL(chL),1))< 0.10 ...
            && min(rCopL(idxL(chL),2)) > -0.15 && max(rCopL(idxL(chL),2))< 0.25)
          
          switch (withShoes(indexTrial,1))            
            case 0
              copLeftBarefoot         = [copLeftBarefoot;rCopL(idxL(chL),:)];
            case 1
              copLeftWithShoes        = [copLeftWithShoes;rCopL(idxL(chL),:)];
            case 2
              copLeftWithHikingShoes  = [copLeftWithHikingShoes;rCopL(idxL(chL),:)];
            otherwise
              assert(0);
          end
        end

        if(flag_plotFootPrintCopPortrait == 1)
          for m=1:1:length(c3dMarkerNames)
            if(isMarkerInSet(c3dMarkerNames{m},markerLeftFootPrefix,markerSet)==1)
              plot(c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,2),...
                   'Color',[1,1,1].*0.5);
              hold on;
              markerLabel = c3dMarkerNames{m};
              idx = strfind(markerLabel,'_');
              markerLabel(1,idx) = ' ';

              text(c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,2),...
                   markerLabel);
              hold on;
              markerExtents(1,1) = min([markerExtents(1,1); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,1)]);
              markerExtents(1,2) = max([markerExtents(1,2); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,1)]);
              markerExtents(1,3) = min([markerExtents(1,3); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,2)]);
              markerExtents(1,4) = max([markerExtents(1,4); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,2)]);

            end
          end
          axis equal
          box off;
        end
      end
      if(flag_plotFootPrintCopPortrait == 1)
        title('Left Foot')
        ylabel('Y (cm)');
        xlabel('X (cm)');
      end

      if(flag_plotFootPrintCopPortrait == 1)
        [row,col] = find(subPlotPanelIndex==2);          
        subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
        subplot('Position',subPlotVec);
      end
      
      if( sum(c3dGrf(idxFpRight).force(:,3)>forceLowerBound*bodyWeight) > 0 )
        trialData(indexTrial).rightFoot.isValid = 1;
        trialData(indexTrial).rightFoot.r0F0    = r0R0;
        trialData(indexTrial).rightFoot.ea321   = ea321R;
        trialData(indexTrial).rightFoot.rFCF    = rCopR;
        trialData(indexTrial).rightFoot.f0      = c3dGrf(idxFpRight).force;
        trialData(indexTrial).rightFoot.bos     = functionalBosRight;

    
        for m=1:1:length(c3dMarkerNames)
          if(contains(c3dMarkerNames{m},'R_'))  
            trialData(indexTrial).rightFoot.markers.(c3dMarkerNames{m}) = ...
              c3dMarkersInFootFrame.(c3dMarkerNames{m});
          end
        end         
      end
      
      idxR = find(c3dGrf(idxFpRight).force(:,3)>forceLowerBound*bodyWeight ...
                  & abs(ea321R(:,3)) < angleXUB ...
                  & abs(ea321R(:,2)) < angleYUB);          
                
      if(length(idxR)>minimumNumberOfPoints)
        [idxRCH,idxRex] = removeCenterOfPressureOutliers(rCopR(idxR,1:2)./[footWidth footLength],...
                                                fractionOfPointsToAccept);
        idxRAll  = idxR;
        idxRMain = idxR(idxRCH);
        idxREx   = idxR(idxRex);
        idxR     = idxRMain;
                                                     
        chR = convhull(rCopR(idxR,1), rCopR(idxR,2));  
        if(flag_plotFootPrintCopPortrait == 1)
          fill(functionalBosRight(:,1),functionalBosRight(:,2),[1,1,1].*0.75,...
               'EdgeColor','none');
          hold on;                     
          plot(rCopR(idxRAll,1), rCopR(idxRAll,2), 'r');
          hold on;
                           
          plot(rCopR(idxREx,1),rCopR(idxREx,2),'oc');
          hold on;
          
          plot(rCopR(idxR(chR),1),rCopR(idxR(chR),2),'k','LineWidth',2);
          hold on;
        end

        if(    min(rCopR(idxR(chR),1)) > -0.10 && max(rCopR(idxR(chR),1))< 0.10 ...
            && min(rCopR(idxR(chR),2)) > -0.15 && max(rCopR(idxR(chR),2))< 0.25)

          switch (withShoes(indexTrial,1))            
            case 0
              copRightBarefoot         = [copRightBarefoot;rCopR(idxR(chR),:)];
            case 1
              copRightWithShoes        = [copRightWithShoes;rCopR(idxR(chR),:)];
            case 2
              copRightWithHikingShoes  = [copRightWithHikingShoes;rCopR(idxR(chR),:)];
            otherwise
              assert(0);
          end
                  
        end          

        if(flag_plotFootPrintCopPortrait == 1)
          for m=1:1:length(c3dMarkerNames)
            if(isMarkerInSet(c3dMarkerNames{m},markerRightFootPrefix,markerSet)==1)
              plot(c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,2),...
                   'Color',[1,1,1].*0.5);
              hold on;
              markerLabel = c3dMarkerNames{m};
              idx = strfind(markerLabel,'_');
              markerLabel(1,idx) = ' ';                
              text(c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,2),...
                   markerLabel);
              hold on;            
              markerExtents(1,1) = min([markerExtents(1,1); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,1)]);
              markerExtents(1,2) = max([markerExtents(1,2); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,1)]);
              markerExtents(1,3) = min([markerExtents(1,3); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,2)]);
              markerExtents(1,4) = max([markerExtents(1,4); ...
                c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,2)]);
            end

          end      
          axis equal;
          box off;
        end
      end
          
      if(flag_plotFootPrintCopPortrait == 1)
        title('Right Foot')
        ylabel('Y (cm)');
        xlabel('X (cm)');

        markerExtents = markerExtents + [-1, 1, -1, 1].*0.01;

        xDelta = 0.01;
        xStart = -floor(abs(markerExtents(1,1)/xDelta))*xDelta;
        xEnd   =  floor(abs(markerExtents(1,2)/xDelta))*xDelta;        
        xMarks = [xStart:xDelta:xEnd];

        xMarkLabels = cell(size(xMarks));
        for i=1:1:length(xMarks)
          xMarkLabels{1,i} = sprintf('%d',round(xMarks(1,i)*100));
        end

        yDelta = 0.01;
        yStart = -floor(abs(markerExtents(1,3)/xDelta))*yDelta;
        yEnd   =  floor(abs(markerExtents(1,4)/xDelta))*yDelta;        
        yMarks = [yStart:yDelta:yEnd];

        yMarkLabels = cell(size(yMarks));
        for i=1:1:length(yMarks)
          yMarkLabels{1,i} = sprintf('%d',round(yMarks(1,i)*100));
        end

        %if(sum(isnan(markerExtents))>0)
        %  markerExtents = [-0.1,0.1,-0.1,0.2];
        %end

        %axis(markerExtents);


        xticks(xMarks);
        xticklabels(xMarkLabels);
        yticks(yMarks);
        yticklabels(yMarkLabels);

        grid on;

        [row,col] = find(subPlotPanelIndex==1);          
        subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
        subplot('Position',subPlotVec);

        %axis(markerExtents);

        xticks(xMarks);
        xticklabels(xMarkLabels);
        yticks(yMarks);
        yticklabels(yMarkLabels);        
        grid on;


        configPlotExporter;
        print('-djpeg',[plotNameAndPath,'.jpeg']);
        print('-dpdf',[plotNameAndPath,'.pdf']);
      end

      %%
      %
      %%
      if(flag_plotRawData == 1)
        figTest = figure;

        for k=1:100:length(c3dTime)
          clf(figTest,'reset');
          for m=1:1:length(c3dMarkerNames)
            markerColor = 'r';
            if(contains(c3dMarkerNames{m},'L_')==1)
              markerColor = 'b';
            end

            plot3(c3dMarkers.(c3dMarkerNames{m})(k,1),...
                  c3dMarkers.(c3dMarkerNames{m})(k,2),...
                  c3dMarkers.(c3dMarkerNames{m})(k,3),...
                  'o','Color',markerColor);
            hold on;
            text(c3dMarkers.(c3dMarkerNames{m})(k,1),...
                  c3dMarkers.(c3dMarkerNames{m})(k,2),...
                  c3dMarkers.(c3dMarkerNames{m})(k,3),...
                  c3dMarkerNames{m});
            hold on;         
          end



          for f=1:1:length(c3dGrf)
            fstart = c3dGrf(f).cop(k,:);
            fend   = c3dGrf(f).cop(k,:)-c3dGrf(f).force(k,:).*0.001;
            fvec = [fstart;fend];
            fColor = 'r';
            if(f==2)
              fColor ='b';
            end

            plot3(fvec(:,1),fvec(:,2),fvec(:,3),fColor);
            hold on;
            plot3(fvec(1,1),fvec(1,2),fvec(1,3),['x',fColor(:)]);
            hold on;

            fpSquare = c3dForcePlates(f).corners';
            fpSquare = [fpSquare;fpSquare(1,:)];

            plot3(fpSquare(:,1),fpSquare(:,2),fpSquare(:,3),...
                  ['-',fColor(:)]);              
            hold on;
          end




          [frameLeft,frameRight] = getFootFrames(k,c3dMarkers,c3dMarkerNames,...
                                            frameLeftOffset,frameRightOffset);

          plot3(frameLeft.r(1,1),frameLeft.r(2,1),frameLeft.r(3,1),'xb');
          hold on;
          vecColor = ['r','g','b'];
          for a=1:1:size(frameLeft.E,2)
            vec = [frameLeft.r(:,1)';(frameLeft.r(:,1)'+0.1.*frameLeft.E(:,a)')];
            plot3(vec(:,1),vec(:,2),vec(:,3),vecColor(1,a));
            hold on;
          end

          plot3(frameRight.r(1,1),frameRight.r(2,1),frameRight.r(3,1),'xr');
          hold on;
          vecColor = ['r','g','b'];
          for a=1:1:size(frameRight.E,2)
            vec = [frameRight.r(:,1)';(frameRight.r(:,1)'+0.1.*frameRight.E(:,a)')];
            plot3(vec(:,1),vec(:,2),vec(:,3),vecColor(1,a));
            hold on;
          end

          xlabel('X');
          ylabel('Y');
          zlabel('Z');

          grid on;
          axis equal;

          hold on;
        end
      end

    
    end

    here=1;

  end
  

  subjectDataFile =[outputPath,'/',outputFolder,'/data',subjectId,'.mat'];
  save( subjectDataFile, 'trialData' );

  assert(length(trialData) <= length(withShoes));
  clear('trialData');
  
  if(flag_plotConvexHullOfAllTrials==1)
    figBos = figure;

    [row,col] = find(subPlotPanelIndex==1);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);  

    if(isempty(copLeftWithHikingShoes)==0)
      hL  = convhull(copLeftWithHikingShoes(:,1),copLeftWithHikingShoes(:,2));      
    end
    sL  = convhull(copLeftWithShoes(:,1),copLeftWithShoes(:,2));  
    bL  = convhull(copLeftBarefoot(:,1),copLeftBarefoot(:,2));

    if(isempty(copLeftWithHikingShoes)==0)
      plot(copLeftWithHikingShoes(hL,1),copLeftWithHikingShoes(hL,2),...
        'Color',[1,1,1].*0.75,'LineWidth',2);
      hold on;    
    end
    plot(copLeftWithShoes(sL,1),copLeftWithShoes(sL,2),...
      'Color',[1,1,1].*0.5);
    hold on;
    plot(copLeftBarefoot(bL,1),copLeftBarefoot(bL,2),...
      'Color',[1,1,1].*0.);
    hold on;

    markerExtents = [-0.05,0.05,-0.05,0.18]; 

    markerExtents = markerExtents;

    xDelta = 0.01;
    xStart = -floor(abs(markerExtents(1,1)/xDelta))*xDelta;
    xEnd   =  floor(abs(markerExtents(1,2)/xDelta))*xDelta;        
    xMarks = [xStart:xDelta:xEnd];

    xMarkLabels = cell(size(xMarks));
    for i=1:1:length(xMarks)
      xMarkLabels{1,i} = sprintf('%d',round(xMarks(1,i)*100));
    end

    yDelta = 0.01;
    yStart = -floor(abs(markerExtents(1,3)/xDelta))*yDelta;
    yEnd   =  floor(abs(markerExtents(1,4)/xDelta))*yDelta;        
    yMarks = [yStart:yDelta:yEnd];

    yMarkLabels = cell(size(yMarks));
    for i=1:1:length(yMarks)
      yMarkLabels{1,i} = sprintf('%d',round(yMarks(1,i)*100));
    end


    axis(markerExtents);

    xticks(xMarks);
    xticklabels(xMarkLabels);
    yticks(yMarks);
    yticklabels(yMarkLabels);

    grid on;


    box off;
    axis equal;
    if(isempty(copLeftWithHikingShoes)==0)
      legend('Hiking shoes','Running shoes','Barefoot');
    else
      legend('Running shoes','Barefoot');      
    end

    xlabel('X (m)');
    ylabel('Y (m)');
    title('Left Foot');

    [row,col] = find(subPlotPanelIndex==2);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);  

    if(isempty(copRightWithHikingShoes)==0)
      hR  = convhull(copRightWithHikingShoes(:,1),copRightWithHikingShoes(:,2));      
    end
    sR  = convhull(copRightWithShoes(:,1),copRightWithShoes(:,2));  
    bR  = convhull(copRightBarefoot(:,1),copRightBarefoot(:,2));

    if(isempty(copRightWithHikingShoes)==0)
      plot(copRightWithHikingShoes(hR,1),copRightWithHikingShoes(hR,2),...
        'Color',[1,1,1].*0.75,'LineWidth',2);
      hold on;    
    end
    plot(copRightWithShoes(sR,1),copRightWithShoes(sR,2),...
      'Color',[1,1,1].*0.5);
    hold on;
    plot(copRightBarefoot(bR,1),copRightBarefoot(bR,2),...
      'Color',[1,1,1].*0.);
    hold on;

    axis(markerExtents);

    xticks(xMarks);
    xticklabels(xMarkLabels);
    yticks(yMarks);
    yticklabels(yMarkLabels);

    grid on;

    box off;
    axis equal;
    if(isempty(copRightWithHikingShoes)==0)
      legend('Hiking shoes','Running shoes','Barefoot');
    else
      legend('Running shoes','Barefoot');      
    end
    

    xlabel('X (cm)');
    ylabel('Y (cm)');
    title('Right Foot');

    netBosPlot=[outputPath,'/',outputFolder,'/bosBarefootShoes'];


    configPlotExporter;
    print('-djpeg',[netBosPlot,'.jpeg']);    
    print('-dpdf',[netBosPlot,'.pdf']);
  end
  
end
