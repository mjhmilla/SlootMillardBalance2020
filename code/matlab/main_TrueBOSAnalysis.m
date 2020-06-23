clc;
close all;
clear all;

flag_plotRawData              = 0;
flag_plotFootPrintCopPortrait = 1;
forceLowerBound = 100;

inputDirRelative = '../../inputData/TrueBOS_SubAnalysis';
outputDirRelative = '../../outputData/TrueBOS_SubAnalysis';

pathToBTK='/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk/';
addpath(pathToBTK);

subjectsToProcess = {'configPilot2'}; %'configPilot1'

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
  
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
  
  numberOfTrials = length(inputC3DFiles);

  disp(['Processing: ', subjectId]);  

  %Get the frame print frame offsets
  c3dFileNameAndPath = [inputPath,'/',inputFolder,'/',inputC3DOffsetFile];
  flag_createC3DFilesForRBDL = 0;
  c3dRbdlPlanarSettings=[]; 
  c3dRbdlPath=[];
  flag_MetersRadians=1; 
  flag_grfDataRecorded=1;
  flag_verbose=1;

  [c3dTime, c3dMarkers,c3dMarkerNames, c3dMarkerUnits, c3dForcePlates,...
   c3dForcePlateInfo, c3dGrf, c3dGrfDataAvailable] = ...
      preprocessC3DData(  c3dFileNameAndPath, flag_createC3DFilesForRBDL,...
                      c3dRbdlPlanarSettings, c3dRbdlPath,...                       
                      flag_MetersRadians, flag_grfDataRecorded,...
                      flag_verbose );  
  
  [frameLeftOffset, frameRightOffset]=...
    getFootOffsetFrames(inputOffsetTimeIndex, c3dMarkers,c3dMarkerNames);    
  
  for indexTrial = 1:1:length(inputC3DFiles)
    fprintf('  Trial %i/%i\n',indexTrial, length(inputC3DFiles));
    
    c3dFileNameAndPath = [inputPath,'/',inputFolder,'/',inputC3DFiles{indexTrial}];
    fname = inputC3DFiles{indexTrial};
    idx = strfind(fname,'.c3d');
    pname = fname;
    pname = ['fig_',fname(1:1:idx),'pdf'];  
    plotNameAndPath = [outputPath,'/',outputFolder,'/',pname];
    
    flag_createC3DFilesForRBDL = 0;
    c3dRbdlPlanarSettings=[]; 
    c3dRbdlPath=[];
    flag_MetersRadians=1; 
    flag_grfDataRecorded=1;
    flag_verbose=1;
    
    [c3dTime, c3dMarkers,c3dMarkerNames, c3dMarkerUnits, c3dForcePlates,...
     c3dForcePlateInfo, c3dGrf, c3dGrfDataAvailable] = ...
        preprocessC3DData(  c3dFileNameAndPath, flag_createC3DFilesForRBDL,...
                        c3dRbdlPlanarSettings, c3dRbdlPath,...                       
                        flag_MetersRadians, flag_grfDataRecorded,...
                        flag_verbose );
    
                      
    idxFpLeft = 1;
    idxFpRight = 2;
    
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
    
    if(    norm(rFpCenter'-rFootL) ...
        >  norm(rFpCenter'-rFootR))
      tmp = idxFpLeft;
      idxFpLeft = idxFpRight;
      idxFpRight=tmp;
    end          
    
    
    if(flag_plotFootPrintCopPortrait == 1)
      figure(figFootPrint);
      clf(figFootPrint,'reset');
      
      %Project COP onto the respective foot print frames
      rCopL = zeros(length(c3dTime),3);
      rCopR = zeros(length(c3dTime),3);
      
      c3dMarkersInFootFrame = c3dMarkers;
      
      markerExtents = [NaN,NaN,NaN,NaN];
      
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
          elseif(contains(c3dMarkerNames{m},'R_'))            
            r = frameRight.r;
            E = frameRight.E;
          else
            assert('Marker naming convention broken');
          end
          
          c3dMarkersInFootFrame.(c3dMarkerNames{m})(k,:) = ...
            (E'*(c3dMarkers.(c3dMarkerNames{m})(k,:)'-r))';
        end
      end
      
      [row,col] = find(subPlotPanelIndex==1);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
      subplot('Position',subPlotVec);
        idxL = find(c3dGrf(idxFpLeft).force(:,3)>forceLowerBound);
        if(length(idxL)>0)        
          plot(rCopL(idxL,1), rCopL(idxL,2), 'b');
          hold on;       
          for m=1:1:length(c3dMarkerNames)
            if(isMarkerInSet(c3dMarkerNames{m},markerLeftFootPrefix,markerSet)==1)
              plot(c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxL,2),...
                   'Color',[1,1,1].*0.5);
              hold on;
              text(c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,2),...
                   c3dMarkerNames{m});
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
        title('Left Foot')
        ylabel('Y (cm)');
        xlabel('X (cm)');
        
      [row,col] = find(subPlotPanelIndex==2);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
      subplot('Position',subPlotVec);
        idxR = find(c3dGrf(idxFpRight).force(:,3)>forceLowerBound);
        if(length(idxR)>0)
          plot(rCopR(idxR,1), rCopR(idxR,2), 'r');
          hold on;
          for m=1:1:length(c3dMarkerNames)
            if(isMarkerInSet(c3dMarkerNames{m},markerRightFootPrefix,markerSet)==1)
              plot(c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(idxR,2),...
                   'Color',[1,1,1].*0.5);
              hold on;
              text(c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,1),...
                   c3dMarkersInFootFrame.(c3dMarkerNames{m})(1,2),...
                   c3dMarkerNames{m});
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
        
        
        axis(markerExtents);
        
        
        xticks(xMarks);
        xticklabels(xMarkLabels);
        yticks(yMarks);
        yticklabels(yMarkLabels);
        
        grid on;
        
        [row,col] = find(subPlotPanelIndex==1);          
        subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
        subplot('Position',subPlotVec);
        
        axis(markerExtents);

        xticks(xMarks);
        xticklabels(xMarkLabels);
        yticks(yMarks);
        yticklabels(yMarkLabels);        
        grid on;

        
        configPlotExporter;
        print('-dpdf',plotNameAndPath);  
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
  
end
