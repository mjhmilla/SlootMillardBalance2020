clc;
close all;
clear all;

%List of the subjects to process
subjectsToProcess = ...
 {'configE01','configE02','configE03','configE05', 'configE06',...
  'configE07','configE08',...
  'configH01','configH02','configH03','configH04','configH05',...
  'configH06','configH07','configH08','configH09','configH10'};


pathToBTK  = '/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk/';
addpath(pathToBTK);

%The model building tool in the RBDL tool chain
pathToModelFactoryRoot ='/home/mjhmilla/dev/ModelFactory';
pathToModelFactory = genpath(pathToModelFactoryRoot);
addpath(pathToModelFactory);

%%
% Setup the input/output directory structure
%%
inputDirRelative = '../../inputData';
outputDirRelative = '../../outputData';

codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);

%%
% Global Variables
%%

gravityVector = [0;0;-9.81];
contactPlanes = [0,0,0];  

%Motion segmentation constants
%
% The segmentData script classfies every data point 
% into one of these six states:
%
% sitting-static
%   on the chair and the fpe close to the CoM ground proj.
% sitting-dynamic
%   on the chair and the fpe far from the CoM ground proj.
%
% crouching-stable
%   off the chair, legs bent, fpe within convex hull of feet
% crouching-unstable
%   off the chair, legs bent, fpe outside convex hull of feet
%
% standing-stable
%   off the chair, legs straight, fpe inside convex hull of feet
% standing-unstable
%   off the chair, legs straight, fpe outside convex hull of feet
%

%This means that above a height of comZMax - 0.25*(comZMax-comZMin)
%is considered standing. For subject E01 standing happens at 76 cm while
%the maximum CoM height is 82 cm.
thresholdNormDistanceBelowStanding     = 0.25;         


%Foot
footContactZMovementTolerance  = 0.01; %m

%%
%Motion sequence constants
%
% This code goes through the segmented data and identifies motion 
% sequences:
%
% quiet-sit-to-quiet-stand: subject progresses from
%   -sitting-static for thresholdQuietTime
%   -standing-stable for thresholdQuietTim
%   -s.t. once seat-off occurs the subject does not contact the seat
%         again until (eventually) a standing-stable pose is reached.
%%
thresholdQuietTime = 0.5;

%%
% Processing Flags
%%


%Preprocessing of C3D data
flag_loadC3DMatFileData         = 1;
flag_useMetersRadiansInC3DData  = 1;
flag_writeC3DDataForMeshup      = 0;
flag_verbose                    = 0;

%FPE processing
flag_loadFpeDataFromFile   = 0;
flag_writeFpeDataForMeshup = 0;

%Capture point processing
flag_loadCapDataFromFile   = 0;
flag_writeCapDataForMeshup = 0;

%Calc distance between key ground points and the convex hull of the feet.
%Here the key ground points are: CoM ground projection, CoP, Fpe, Cap
flag_loadKeyPointDistanceToFootConvexHull = 0;

%Motion segmentation
flag_loadSegmentedMotionDataFromFile = 0;

%Motion sequence identification
flag_loadMovementSequenceFromFile =0;


%Plotting
flag_visualize           = 0;
    numberOfFramesToDraw = 50;
    scaleForceToDistance = 1/1000;

%Writing of files for RBDL's visualization tools (puppeteer/meshup)
flag_writeRBDLToolChainFiles = 1;


for indexSubject = 1:1:length(subjectsToProcess)
  
  %%
  % 1. Configure the list of input/output files for each trial
  %%
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
  

  
  numberOfTrials = length(inputC3DFiles);

  disp(['Processing: ', subjectId]);

  %Echo the input files, folders, and output folders
  disp('  file and folder layout');
  for i=1:1:length(inputC3DFolders)
    tmp = inputC3DFolders{i};
    idx0 = strfind(tmp,'/');
    inFolder = tmp((idx0(end-2):1:idx0(end)));
    
    inFile   = inputC3DFiles{i};    
    
    tmp = outputTrialFolders{i};
    idx0 = strfind(tmp,'/');
    outFolder = tmp((idx0(end-2):1:idx0(end)));
    
    if(length(inFile) < 30)
      for k=length(inFile):1:30
        inFile = [inFile,' '];
      end
    end
    
    fprintf('    %s\t%s\n',inFile,outFolder);
    here=1;
  end  
  
  %%
  %Process each trial
  %%

  disp('  Processing data and metrics:');
  for indexTrial = 1:1:numberOfTrials
    
    disp(['    :',inputC3DFiles{indexTrial}]);
        
    
    inputC3DFolder = inputC3DFolders{indexTrial};
    c3dFileName         = inputC3DFiles{indexTrial};
    

    inputV3DFolder = inputV3DFolders{indexTrial};    
    wholeBodyFileName   = inputWholeBodyFiles{indexTrial};
    
    inputAnthroFolder  = inputV3DFolder;    
    anthroFileName     = inputAnthroFiles{indexTrial};
    
    if(indexTrial == index_Static)
      if(indexTrial < numberOfTrials)
        inputAnthroFolder = inputV3DFolders{indexTrial+1};
        anthroFileName     = inputAnthroFiles{indexTrial+1};
      else
        inputAnthroFolder = inputV3DFolders{indexTrial-1};
        anthroFileName     = inputAnthroFiles{indexTrial-1};
      end
    end    
    
    
    outputTrialFolder     = outputTrialFolders{indexTrial};    
    meshupGrfFileName     = outputMeshupGrfFiles{indexTrial};
    meshupFpeFileName     = outputMeshupFpeFiles{indexTrial};  
    meshupCapFileName     = outputMeshupCapFiles{indexTrial};  
    outputFpeFileName     = outputFpeFileNames{indexTrial};
    outputCapFileName     = outputCapFileNames{indexTrial};
    outputSegmentationFileName  = outputSegmentationFileNames{indexTrial};
    outputSequenceFileName      = outputMovementSequenceFileNames{indexTrial};    
    
    outputFpeToFootHullDistanceFileName ...
      = outputFpeToFootHullDistanceFileNames{indexTrial}   ;
    outputCapToFootHullDistanceFileName ...
      = outputCapToFootHullDistanceFileNames{indexTrial}   ;
    outputComGPToFootHullDistanceFileName ...
      = outputComGPToFootHullDistanceFileNames{indexTrial} ;
    outputCopToFootHullDistanceFileName ...
      = outputCopToFootHullDistanceFileNames{indexTrial}   ;    
    
    if(exist(outputTrialFolder,'dir') == 0)
      idx = strfind(outputTrialFolder,'/');
      parent = outputTrialFolder(1:1:idx(end-1));
      folder = outputTrialFolder((idx(end-1)+1):1:(end-1));
      mkdir(parent,folder);
    end

    %%
    % Get the trial data
    %%
    headerRows          = 4; 
    textInRowBeforeData = 'ITEM';


    [c3dTime, ...
     c3dMarkers, ... 
     c3dMarkerNames,...
     c3dMarkerUnits,...
     c3dForcePlates, ... 
     c3dForcePlateInfo, ...
     c3dGrf] = getC3DTrialData( inputC3DFolder, ...
                        outputTrialFolder,...
                        c3dFileName,...
                        flag_loadC3DMatFileData, ...
                        flag_useMetersRadiansInC3DData, ...
                        flag_verbose);
                      
            
    [anthroData, anthroColNames] = ...
        getAnthropometryData( inputAnthroFolder,...
                              outputTrialFolder,...
                              anthroFileName,...
                              headerRows,...
                              textInRowBeforeData,... %nanNumberCode,...
                              flag_loadC3DMatFileData, ...
                              flag_verbose);   
                      
    %Anthropometry
    colMass     = getColumnIndex({'MASS';'METRIC';'PROCESSED';'X'},...
                                  headerRows,anthroColNames);
    mass        = anthroData(1,colMass);

    colHeight   = getColumnIndex({'HEIGHT';'METRIC';'PROCESSED';'X'},...
                                headerRows,anthroColNames);
    height      = anthroData(1,colHeight);
      
    if(indexTrial ~= index_Static)
        [ wholeBodyData, ...
          wholeBodyColNames] = ...
            getWholeBodyTrialData( inputV3DFolder,...
                                    outputTrialFolder,...
                                    wholeBodyFileName,...
                                    headerRows,...
                                    textInRowBeforeData, ...%nanNumberCode,...
                                    flag_loadC3DMatFileData,  ...
                                    flag_verbose);

        assert( size(wholeBodyData,1) == size(c3dTime,1),...
                ['Error: number of items in c3d marker data',...
                ' and whole body data should match']);                        

        %%
        % Extract often used data
        %%                        
        %

        %Whole body quantities
        colItem    = getColumnIndex({[];[];[];'ITEM'},headerRows,wholeBodyColNames);

        colComPos = zeros(1,3);
        colComPos(1,1)  = getColumnIndex(...
                      {'LBody_CoM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                      headerRows,wholeBodyColNames);
        colComPos(1,2) = colComPos(1,1)+1;
        colComPos(1,3) = colComPos(1,2)+1;

        colComVel = zeros(1,3);
        colComVel(1,1)  = getColumnIndex(...
                      {'LBody_CoM_vel';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                      headerRows,wholeBodyColNames);
        colComVel(1,2) = colComVel(1,1)+1;
        colComVel(1,3) = colComVel(1,2)+1;

        colHo = zeros(1,3);
        colHo(1,1)  = getColumnIndex(...
                      {'LBody_ANGMOM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                      headerRows,wholeBodyColNames);
        colHo(1,2) = colHo(1,1)+1;
        colHo(1,3) = colHo(1,2)+1;

        colJo       = zeros(1,3);
        colJo(1,1)  = getColumnIndex(...
                    {'LBody_MOMINERT';'LINK_MODEL_BASED';'PROCESSED_MATT';'0'},...
                    headerRows,wholeBodyColNames);
        for i=2:1:9
          colJo(1,i) = colJo(1,i-1)+1;
        end

        %%
        % Check the inertia matrix
        %   Go through the inertia matrix and evaluate its eigen values and
        %   vectors. Do this for 2 reasons
        %     1. To ensure that the inertia matrix is valid: all eigen values > 0
        %     2. To plot/animate the pricipal axis to give the results a check
        %%

        %Only perform this check when processing the data for the first time.
        if(flag_loadC3DMatFileData==0)
          [ JcmEigenValuesDiag, ...
            JcmEigenVectorsRowWise] = ...
            decomposeInertiaMatrices(wholeBodyData(:,colJo(1,:)));  

          %For plotting/debugging purposes compute the whole-body angular velocity
          %and the radius of gyration

          JcmRadiusOfGyration = sqrt(JcmEigenValuesDiag./mass);
          WcmAngularVelocity = zeros(size(JcmEigenValuesDiag,1),3);
          for i=1:1:size(JcmEigenValuesDiag,1)
            JC0 = [wholeBodyData(i,colJo(1,1:3));...
                 wholeBodyData(i,colJo(1,4:6));...
                 wholeBodyData(i,colJo(1,7:9))];
            HC0 = [wholeBodyData(i,colHo(1,1:3))]';

            WcmAngularVelocity(i,:) = zeros(1,3).*NaN;

            if( sum(sum(isnan(JC0))) == 0 && sum(isnan(HC0))==0)
              WcmAngularVelocity(i,:) = (JC0\HC0)';
            end
          end
        end
        %%
        %
        % Evaluate the FPE / Retreive Data
        %
        %%

        %Planes to evaluate the FPE  
        [valMaxRFax, idxMaxRFax] = max(c3dMarkers.('R_FAX')(:,3));
        [valMaxLFax, idxMaxLFax] = max(c3dMarkers.('L_FAX')(:,3));
        idxSeat = 1;
        idxFloor= 2;


        %Numerical tolerances on the solution               
        tol     = 1e-9;
        iterMax = 100;                

        %File to save/laod the file
        flag_fpeVerbose = 0;
        flag_fpeEvaluateDerivatives = 0;
        fpeData = process3DFootPlacementEstimator(...
                    mass,...
                    wholeBodyData,...
                    colComPos,...
                    colComVel,...
                    colJo,...
                    colHo,...
                    gravityVector,...
                    contactPlanes,...
                    tol,...
                    iterMax,...
                    flag_fpeEvaluateDerivatives,...
                    flag_fpeVerbose,...
                    [outputTrialFolder,outputFpeFileName],...
                    flag_loadFpeDataFromFile)  ;

        %%
        %
        % Evaluate the CAP / Retreive Data
        %
        %%
        capData = processCapturePoint(...
                    mass,...
                    wholeBodyData,...
                    colComPos,...
                    colComVel,...
                    gravityVector,...
                    contactPlanes,...
                    [outputTrialFolder,outputCapFileName],...
                    flag_loadCapDataFromFile);

      %%
      %
      % Calculate the distance between key points and the convex hull
      % of the feet
      %
      %%
      if(flag_loadKeyPointDistanceToFootConvexHull==0)
        %Do this both with and excluding toes?        
        for indexFootHull = 1:1:2
        
          c3dFootMarkerRightNames =...
            {'R_FM1','R_FM2','R_FM5','R_FAL','R_FCC','R_TAM'};
          c3dFootMarkerLeftNames =...
            {'L_FM1','L_FM2','L_FM5','L_FAL','L_FCC','L_TAM'};
          
          fpeFileName = outputFpeToFootHullDistanceFileName;
          capFileName = outputCapToFootHullDistanceFileName;
          copFileName = outputCopToFootHullDistanceFileName;
          comFileName = outputComGPToFootHullDistanceFileName;
          
          if(indexFootHull ==2)
            c3dFootMarkerRightNames =...
              {'R_FM1','R_FM5','R_FAL','R_FCC','R_TAM'};
            c3dFootMarkerLeftNames =...
              {'L_FM1','L_FM5','L_FAL','L_FCC','L_TAM'};
            
            fpeFileName = ...
              [outputFpeToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat'];
            capFileName = ...
              [outputCapToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat'];
            comFileName = ...
              [outputComGPToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat'];            
            copFileName = ...
              [outputCopToFootHullDistanceFileName(1:1:(end-4)),'_NoToes.mat'];

          end

          c3dFootMarkerNames = {c3dFootMarkerRightNames{:},...
                            c3dFootMarkerLeftNames{:}};

          fpe2FootConvexHullDist = processDistanceToConvexHull(...
                                      fpeData.r0F0,...
                                      c3dMarkers,...
                                      c3dFootMarkerNames);          
           save([outputTrialFolder,fpeFileName],...
              'fpe2FootConvexHullDist');

          cap2FootConvexHullDist = processDistanceToConvexHull(...
                                      capData.r0F0,...
                                      c3dMarkers,...
                                      c3dFootMarkerNames);          
           save([outputTrialFolder,capFileName],...
              'cap2FootConvexHullDist');

          comgp2FootConvexHullDist = processDistanceToConvexHull(...
                                      wholeBodyData(:,colComPos),...
                                      c3dMarkers,...
                                      c3dFootMarkerNames);          
           save([outputTrialFolder,comFileName],...
              'comgp2FootConvexHullDist');

          cop2FootConvexHullDist = processDistanceToConvexHull(...
                    c3dGrf(index_FeetForcePlate).cop,...
                    c3dMarkers,...
                    c3dFootMarkerNames);

           save([outputTrialFolder,copFileName],...
              'cop2FootConvexHullDist');          
        end
      else       
        load([outputTrialFolder, outputFpeToFootHullDistanceFileName]);
        load([outputTrialFolder, outputCapToFootHullDistanceFileName]);
        load([outputTrialFolder, outputComGPToFootHullDistanceFileName]);
        load([outputTrialFolder, outputCopToFootHullDistanceFileName]);        
      end
      
                  
      
      
      %%
      %
      % Classify the state of the sit-to-stand process
      %
      %%
      
      maxComHeight = max(wholeBodyData(:,colComPos(1,3)));
      minComHeight = min(wholeBodyData(:,colComPos(1,3)));
            
      thresholdStanding = maxComHeight ...
        - thresholdNormDistanceBelowStanding*(maxComHeight-minComHeight);      
 
      c3dFootMarkerRightNames =...
        {'R_FM1','R_FM2','R_FM5','R_FAL','R_FCC','R_TAM'};
      c3dFootMarkerLeftNames =...
        {'L_FM1','L_FM2','L_FM5','L_FAL','L_FCC','L_TAM'};

      [segmentedData,segmentInfo] = segmentDataKMeans(c3dTime,...
                      c3dMarkers,...
                      c3dFootMarkerRightNames,...
                      c3dFootMarkerLeftNames,...
                      wholeBodyData(:,colComPos),...                    
                      c3dGrf(index_ChairForcePlate), ...
                      c3dGrf(index_FeetForcePlate), ...
                      fpeData(1),...
                      thresholdStanding,footContactZMovementTolerance,...
                      [outputTrialFolder,outputSegmentationFileName],...
                      flag_loadSegmentedMotionDataFromFile);
                      
                      
%       [segmentedData,segmentInfo] = segmentData(c3dTime,...
%                       c3dMarkers,...
%                       c3dFootMarkerNames,...
%                       wholeBodyData(:,colComPos),...                    
%                       c3dGrf(index_ChairForcePlate), ...
%                       c3dGrf(index_FeetForcePlate), ...
%                       fpeData(1), ...
%                       thresholdLowForce,...                      
%                       thresholdStaticSittingBalanceDistance,...
%                       thresholdStanding,...
%                       thresholdSitting,...
%                       [outputTrialFolder,outputSegmentationFileName],...
%                       flag_loadSegmentedMotionDataFromFile);
      %thresholdHighForce,...
  
      


      flag_verboseMotionSequence = 1;
      [sitToStandQuietSequence, sitToStandDynamicSequence] = ...
        extractMovementSequence(c3dTime, segmentedData, segmentInfo,...
            thresholdQuietTime, ...
            [outputTrialFolder,outputSequenceFileName],...
            flag_loadMovementSequenceFromFile, ...
            flag_verboseMotionSequence);    

      here=1;
    end



    if(flag_writeRBDLToolChainFiles == 1)        
      if(indexTrial == index_Static)
        %
        % Write ModelFactory File
        %   Note: this is used to generate an RBDL Lua model of the human
        %         subject which can be used to for visualization purposes
        %         and/or to check the results.


        pelvisWidth = ...
          norm(c3dMarkers.('R_IAS')(1,:)-c3dMarkers.('L_IAS')(1,:));

        hipCenterWidth = 2*0.38*pelvisWidth; 
        % Leardini A, Cappozzo A, Catani F, Toksvig-Larsen S, Petitto A, Sforza V, 
        % Cassanelli G, Giannini S. Validation of a functional method for the 
        % estimation of hip joint centre location. Journal of biomechanics. 1999 
        % Jan 1;32(1):99-103.

        shoulderWidth = ...
          0.8*norm(c3dMarkers.('R_SAE')(1,:)-c3dMarkers.('L_SAE')(1,:));

        heelLength = ...
          (norm( 0.5*(c3dMarkers.('R_FAL')(1,:)+c3dMarkers.('R_TAM')(1,:)) ...
                -c3dMarkers.('R_FCC')(1,:)));

        heelHeight = ...
          0.5*(c3dMarkers.('R_FAL')(1,3)+c3dMarkers.('L_FAL')(1,3));

        footWidth = ...
            norm(c3dMarkers.('R_FM1')(1,:)-c3dMarkers.('R_FM5')(1,:));

        shoulderToC7Length = ...
            norm( 0.5.*(c3dMarkers.('R_SAE')(1,3)+c3dMarkers.('L_SAE')(1,3)) ...
                      - c3dMarkers.('CV7')(1,3) );
        humerusHeadRadius = 0.03;
        shoulderToC7Length = shoulderToC7Length + humerusHeadRadius; 

        success = writeModelFactoryModelFile(...
                  [outputTrialFolder, outputModelFactoryAnthropometryFile],...
                  subjectAge,...
                  height,...
                  mass, ...
                  subjectGender1Male0Female, ...
                  pelvisWidth, ...
                  hipCenterWidth, ...
                  shoulderWidth, ...
                  heelLength,...
                  heelHeight,...
                  shoulderToC7Length,...
                  footWidth);
        
                
        success  = writeModelFactoryEnvironmentFile(...
                      [outputTrialFolder, outputModelFactoryEnvironmentFile], ...
                      flag_modelFactoryAddMarkers, ...
                       modelFactoryLuaModel,...
                       modelFactoryDescriptionFile,...                      
                       modelFactoryScalingAlgorithm,...
                       outputModelFactoryAnthropometryFile);
               
        flag_plotModelFactoryModel = 0;
        flag_verboseModelFactoryOutput =0;
        success = ...
          createModel([outputTrialFolder, outputModelFactoryEnvironmentFile],...
                       flag_plotModelFactoryModel,...
                       flag_verboseModelFactoryOutput);
                     
        [success,message,messageId] = ...
          copyfile([outputTrialFolder,modelFactoryLuaModel],...
                   [inputC3DFolders{indexTrial},modelFactoryLuaModel],'f');
                     
      end
      
      if(indexTrial ~= index_Static)
        %
        % Write the ground forces file
        %
        
        nanVal = 0.;
        if(flag_writeC3DDataForMeshup == 1)
          success = writeMeshupGroundForcesFile(...
                      [outputTrialFolder,meshupGrfFileName],...
                      c3dTime,  c3dGrf, nanVal);
        end
                  
        if(flag_writeFpeDataForMeshup == 1)
          success = write3DFootPlacementDataFile(...
                      [outputTrialFolder,meshupFpeFileName],...
                      c3dTime,  fpeData, nanVal);
        end
        if(flag_writeCapDataForMeshup == 1)
          success = writeCapturePointDataFile(...
                      [outputTrialFolder,meshupCapFileName],...
                      c3dTime,  capData, nanVal);
        end        
      end
    end

    %%
    % Vizualization
    %%

    if(flag_visualize == 1 && indexTrial ~= index_Static)

      %%
      % Plot the whole body quantities
      %%  
      fig_wholebody = figure;

      timeInterval = [1:1:100];

      lineColor3 = [1,0,0; 0,1,0; 0,0,1];
      lineColor9 = [1, 0, 0; 0.66,0.33,   0; 0.33, 0.66, 0;...
                    0, 1, 0;    0,0.66,0.33;    0, 0.33, 0.66;...
                    0, 0, 1;    0,  0,    1;    0,    0,    1];

        [ JcmEigenValuesDiag, ...
          JcmEigenVectorsRowWise] = ...
          decomposeInertiaMatrices(wholeBodyData(:,colJo(1,:)));                    
                  
      figure(fig_wholebody);
      subplot(2,2,1);
      for i=1:1:3
        plot( c3dTime(timeInterval,1),...
              wholeBodyData(timeInterval,colComPos(1,i)),...
              'Color',lineColor3(i,:));
        hold on;        
      end
      xlabel('Time (s)');
      ylabel('Distance (m)');
      title('CoM Position');

      subplot(2,2,2);
      for i=1:1:3
        plot( c3dTime(timeInterval,1),...
              wholeBodyData(timeInterval,colComVel(1,i)),...
              'Color',lineColor3(i,:));
        hold on;        
      end
      xlabel('Time (s)');
      ylabel('Distance (m/s)');
      title('CoM Velocity');

      subplot(2,2,3);
      for i=1:1:3
        plot( c3dTime(timeInterval,1),...
              wholeBodyData(timeInterval,colHo(1,i)),...
            'Color',lineColor3(i,:));
        hold on;        
      end
      xlabel('Time (s)');
      ylabel('Angular Momentrum (kg m^2/s)');
      title('Hcm');

      subplot(2,2,4);
      for i=1:1:3
        plot( c3dTime(timeInterval,1),...
              JcmEigenValuesDiag(:,i),...
            'Color',lineColor3(i,:));
        hold on;        
      end
      xlabel('Time (s)');
      ylabel('Inertia (kg m^2)');
      title('Inertia Matrix Eigen Values');

      %%
      % Plot forces
      %%
      fig_forces = figure;
      subplot(2,2,1);
      %plot( c3dTime(timeInterval,1),...



      %%
      % Animate the markers and CoM position
      %%  
      fig_markers3D = figure;

      c3dMarkerFaceColor = [1,1,1].*0.5;
      c3dMarkerSize = 2;
      comMarkerFaceColor = [1,0,0];
      comMarkerSize = 10;

      comVelColor = [0,0,0];
      comAngVelColor = [0,0,1];

      extent = 2000;
      if(flag_useMetersRadiansInC3DData==1)
        extent = extent.*0.001;
      end
      xAxisExtents = [   0,  1  ].*extent;
      yAxisExtents = [-0.5,  0.5].*extent;
      zAxisExtents = [   0,  1  ].*extent;

      numberOfFramesToDraw = 50;
      numberOfFramesToSkip = 5;

      scaleAngularVelocityVector = 1;
      scaleVelocityVector = 1;
      scaleInertiaFrame = 1;

      frameRate     = 25;
      timeToAnimate = 2; 

      nFrames = 1 + numberOfFramesToSkip*frameRate*timeToAnimate;

      %writerObj = VideoWriter('fpeVideo.avi');
      %writerObj.FrameRate = frameRate;
      %open(writerObj);

      for i = 1:numberOfFramesToSkip:(length(c3dTime))

        %Plot the Mocap Markers
        for j=1:1:length(c3dMarkerNames)
          figure(fig_markers3D);

          plot3( c3dMarkers.(c3dMarkerNames{j})(i,1),...
                 c3dMarkers.(c3dMarkerNames{j})(i,2),...
                 c3dMarkers.(c3dMarkerNames{j})(i,3),'o',...
                 'Color',c3dMarkerFaceColor,...
                 'MarkerFaceColor',c3dMarkerFaceColor,...
                 'MarkerSize', 2);
          hold on;      


        end

        %Plot the Com location
        figure(fig_markers3D);

        plot3( wholeBodyData(i,colComPos(1,1)),...
               wholeBodyData(i,colComPos(1,2)),...
               wholeBodyData(i,colComPos(1,3)),...
               'o','Color',comMarkerFaceColor,...
               'MarkerFaceColor',comMarkerFaceColor,...
               'MarkerSize', comMarkerSize);
        hold on;      

        %Plot the Com velocity 
        pointA = wholeBodyData(i,colComPos(1,:))';
        pointB = pointA ...
              + (wholeBodyData(i,colComVel(1,:))').*scaleVelocityVector;
        plot3( [pointA(1,1);pointB(1,1)],...
               [pointA(2,1);pointB(2,1)],...
               [pointA(3,1);pointB(3,1)], ...
               'Color',comVelColor );         
        hold on;
        plot3( [pointB(1,1)],...
               [pointB(2,1)],...
               [pointB(3,1)], ...
               'x','Color',comVelColor,...
               'MarkerFaceColor',[1,1,1],...
               'MarkerSize', comMarkerSize*0.5);         
        hold on;

        %Plot the angular velocity vector    
        pointA = wholeBodyData(i,colComPos(1,:))';
        pointB = pointA ...
              + (WcmAngularVelocity(i,:))'.*(scaleAngularVelocityVector);
        plot3( [pointA(1,1);pointB(1,1)],...
               [pointA(2,1);pointB(2,1)],...
               [pointA(3,1);pointB(3,1)], ...
               'Color',comAngVelColor );         
        hold on;
        plot3( [pointB(1,1)],...
               [pointB(2,1)],...
               [pointB(3,1)], ...
               '+','Color',comAngVelColor,...
               'MarkerFaceColor',[1,1,1],...
               'MarkerSize', comMarkerSize*0.5);         
        hold on;

        %Plot the priciple axis of inertia matrix and scale each by the
        %size of the respective eigen value
        for j=1:1:3
          pointA      = wholeBodyData(i,colComPos(1,:))';
          JcmColIndex = [1,4,7] + (j-1);
          dirAB       = JcmEigenVectorsRowWise(i,JcmColIndex)';
          pointB      = pointA + ...
                   ((dirAB.*JcmRadiusOfGyration(i,j))).*scaleInertiaFrame;

          plot3( [pointA(1,1);pointB(1,1)],...
                 [pointA(2,1);pointB(2,1)],...
                 [pointA(3,1);pointB(3,1)], ...
                 'Color',lineColor3(j,:) );
          hold on;
          plot3( [pointB(1,1)],...
                 [pointB(2,1)],...
                 [pointB(3,1)], ...
                 'o','Color',lineColor3(j,:),...
               'MarkerFaceColor',lineColor3(j,:),...
               'MarkerSize', comMarkerSize*0.25);
          hold on;

        end

        %%
        %Draw the force plates and forces
        %%

        for j=1:1:length(c3dForcePlates)
          plot3( c3dForcePlates(j).corners(1,:),...
                 c3dForcePlates(j).corners(2,:),...
                 c3dForcePlates(j).corners(3,:),...
                 'k','LineWidth',1);
          hold on;
          plot3( c3dForcePlates(j).corners(1,1),...
                 c3dForcePlates(j).corners(2,1),...
                 c3dForcePlates(j).corners(3,1),...
                 'ko','LineWidth',1,...
                 'MarkerSize',10,...
                 'MarkerFaceColor',[1,1,1]);
          hold on;

          fa = (c3dGrf(j).cop(i,:));
          fb = (c3dGrf(j).cop(i,:) ...
             + (c3dGrf(j).force(i,:).*scaleForceToDistance));

          f = [fa;fb];
          plot3(f(:,1),f(:,2),f(:,3),...
                'r','LineWidth',1);
          hold on;
          plot3(f(1,1),f(1,2),f(1,3),...
                'ro','LineWidth',1,...
                'MarkerSize',10,...
                'MarkerFaceColor',[1,1,1]);
          hold on;

          text(fb(1,1),fb(1,2),fb(1,3),...
               sprintf('%1.1f',norm(c3dGrf(j).force(i,:))));

          hold on;

        end

        %%
        % Plot the foot-placement estimator for the seat & floor
        %%
        for k=1:1:length(fpeData)
          plot3( fpeData(k).r0F0(i,1),...
                 fpeData(k).r0F0(i,2),...
                 fpeData(k).r0F0(i,3),...
                 'om','LineWidth',1,...
                 'MarkerSize',10,...
                 'MarkerFaceColor',[1,1,1]);    

          hold on;
%         plot3( fpeData(idxFloor).r0F0(i,1),...
%                fpeData(idxFloor).r0F0(i,2),...
%                fpeData(idxFloor).r0F0(i,3),...
%                'om','LineWidth',1,...
%                'MarkerSize',10,...
%                'MarkerFaceColor',[1,1,1]);    
        %hold on;
        end

        xlim(xAxisExtents);
        ylim(yAxisExtents);
        zlim(zAxisExtents);
        view(26,22);     
        grid on;

        xlabel('X');
        ylabel('Y');
        zlabel('Z');

        %writeVideo(writerObj, getframe(gcf));
        pause(0.05);
        clf(fig_markers3D);
      end
      %close(writerObj);  

    end
  end
end
            
%%
% External folders/libraries: remove
%%          
rmpath(pathToBTK);
rmpath(pathToModelFactory);