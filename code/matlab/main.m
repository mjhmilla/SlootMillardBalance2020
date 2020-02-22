clc;
close all;
clear all;


%%
% Subject Meta Data
%%
subjectAge                  = 65; %for now
subjectGender1Male0Female   = 0;
%%
% Processing Flags
%%
flag_loadMatFileData     = 1;
flag_verbose             = 0;
flag_loadFpeDataFromFile = 1;
flag_visualize           = 0;
    numberOfFramesToDraw = 50;
    scaleForceToDistance = 1/1000;

%%
% Input Data 
%%
pathToBTK  = '/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk';

dataFolderC3D              = '../../data/c3d/';
dataFolderV3D              = '../../data/v3d/';

dataFolderMat              = '../../data/mat/';

c3dStaticFileName       = 'static2.c3d';
c3dFileName             = 'sts_0001_Side.c3d';
wholeBodyFileName       = 'DATA.txt';
anthroFileName          = 'Metrics.txt';

headerRows          = 4; 
textInRowBeforeData = 'ITEM';
nanNumberCode       = 1234567890;



modelFactoryFileName = ['../../models/3DHumanDefault_',...
        c3dStaticFileName(1:1:(strfind( c3dStaticFileName,'.')-1)),...
        '.txt'];
%%
% Constants
%%

m2mm         = 1000;
mm2m         = 1/m2mm;
indexPadding = 3; %Numerical differentiation has made the 
                  %indices in velocity-related quantites undefined at the
                  %the ends

%%
% External libraries/folders
%%
addpath(pathToBTK);
                                    
%%
% Get the trial data
%%

[c3dTime, ...
 c3dMarkers, ... 
 c3dMarkerNames,...
 c3dForcePlates, ... 
 c3dForcePlateInfo, ...
 c3dGrf] = getC3DTrialData( dataFolderC3D,dataFolderMat, c3dFileName,...
                            flag_loadMatFileData, flag_verbose);

[ wholeBodyData, ...
  wholeBodyColNames] = ...
    getWholeBodyTrialData( dataFolderV3D,dataFolderMat,wholeBodyFileName,...
                          headerRows,textInRowBeforeData, nanNumberCode,...
                          flag_loadMatFileData, flag_verbose);

assert( size(wholeBodyData,1) == size(c3dTime,1),...
        ['Error: number of items in c3d marker data',...
        ' and whole body data should match']);                        
                        
 [anthroData, anthroColNames] = ...
    getAnthropometryData( dataFolderV3D,dataFolderMat,anthroFileName,...
                          headerRows,textInRowBeforeData, nanNumberCode,...
                          flag_loadMatFileData, flag_verbose);                        
  


                        
%%
% Extract often used data
%%                        
%

%Anthropometry
colMass     = getColumnIndex({'MASS';'METRIC';'PROCESSED';'X'},...
                              headerRows,anthroColNames);
mass        = anthroData(1,colMass);

colHeight   = getColumnIndex({'HEIGHT';'METRIC';'PROCESSED';'X'},...
                            headerRows,anthroColNames);
height      = anthroData(1,colHeight);

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

  WcmAngularVelocity(i,:) = (JC0\HC0)';
end

%%
%
% Evaluate the FPE / Retreive Data
%
%%

gravityVector = [0;0;-9.81];
%Planes to evaluate the FPE  
[valMaxRFax, idxMaxRFax] = max(c3dMarkers.('R_FAX')(:,3));
[valMaxLFax, idxMaxLFax] = max(c3dMarkers.('L_FAX')(:,3));
idxSeat = 1;
idxFloor= 2;

contactPlanes = [ 0,0, mean([valMaxRFax,valMaxLFax])*mm2m; ...
                  0,0,0];  

%Numerical tolerances on the solution               
tol     = 1e-9;
iterMax = 100;                

%File to save/laod the file
idx = strfind(c3dFileName,'.');
fpeMatFileName = [c3dFileName(1:1:(idx-1)),'_FPE.mat'];  

flag_fpeVerbose = 0;
flag_fpeEvaluateDerivatives = 0;
fpeData = process3DFootPlacementEstimator(...
            mass,...
            wholeBodyData,...
            colComPos,colComVel,...
            colJo,colHo,...
            indexPadding,...
            gravityVector,...
            contactPlanes,...
            tol,...
            iterMax,...
            flag_fpeEvaluateDerivatives,...
            flag_fpeVerbose,...
            [dataFolderMat,fpeMatFileName],...
            flag_loadFpeDataFromFile)  ;


%%
%
% Write ModelFactory File
%   Note: this is used to generate an RBDL Lua model of the human
%         subject which can be used to for visualization purposes
%         and/or to check the results.
%%

pelvisWidth = ...
  norm(c3dMarkers.('R_IAS')(1,:)-c3dMarkers.('L_IAS')(1,:))*mm2m;

hipCenterWidth = 2*0.38*pelvisWidth; 
% Leardini A, Cappozzo A, Catani F, Toksvig-Larsen S, Petitto A, Sforza V, 
% Cassanelli G, Giannini S. Validation of a functional method for the 
% estimation of hip joint centre location. Journal of biomechanics. 1999 
% Jan 1;32(1):99-103.

shoulderWidth = ...
  norm(c3dMarkers.('R_SAE')(1,:)-c3dMarkers.('L_SAE')(1,:))*mm2m;

heelLength = ...
  (norm( 0.5*(c3dMarkers.('R_FAL')(1,:)+c3dMarkers.('R_TAM')(1,:)) ...
        -c3dMarkers.('R_FCC')(1,:)))*mm2m;

heelHeight = ...
  0.5*(c3dMarkers.('R_FAL')(1,3)+c3dMarkers.('L_FAL')(1,3))*mm2m;

footWidth = ...
    norm(c3dMarkers.('R_FM1')(1,:)-c3dMarkers.('R_FM5')(1,:))*mm2m;

shoulderToC7Length = ...
    norm( 0.5.*(c3dMarkers.('R_SAE')(1,:)+c3dMarkers.('L_SAE')(1,:)) ...
              - c3dMarkers.('CV7')(1,:) ).*mm2m;

success = writeModelFactoryModelFile(modelFactoryFileName,...
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

%%
% Vizualization
%%
          
if(flag_visualize == 1)

  %%
  % Plot the whole body quantities
  %%  
  fig_wholebody = figure;
  
  timeInterval = [indexPadding:1:100];
  
  lineColor3 = [1,0,0; 0,1,0; 0,0,1];
  lineColor9 = [1, 0, 0; 0.66,0.33,   0; 0.33, 0.66, 0;...
                0, 1, 0;    0,0.66,0.33;    0, 0.33, 0.66;...
                0, 0, 1;    0,  0,    1;    0,    0,    1];
  
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
          JcmEigenValuesDiag(timeInterval,i),...
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
  fig_markers3d = figure;
  
  c3dMarkerFaceColor = [1,1,1].*0.5;
  c3dMarkerSize = 2;
  comMarkerFaceColor = [1,0,0];
  comMarkerSize = 10;
  
  comVelColor = [0,0,0];
  comAngVelColor = [0,0,1];
  
  extent = 2000;
  xAxisExtents = [   0,  1  ].*extent;
  yAxisExtents = [-0.5,  0.5].*extent;
  zAxisExtents = [   0,  1  ].*extent;
  
  numberOfFramesToDraw = 50;
  numberOfFramesToSkip = 5;%floor((length(c3dTime)-indexPadding*2)/numberOfFramesToDraw);
  
  scaleAngularVelocityVector = 1;
  scaleVelocityVector = 1;
  scaleInertiaFrame = 1;
  
  frameRate     = 25;
  timeToAnimate = 2; 
  
  nFrames = indexPadding + numberOfFramesToSkip*frameRate*timeToAnimate;
  
  writerObj = VideoWriter('fpeVideo.avi');
  writerObj.FrameRate = frameRate;
  open(writerObj);
  
  for i = indexPadding:numberOfFramesToSkip:(length(c3dTime)-indexPadding)
       
    %Plot the Mocap Markers
    for j=1:1:length(c3dMarkerNames)
      figure(fig_markers3d);
      
      plot3( c3dMarkers.(c3dMarkerNames{j})(i,1),...
             c3dMarkers.(c3dMarkerNames{j})(i,2),...
             c3dMarkers.(c3dMarkerNames{j})(i,3),'o',...
             'Color',c3dMarkerFaceColor,...
             'MarkerFaceColor',c3dMarkerFaceColor,...
             'MarkerSize', 2);
      hold on;      
            

    end
    
    %Plot the Com location
    figure(fig_markers3d);
      
    plot3( wholeBodyData(i,colComPos(1,1)).*m2mm,...
           wholeBodyData(i,colComPos(1,2)).*m2mm,...
           wholeBodyData(i,colComPos(1,3)).*m2mm,...
           'o','Color',comMarkerFaceColor,...
           'MarkerFaceColor',comMarkerFaceColor,...
           'MarkerSize', comMarkerSize);
    hold on;      

    %Plot the Com velocity 
    pointA = wholeBodyData(i,colComPos(1,:))'.*m2mm;
    pointB = pointA ...
          + (wholeBodyData(i,colComVel(1,:))'.*m2mm).*scaleVelocityVector;
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
    pointA = wholeBodyData(i,colComPos(1,:))'.*m2mm;
    pointB = pointA ...
          + (WcmAngularVelocity(i,:))'.*(m2mm.*scaleAngularVelocityVector);
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
      pointA      = wholeBodyData(i,colComPos(1,:))'.*m2mm;
      JcmColIndex = [1,4,7] + (j-1);
      dirAB       = JcmEigenVectorsRowWise(i,JcmColIndex)';
      pointB      = pointA + ...
               ((dirAB.*JcmRadiusOfGyration(i,j)).*m2mm).*scaleInertiaFrame;
             
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
         + (c3dGrf(j).forces(i,:).*scaleForceToDistance).*m2mm);
      
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
           sprintf('%1.1f',norm(c3dGrf(j).forces(i,:))));
         
      hold on;
      
    end
    
    %%
    % Plot the foot-placement estimator for the seat & floor
    %%
    plot3( fpeData(idxSeat).r0F0(i,1).*m2mm,...
           fpeData(idxSeat).r0F0(i,2).*m2mm,...
           fpeData(idxSeat).r0F0(i,3).*m2mm,...
           'om','LineWidth',1,...
           'MarkerSize',10,...
           'MarkerFaceColor',[1,1,1]);    
    
    hold on;
    plot3( fpeData(idxFloor).r0F0(i,1).*m2mm,...
           fpeData(idxFloor).r0F0(i,2).*m2mm,...
           fpeData(idxFloor).r0F0(i,3).*m2mm,...
           'om','LineWidth',1,...
           'MarkerSize',10,...
           'MarkerFaceColor',[1,1,1]);    
    %hold on;
         
         
    xlim(xAxisExtents);
    ylim(yAxisExtents);
    zlim(zAxisExtents);
    view(26,22);     
    grid on;
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    writeVideo(writerObj, getframe(gcf));
    pause(0.05);
    clf(fig_markers3d);
  end
  close(writerObj);  

end

            
%%
% External folders/libraries: remove
%%          
rmpath(pathToBTK);