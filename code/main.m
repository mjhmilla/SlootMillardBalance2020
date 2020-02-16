clc;
close all;
clear all;

%%
% Process control flags
%%
flag_readData   = 0;
flag_verbose    = 0;
flag_visualize  = 1;
    numberOfFramesToDraw = 50;
    scaleForceToDistance = 1/1000;

%%
% Configuration variables
%%
pathToBTK  = '/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk';

dataFolderRaw              = '../data/raw/';
dataFolderMat              = '../data/mat/';

c3dFileName             = 'sts_0001_Side.c3d';
wholeBodyFileName       = 'DATA.txt';
anthroFileName          = 'Metrics.txt';

textInRowBeforeData = 'ITEM';
nanNumber           = 1234567890;
headerRows          = 4;

indexPadding = 3; %Numerical differentiation has made the 
                  %indices in velocity-related quantites undefined at the
                  %the ends
%%
% Constants
%%

m2mm = 1000;


%%
% External libraries/folders
%%
addpath(pathToBTK);

%%
% Get the c3d data
%%

%Markers
[c3dH, c3dByteOrder, c3dStorageFormat] = ...
    btkReadAcquisition([dataFolderRaw,c3dFileName]);

[c3dMarkers, c3dMarkersInfo, c3dMarkersResidual] = btkGetMarkers(c3dH);
c3dMarkerNames = fieldnames(c3dMarkers);

numberOfItems = size(c3dMarkers.(c3dMarkerNames{1}),1);
c3dTime = [0:1:(numberOfItems-1)]'./c3dMarkersInfo.frequency;

%Force plates
inGlobalFrame=1;
fpw = btkGetForcePlatformWrenches(c3dH, inGlobalFrame);
[forceplates, forceplatesInfo] = btkGetForcePlatforms(c3dH);

numberSkip = size(fpw(1).P,1)/numberOfItems;
timeForcePlate = (size(fpw(1).P,1)-1)/forceplatesInfo(1).frequency;

assert( abs(max(c3dTime) - timeForcePlate) < 1/c3dMarkersInfo.frequency,...
    'Error: Markers and Forceplates have different collection durations.');

c3dGrf(length(forceplates)) = struct('cop',zeros(numberOfItems,3),...
                                     'forces',zeros(numberOfItems,3),...
                                     'moment',zeros(numberOfItems,3));
                                   
for i=1:1:length(forceplates)
  origin = forceplates(i).origin';
  for j=1:1:numberOfItems
    k = numberSkip*(j-1)+1;
    p = fpw(i).P(k,:) + origin;
    f = fpw(i).F(k,:);
    m = fpw(i).M(k,:);
    
    fx = f(1,1);
    fy = f(1,2);
    fz = f(1,3);

    mx = m(1,1);
    my = m(1,2);
    mz = m(1,3);    

    dx = -my/fz;
    dy =  mx/fz;
    dz = 0;
    
    xR = [   0,  -dz,  dy;...
            dz,   0, -dx;...
           -dy,  dx,  0];
    dm = (xR * f')';         
    
    c3dGrf(i).forces(j,:) = f;
    c3dGrf(i).cop(j,:)    = p + [dx,dy,dz];
    c3dGrf(i).moment(j,:) = m - dm;
    
  end
end
 
if(flag_verbose==1)
  

  disp('Mocap Marker Information:');
  fprintf('\t%s:\t%d Hz\n','Freq',c3dMarkersInfo.frequency);
  fprintf('\t%s:\t%s\n'   ,'Units',c3dMarkersInfo.units.ALLMARKERS);
  
  for i=1:1:length(c3dMarkerNames)
    fprintf('\t%d. %s\n',i,c3dMarkerNames{i,1});
  end
end
btkCloseAcquisition(c3dH);


%%
% Get the subject's height and weight
%%
anthroData      = [];
anthroColNames  = [];

if(flag_readData == 1)
  if(flag_verbose)
    disp('Anthropometry Data');
  end
  
  [anthroData, anthroColNames] = ...
    getFileAndColumnNames(dataFolderRaw, anthroFileName, ...
                          textInRowBeforeData, headerRows, nanNumber,...
                          flag_verbose);

  i=strfind(anthroFileName,'.');
  fname = [anthroFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'anthroData','anthroColNames');

else
  if(flag_verbose)
    disp('Anthropometry Data: reading in saved mat structures');
  end
  
  i=strfind(anthroFileName,'.');
  fname = [anthroFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  anthroData = data.anthroData;
  anthroColNames = data.anthroColNames;
  
end

colHeight = getColumnIndex(...
              {'HEIGHT';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
colMass = getColumnIndex(...
              {'MASS';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
colID = getColumnIndex(...
              {'SUBJECT_ID';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
    
height = anthroData(1,colHeight);
mass   = anthroData(1,colMass);
id     = anthroData(1,colID);            

                       
%%
% Get the whole-body data 
%%
wholeBodyData = [];
wholeBodyColNames = [];

if(flag_readData == 1)

  if(flag_verbose)
    disp('Wholebody Data');
  end  
  [wholeBodyData, wholeBodyColNames] = ...
    getFileAndColumnNames(dataFolderRaw, wholeBodyFileName,...
                          textInRowBeforeData, headerRows, nanNumber,...
                          flag_verbose);

  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'wholeBodyData','wholeBodyColNames');
  
 
else
  
  if(flag_verbose)
    disp('Wholebody Data: reading in saved mat structures');
  end 
  
  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  wholeBodyData = data.wholeBodyData;
  wholeBodyColNames = data.wholeBodyColNames;
    
end
                      
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

%Go through the inertia matrix and evaluate its eigen values and
%vectors. Do this for 2 reasons
% 1. To ensure that the inertia matrix is valid: all eigen values > 0
% 2. To plot/animate the pricipal axis to give the results a check

assert( size(wholeBodyData,1) == size(c3dTime,1),...
        ['Error: number of items in c3d marker data',...
        ' and whole body data should match']);

JcmEigenValuesDiag   = zeros(length(c3dTime),3);
JcmRadiusOfGyration  = zeros(length(c3dTime),3);
JcmEigenVectorsRowWise  = zeros(length(c3dTime),9);
WcmAngularVelocity   = zeros(length(c3dTime),3);

for i=1:1:length(c3dTime)

  jrow = wholeBodyData(i,colJo(1,:));
  Jcmi = [jrow(1,1),jrow(1,2),jrow(1,3);...
          jrow(1,4),jrow(1,5),jrow(1,6);...
          jrow(1,7),jrow(1,8),jrow(1,9)];  
  Hcmi = wholeBodyData(i,colHo(1,:))';
  
  Wcmi = Jcmi\Hcmi;
        
  [v,d] = eig(Jcmi);
  

  
  
  JcmEigenValuesDiag(i,:)     = [d(1,1),d(2,2),d(3,3)];
  JcmRadiusOfGyration(i,:)    = [sqrt(d(1,1)/mass),sqrt(d(2,2)/mass),sqrt(d(3,3)/mass)];
  JcmEigenVectorsRowWise(i,:) = [v(1,:),v(2,:),v(3,:)];
  WcmAngularVelocity(i,:)     = Wcmi';

  assert( d(1,1) > 0 && d(2,2) > 0 && d(3,3) > 0,...  
          ['Error: whole body inertia matrix has a negative eigen value,'...
          'which is physically impossible. ']);
end

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
  xAxisExtents = [-0.5,  0.5].*extent;
  yAxisExtents = [-0.5,  0.5].*extent;
  zAxisExtents = [   0,  1  ].*extent;
  
  numberOfFramesToDraw = 50;
  numberOfFramesToSkip = floor((length(c3dTime)-indexPadding*2)/numberOfFramesToDraw);
  
  scaleAngularVelocityVector = 1;
  scaleVelocityVector = 1;
  scaleInertiaFrame = 1;
  
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
    
    for j=1:1:length(forceplates)
      plot3( forceplates(j).corners(1,:),...
             forceplates(j).corners(2,:),...
             forceplates(j).corners(3,:),...
             'k','LineWidth',1);
      hold on;
      plot3( forceplates(j).corners(1,1),...
             forceplates(j).corners(2,1),...
             forceplates(j).corners(3,1),...
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
    
    xlim(xAxisExtents);
    ylim(yAxisExtents);
    zlim(zAxisExtents);
    view(26,22);     
    grid on;
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    pause(0.05);
    clf(fig_markers3d);
  end
  
end

            
%%
% External folders/libraries: remove
%%          
rmpath(pathToBTK);