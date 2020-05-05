function [c3dTime, c3dMarkers,c3dMarkerNames, c3dMarkerUnits,...
          c3dForcePlates, c3dForcePlateInfo, c3dGrf, c3dGrfDataAvailable,c3dPlanar] = ...
    preprocessC3DData( c3dFileNameAndPath, c3dPlanarSettings, ...
                       flag_MetersRadians, flag_grfDataRecorded,...
                       flag_verbose )
%%
% This function will read in a c3d file, resolve force plate readings to
% forces and moments acting at the center of pressure, and give everything 
%a common sense of time 
%
% @param c3dFileNameAndPath : full path to the c3d file
%
% @return 
%   c3dTime: an n-by-1 vector of time samples for both marker and grf data
%   c3dMarkers: a struct with fields of marker names containing n-by-3 data
%   c3dMarkerNames : a cell array of the marker names
%   c3dGrf: an array of structs (one per force plate) that contain fields
%%

%Markers
[c3dH, c3dByteOrder, c3dStorageFormat] = ...
    btkReadAcquisition(c3dFileNameAndPath);

c3dGrfDataAvailable = flag_grfDataRecorded;

unitMarker = btkGetPointsUnit(c3dH,'MARKER');
unitAngle  = btkGetPointsUnit(c3dH,'ANGLE');
unitMoment = btkGetPointsUnit(c3dH,'MOMENT');

scaleDistance = 1.;
scaleAngles   = 1.;
scaleMoments  = 1.;

if(flag_MetersRadians == 1)
  if(strcmp(unitMarker,'mm')==1)
    scaleDistance = 0.001;
  end
  if(strcmp(unitMoment,'Nmm')==1)
    scaleMoments = 0.001;
  end
  if(strcmp(unitAngle,'deg')==1)
    scaleAngles = 180/pi;
  end
    
  unitMarker ='m';
  unitAngle  = 'rad';
  unitMoment = 'Nm';
end

if(flag_verbose == 1)  
  disp('C3D Units');
  disp([unitMarker, ' : Points Distance Unit']);
  disp([unitAngle,  ' : Points Angle Unit']);
  disp([unitMoment, ' : Points Moment Unit']);
end

c3dMarkerUnits = struct('marker',unitMarker,...
                        'angle' ,unitAngle,...
                        'moment',unitMoment);

[c3dMarkers, c3dMarkersInfo, c3dMarkersResidual] = btkGetMarkers(c3dH);
c3dMarkerNames = fieldnames(c3dMarkers);

for i=1:1:length(c3dMarkerNames)
  c3dMarkers.(c3dMarkerNames{i}) = ...
    c3dMarkers.(c3dMarkerNames{i}).*scaleDistance;
end


numberOfItems = size(c3dMarkers.(c3dMarkerNames{1}),1);
c3dTime = [0:1:(numberOfItems-1)]'./c3dMarkersInfo.frequency;


%Force plates

inGlobalFrame=1;
fpw = btkGetGroundReactionWrenches(c3dH, inGlobalFrame);


%Check to see if there is ground force data to process
if(length(fpw) == 0)
  c3dGrfDataAvailable = 0;
end

flag_dataIsAllIdentical = 1;
if(c3dGrfDataAvailable == 1)
    for i=1:1:length(fpw)
      minF = min(min(fpw(i).F));
      maxF = max(max(fpw(i).F));

      if(abs(maxF-minF) > sqrt(eps))
        flag_dataIsAllIdentical = 0;
      end
    end
end

if(flag_dataIsAllIdentical==1)
  c3dGrfDataAvailable = 0;
end

if(c3dGrfDataAvailable==0)
  c3dForcePlates    = []; 
  c3dForcePlateInfo = []; 
  c3dGrf            = [];
end

%If there is ground force data to process, process it.
if(c3dGrfDataAvailable ==1)
  [c3dForcePlates, c3dForcePlateInfo] = btkGetForcePlatforms(c3dH);


  numberSkip = size(fpw(1).P,1)/numberOfItems;
  timeForcePlate = (size(fpw(1).P,1)-1)/c3dForcePlateInfo(1).frequency;

  assert( abs(max(c3dTime) - timeForcePlate) < 1/c3dMarkersInfo.frequency,...
      'Error: Markers and Forceplates have different collection durations.');

  c3dGrf(length(c3dForcePlates)) = struct('cop',zeros(numberOfItems,3),...
                                       'force',zeros(numberOfItems,3),...
                                       'moment',zeros(numberOfItems,3));


  for i=1:1:length(c3dForcePlates)
    origin = c3dForcePlates(i).origin';
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

      c3dGrf(i).force(j,:) = f;        
      c3dGrf(i).cop(j,:)    = (p + [dx,dy,dz]).*scaleDistance;
      c3dGrf(i).moment(j,:) = (m - dm).*scaleDistance;

    end
  end

  for i=1:1:length(c3dForcePlates)
    c3dForcePlates(i).corners = c3dForcePlates(i).corners.*scaleDistance;
    c3dForcePlates(i).origin  = c3dForcePlates(i).origin.*scaleDistance;  
  end
end

if(c3dPlanarSettings.generate == 1)
  pn = btkGetPointNumber(c3dH);       %number of points
  fn = btkGetPointFrameNumber(c3dH);  %number of frames;
  an = btkGetAnalogNumber(c3dH); %number of analog channels
  %r  = btkGetAnalogSampleNumberPerFrame(c3dH); %number of samples per video frame.
  r = 1;
  c3dPlanar = btkNewAcquisition(pn,fn,an,r);
  btkSetFrequency(c3dPlanar, c3dMarkersInfo.frequency)
  assert( abs( norm(c3dPlanarSettings.normal)-1) < eps*10);
  
  projDir = c3dPlanarSettings.normal;
  
  for i=1:1:length(c3dMarkerNames)
    [points,pointInfo] = btkSetPointLabel(c3dPlanar,i,c3dMarkerNames{i});
    
    [values,residuals,info] = btkGetPoint(c3dH,i);
    
    for j=1:1:fn
      values(j,:) = values(j,:) - sum(values(j,:).*projDir).*projDir;
    end        
    [points,pointInfo]=btkSetPoint(c3dPlanar,i,values);    
  end
  
  %if(c3dGrfDataAvailable ==1)
  %  for i=1:1:size(c3dGrf)      
  %  end
  %end
    
else
  c3dPlanar = [];
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
