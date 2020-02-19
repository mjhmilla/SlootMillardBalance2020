function [c3dTime, c3dMarkers,c3dMarkerNames,...
          c3dForcePlates, c3dForcePlateInfo, c3dGrf] = ...
    preprocessC3DData( c3dFileNameAndPath, flag_verbose )
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

[c3dMarkers, c3dMarkersInfo, c3dMarkersResidual] = btkGetMarkers(c3dH);
c3dMarkerNames = fieldnames(c3dMarkers);

numberOfItems = size(c3dMarkers.(c3dMarkerNames{1}),1);
c3dTime = [0:1:(numberOfItems-1)]'./c3dMarkersInfo.frequency;

%Force plates
inGlobalFrame=1;
fpw = btkGetForcePlatformWrenches(c3dH, inGlobalFrame);
[c3dForcePlates, c3dForcePlateInfo] = btkGetForcePlatforms(c3dH);

numberSkip = size(fpw(1).P,1)/numberOfItems;
timeForcePlate = (size(fpw(1).P,1)-1)/c3dForcePlateInfo(1).frequency;

assert( abs(max(c3dTime) - timeForcePlate) < 1/c3dMarkersInfo.frequency,...
    'Error: Markers and Forceplates have different collection durations.');

c3dGrf(length(c3dForcePlates)) = struct('cop',zeros(numberOfItems,3),...
                                     'forces',zeros(numberOfItems,3),...
                                     'moment',zeros(numberOfItems,3));

%Resolve the force plate moments into a CoP location                                   
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
