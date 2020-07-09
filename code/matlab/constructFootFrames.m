function [footFrameLeft, footFrameRight, ...
          footFrameOffsetLeft, footFrameOffsetRight, ...
          footGeometry] = ...
  constructFootFrames(c3dTime,c3dMarkers,c3dMarkerNames,comPos,comVel)

footFrameLeft   = []; 
footFrameRight  = [];
footGeometry = struct('length',0,'width',0,'lengthMidFoot',0);
comPosZLabels = kmeans(comPos(:,1),3);

hMax = 0;
hMed = 0;
idxH = [];
labelH = 0;

for i=1:1:3
  idxSet = find(comPosZLabels == i);
  hSet = mean(comPos(idxSet,1)); 
  if(hSet > hMax)
    hMax = hSet;
    hMed = median(comPos(idxSet,1));
    idxH = idxSet;
    labelH = i;
  end
end

idxRef = 0;
velRef = Inf;

%Pull out the point that is in the standing cluster that is within
%1cm of the median height of the standing cluster and has the lowest
%velocity
for i=1:1:length(c3dTime)
  
  %Standing cluster
  if(comPosZLabels(i,1) == labelH)
    hErr = abs(comPos(i,1)-hMed);
    if(hErr < 0.01)
      vPoint = sqrt(sum(comVel(i,:).^2));
      
      if(vPoint < velRef)
        idxRef = i;
        velRef = vPoint;
      end      
    end    
  end  
end

%Now go and construct the foot frame treating both feet as if they
%are flat on the ground.

[footFrameOffsetLeft, footFrameOffsetRight]=...
          getFootOffsetFrames(idxRef,c3dMarkers,c3dMarkerNames);
        
[footFrameLeft, footFrameRight] = ...
          getFootFrames(idxRef,c3dMarkers,c3dMarkerNames,...
                        footFrameOffsetLeft, footFrameOffsetRight);       
                             
                             
% Measure the foot length and width

rLenL = footFrameLeft.E'*(...
                   0.5.*(c3dMarkers.('L_FM1')(idxRef,:)' ...
                        +c3dMarkers.('L_FM5')(idxRef,:)')...
                  -0.5.*(c3dMarkers.('L_FAL')(idxRef,:)' ...
                        +c3dMarkers.('L_TAM')(idxRef,:)'));         
midFootLengthL = rLenL(2,1);

rLenR = footFrameRight.E'*(...
         0.5.*(c3dMarkers.('R_FM1')(idxRef,:)' ...
              +c3dMarkers.('R_FM5')(idxRef,:)') ...
        -0.5.*(c3dMarkers.('R_FAL')(idxRef,:)' ...
              +c3dMarkers.('R_TAM')(idxRef,:)'));         

midFootLengthR = rLenR(2,1);

rLenL = footFrameLeft.E'*(c3dMarkers.('L_FM2')(idxRef,:)' ...
                         -c3dMarkers.('L_FCC')(idxRef,:)');         
footLengthL = rLenL(2,1);
rLenR = footFrameRight.E'*(c3dMarkers.('R_FM2')(idxRef,:)' ...
                         -c3dMarkers.('R_FCC')(idxRef,:)');         
footLengthR = rLenR(2,1);

rWL = 0.5.*footFrameLeft.E'*(c3dMarkers.('L_FM1')(idxRef,:)'...
                             -c3dMarkers.('L_FM5')(idxRef,:)') ...
     +0.5.*footFrameLeft.E'*(c3dMarkers.('L_TAM')(idxRef,:)'...
                             -c3dMarkers.('L_FAL')(idxRef,:)');          
footWidthL = rWL(1,1);

rWR = 0.5.*footFrameRight.E'*(c3dMarkers.('R_FM5')(idxRef,:)'...
                             -c3dMarkers.('R_FM1')(idxRef,:)') ...
     +0.5.*footFrameRight.E'*(c3dMarkers.('R_FAL')(idxRef,:)'...
                             -c3dMarkers.('R_TAM')(idxRef,:)');         
footWidthR = rWR(1,1);

footGeometry.length         = 0.5*(footLengthL+footLengthR);
footGeometry.width          = 0.5*(footWidthL+footWidthR);
footGeometry.lengthMidFoot  = 0.5*(midFootLengthL+midFootLengthR);




