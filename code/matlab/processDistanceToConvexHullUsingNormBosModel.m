function output = processDistanceToConvexHullUsingNormBosModel(...
                      r0P0, c3dMarkers, c3dMarkerNames, ...
                      footFrameOffsetLeft,footFrameOffsetRight,...
                      normFootBosLeft,footGeometry) 
                                    
%assert(~isempty(c3dConvexHullMarkerNames));                                    
%assert(size(r0P0,1) == size(c3dMarkers.(c3dConvexHullMarkerNames{1}),1));                                    
output = struct('distance', zeros(size(r0P0,1),1),...
                  'normal', zeros(size(r0P0,1),3),...
                   'point', zeros(size(r0P0,1),3));
                 
r0L0 = zeros(size(normFootBosLeft.shod.convhullNorm,1),3);
r0R0 = zeros(size(normFootBosLeft.shod.convhullNorm,1),3);

footLength = footGeometry.length;
footWidth = footGeometry.width;

flag_debug = 0;
figDebug=[];
if(flag_debug==1)
  figDebug=figure;
end

for i=1:1:size(r0P0,1)
  
  [footFrameLeft, footFrameRight] = ...
          getFootFrames(i,c3dMarkers,c3dMarkerNames,...
                        footFrameOffsetLeft, footFrameOffsetRight);  

    exL = footFrameLeft.E(:,1);
    eyL = footFrameLeft.E(:,2);
    rL  = footFrameLeft.r;
    
    exR = footFrameRight.E(:,1);
    eyR = footFrameRight.E(:,2);
    rR  = footFrameRight.r;                      
  for j=1:1:size(r0L0,1)
    
    r0L0(j,:) = rL ...
      + exL.*( footWidth*normFootBosLeft.shod.convhullNorm(j,1)) ...
      + eyL.*(footLength*normFootBosLeft.shod.convhullNorm(j,2));
    
    r0R0(j,:) = rR ...
      - exR.*( footWidth*normFootBosLeft.shod.convhullNorm(j,1)) ...
      + eyR.*(footLength*normFootBosLeft.shod.convhullNorm(j,2));
    
  end
  
  r0Q0V = [r0L0;r0R0];
  %for j=1:1:length(c3dConvexHullMarkerNames)
  %  r0Q0V(j,:) = c3dMarkers.(c3dConvexHullMarkerNames{j})(i,:);
  %end
  
  idxConvHull = convhull(r0Q0V(:,1:2));
  [dist,normal,point]= calcDistanceToConvexHull(r0P0(i,1:2), ...
                                                r0Q0V(idxConvHull,1:2));
  output.distance(i,1) = dist;
  output.normal(i,:) = [normal,0];
  output.point(i,:)  = [point,0];
  
  if(flag_debug==1)
     c3dFootMarkerNames = {'R_FM1','R_FM2','R_FM5','R_FAL','R_FCC','R_TAM','L_FM1','L_FM2','L_FM5','L_FAL','L_FCC','L_TAM'};
     markerColors = ['r','r','r','r','r','r', ...
                     'b','b','b','b','b','b'];
     for z=1:1:length(c3dFootMarkerNames)
       plot3( c3dMarkers.(c3dFootMarkerNames{z})(i,1),...
              c3dMarkers.(c3dFootMarkerNames{z})(i,2),...
              c3dMarkers.(c3dFootMarkerNames{z})(i,3),...
              ['o',markerColors(1,z)],'MarkerFaceColor',markerColors(1,z));
       hold on;
       mkrLabel = c3dFootMarkerNames{z};
       idx = strfind(mkrLabel,'_');
       mkrLabel(1,idx) = ' ';
       text(c3dMarkers.(c3dFootMarkerNames{z})(i,1),...
              c3dMarkers.(c3dFootMarkerNames{z})(i,2),...
              c3dMarkers.(c3dFootMarkerNames{z})(i,3),...
              mkrLabel);
       hold on;
     end
     
     plot3(r0R0(:,1),r0R0(:,2),r0R0(:,3),'r');
     hold on;     
     plot3(r0L0(:,1),r0L0(:,2),r0L0(:,3),'b');
     hold on;
     plot3(r0Q0V(idxConvHull,1),...
           r0Q0V(idxConvHull,2),...
           r0Q0V(idxConvHull,3),'k');
     hold on;
     
     plot3(r0P0(i,1),r0P0(i,2),0,'om','MarkerFaceColor','m');
     hold on;

     plot3([r0P0(i,1);point(1,1)],[r0P0(i,2);point(1,2)],[0,0],'-k');
     hold on;     
     
     plot3(point(1,1),point(1,2),0,'ok','MarkerFaceColor','k');
     hold on;
     
     plot3(footFrameLeft.r(1,1),footFrameLeft.r(2,1),footFrameLeft.r(3,1),'xb');
     hold on;
     vecColor = ['r','g','b'];
     for a=1:1:size(footFrameLeft.E,2)
       vec = [footFrameLeft.r(:,1)';(footFrameLeft.r(:,1)'+0.1.*footFrameLeft.E(:,a)')];
       plot3(vec(:,1),vec(:,2),vec(:,3),vecColor(1,a));
       hold on;
     end

     plot3(footFrameRight.r(1,1),footFrameRight.r(2,1),footFrameRight.r(3,1),'xr');
     hold on;
     vecColor = ['r','g','b'];
     for a=1:1:size(footFrameRight.E,2)
       vec = [footFrameRight.r(:,1)';(footFrameRight.r(:,1)'+0.1.*footFrameRight.E(:,a)')];
       plot3(vec(:,1),vec(:,2),vec(:,3),vecColor(1,a));
       hold on;
     end
     
     
     
     xlabel('X');
     ylabel('Y');
     zlabel('Z');
     grid on;
     axis equal;
     axis square;
     
     here=1;
  end
  
  
end

if(flag_debug==1)
  close(figDebug);
end
