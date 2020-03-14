function output = processDistanceToConvexHull(r0P0, c3dMarkers, ...
                                      c3dConvexHullMarkerNames) 
                                    
assert(~isempty(c3dConvexHullMarkerNames));                                    
assert(size(r0P0,1) == size(c3dMarkers.(c3dConvexHullMarkerNames{1}),1));                                    
output = struct('distance', zeros(size(r0P0,1),1),...
                  'normal', zeros(size(r0P0,1),3),...
                   'point', zeros(size(r0P0,1),3));
                 
r0Q0V = zeros(length(c3dConvexHullMarkerNames),3);

for i=1:1:size(r0P0,1)
  
  for j=1:1:length(c3dConvexHullMarkerNames)
    r0Q0V(j,:) = c3dMarkers.(c3dConvexHullMarkerNames{j})(i,:);
  end
  
  idxConvHull = convhull(r0Q0V(:,1:2));
  [dist,normal,point]= calcDistanceToConvexHull(r0P0(i,1:2), ...
                                                r0Q0V(idxConvHull,1:2));
  output.distance(i,1) = dist;
  output.normal(i,:) = [normal,0];
  output.point(i,:)  = [point,0];
  
  
end
                 