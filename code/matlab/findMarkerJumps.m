function indicesOfMarkerJumps = findMarkerJumps(timeVec, markers, ...
                                                velThresh)

indicesOfMarkerJumps = [];                                              
nameOfFields = fieldnames(markers);
velThreshSq = velThresh*velThresh;

for i=1:1:length(nameOfFields)
  dx = calcCentralDifference(timeVec,markers.(nameOfFields{i})(:,1));
  dy = calcCentralDifference(timeVec,markers.(nameOfFields{i})(:,2));
  dz = calcCentralDifference(timeVec,markers.(nameOfFields{i})(:,3));
  
  vel = (dx.*dx + dy.*dy + dz.*dz);
  idx = find(vel >= velThreshSq);
  if(isempty(idx)==0)
    indicesOfMarkerJumps = [indicesOfMarkerJumps;idx];
  end
  
end