function isInSet = isMarkerInSet(markerName, markerPrefix, markerSet)

isInSet = 0;

if(contains(markerName,markerPrefix))  
  for i=1:1:length(markerSet)
    if(contains(markerName,markerSet{i}))
      isInSet=1;
    end
  end  
end