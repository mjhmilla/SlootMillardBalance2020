function output = calcTimeToMaximum(movementSequence, c3dTime,data)

output = zeros(size(c3dTime));
 
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
 
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;
    
    [valSeqMax, idxSeqMax] = max(data(idx0:1:idx2,1));
    
    idxMax = idxSeqMax+idx0-1;
    
    output(idx0:1:idx2,1) = c3dTime(idx0:1:idx2,1)-c3dTime(idxMax,1);
    
  end
  
end