function success = write3DFootPlacementDataFile(fpeFileName,c3dTime,indexPadding, fpeData)

success = 0;

for n=1:1:length(fpeData)
  idxEoN = strfind( fpeFileName,'.')-1;
  idxEoN = max(idxEoN);
  fname = [fpeFileName(1:1:idxEoN),'_',num2str(n),...
           fpeFileName((idxEoN+1):1:end)]; 
  fid = fopen(fname,'w');
  
  for i=indexPadding:1:(length(c3dTime)-indexPadding)
    angleN = atan2(-fpeData(n).n(i,1),fpeData(n).n(i,2));
    
    fprintf(fid,'%1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f\n',...
      c3dTime(i,1),...
      fpeData(n).r0F0(i,1),fpeData(n).r0F0(i,2),fpeData(n).r0F0(i,3),...
      fpeData(n).r0G0(i,1),fpeData(n).r0G0(i,2),angleN);
  end
  
  fclose(fid);
end

success = 1;