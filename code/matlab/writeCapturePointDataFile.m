function success = writeCapturePointDataFile(capFileName,c3dTime, capData, nanVal)

success = 0;

for n=1:1:length(capData)
  idxEoN = strfind( capFileName,'.')-1;
  idxEoN = max(idxEoN);
  fname = [capFileName(1:1:idxEoN),'_',num2str(n),...
           capFileName((idxEoN+1):1:end)]; 
  fid = fopen(fname,'w');
  
  for i=1:1:(length(c3dTime))
    
    if(    sum(isnan(capData(n).u(i,:))) == 0 ...
      && sum(isnan(capData(n).r0F0(i,:))) == 0 ...
      && sum(isnan(capData(n).r0G0(i,:))) == 0)
    
      angleN = atan2(capData(n).u(i,2),capData(n).u(i,1));

      fprintf(fid,'%1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f\n',...
        c3dTime(i,1),...
        capData(n).r0F0(i,1),capData(n).r0F0(i,2),capData(n).r0F0(i,3),...
        capData(n).r0G0(i,1),capData(n).r0G0(i,2),angleN);
    
    else
      fprintf(fid,'%1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f\n',...
        c3dTime(i,1),...
        nanVal,nanVal,nanVal,...
        nanVal,nanVal,nanVal);      
    end
  end
  
  fclose(fid);
end

success = 1;