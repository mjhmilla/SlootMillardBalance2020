function success = writeMeshupGroundForcesFile(fileName,c3dTime, c3dGrf, nanVal)

success = 0;

fid = fopen(fileName,'w');     



for i=1:1:(length(c3dTime))
  fprintf(fid,'%1.6f',c3dTime(i,1));
  

  
  for j=1:1:length(c3dGrf)
    if(    sum(isnan(c3dGrf(j).cop(:,1))) == 0 ...
        && sum(isnan(c3dGrf(j).force(:,1))) == 0 ...
        && sum(isnan(c3dGrf(j).moment(:,1))) == 0)      
      
      fprintf(fid,', %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f',...
         c3dGrf(j).cop(i,1),   c3dGrf(j).cop(i,2),   c3dGrf(j).cop(i,3),...
       c3dGrf(j).force(i,1),   c3dGrf(j).force(i,2), c3dGrf(j).force(i,3),...
       c3dGrf(j).moment(i,1), c3dGrf(j).moment(i,2), c3dGrf(j).moment(i,3));
    else
      fprintf(fid,', %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f',...
         nanVal,  nanVal,   nanVal,...
         nanVal,  nanVal,   nanVal,...
         nanVal,  nanVal,   nanVal);
    end
  end
  fprintf(fid,'\n');
    
end
fclose(fid);

success = 1;