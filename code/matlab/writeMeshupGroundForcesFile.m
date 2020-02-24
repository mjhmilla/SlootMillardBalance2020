function success = writeMeshupGroundForcesFile(fileName,c3dTime,indexPadding, c3dGrf)

success = 0;

fid = fopen(fileName,'w');     



for i=indexPadding:1:(length(c3dTime)-indexPadding)
  fprintf(fid,'%1.6f',c3dTime(i,1));
  
  for j=1:1:length(c3dGrf)
    fprintf(fid,', %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f',...
       c3dGrf(j).cop(i,1),   c3dGrf(j).cop(i,2),   c3dGrf(j).cop(i,3),...
     c3dGrf(j).force(i,1),   c3dGrf(j).force(i,2), c3dGrf(j).force(i,3),...
     c3dGrf(j).moment(i,1), c3dGrf(j).moment(i,2), c3dGrf(j).moment(i,3));
  end
  fprintf(fid,'\n');
    
end
fclose(fid);

success = 1;