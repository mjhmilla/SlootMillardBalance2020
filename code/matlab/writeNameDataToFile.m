function fid = writeNameDataToFile(fid,varName,varData)

fprintf(fid,'%s\n',varName);
for i=1:1:size(varData,1)
  for j=1:1:size(varData,2) 
    fprintf(fid,' %1.15e, ',varData(i,j));
  end
  fprintf(fid,'\n');
end