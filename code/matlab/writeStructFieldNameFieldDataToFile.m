function fid = writeStructFieldNameFieldDataToFile(fid,dataStruct)

fieldNames = fields(dataStruct);

for z=1:1:length(fieldNames)
  fprintf(fid,'%s\n',fieldNames{z,1});
  for i=1:1:size(dataStruct.(fieldNames{z,1}),1)
    for j=1:1:size(dataStruct.(fieldNames{z,1}),2) 
      fprintf(fid,' %1.15e, ',dataStruct.(fieldNames{z,1})(i,j));
    end
    fprintf(fid,'\n');
  end
end