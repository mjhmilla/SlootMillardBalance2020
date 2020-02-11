function columnIndex = getColumnIndex(columnData, headerRows, headerData)

columnIndex = -1;

for j=1:1:size(headerData,2)
  flag_allMatch = 1;
  for i=1:1:headerRows  
    if( isempty(columnData{i,1})==0 && isempty(headerData{i,j})==0  ) 
      if(strcmp(columnData{i,1},headerData{i,j}) == 0)
          flag_allMatch = 0;
      end  
    end
  end
  if(flag_allMatch ==1)
    columnIndex = j;
  end
end

