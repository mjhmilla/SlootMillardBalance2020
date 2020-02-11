function [data, headerData] = ...
    getFileAndColumnNames(folderName,dataFileName,delimiter,...
                          keywordInRowBeforeData, headerRows, nanNumber) 

fid = fopen([folderName,dataFileName]);

tline = fgetl(fid);

headerText = cell(headerRows,1);


flag_rowBeforeData  = 0;
flag_dataRow        = 0;
dataColumns         = 0;
dataFormat          = '';
dataRow = [];
data    = [];
while(ischar(tline))
  
  %Store header text: Keep track of the last n lines for files with multiple header lines
  if(headerRows > 1 && flag_rowBeforeData == 0)
    for i=1:1:(headerRows-1)
      headerText{i} = headerText{i+1};
    end
    headerText{headerRows} = tline;
  else
    headerText{1} = tline;
  end

  %Header processing: Look for the pre-data key word
  if(flag_rowBeforeData == 0)
    loc = strfind(tline,keywordInRowBeforeData);
    if(isempty(loc) == 0)
      flag_rowBeforeData = 1;
    end
  end
  
  %Process Header
  if(flag_rowBeforeData == 1 && flag_dataRow == 0)
    dcount = zeros(headerRows,1);
    for i=1:1:headerRows
      entries = textscan(headerText{i},'%s');
      dcount(i,1) = length(entries{1});
    end
    headerData = cell(headerRows, max(dcount));
    for i=1:1:headerRows
      entries = textscan(headerText{i},'%s');
      idxStart = max(dcount)-length(entries{1})+1;
      headerData(i,idxStart:1:max(dcount)) = entries{1}';
    end
    dataColumns = max(dcount);
    dataRow = zeros(1,dataColumns);
    dataFormat = '%d';
    for i=2:1:dataColumns
      dataFormat = [dataFormat,'\t%f'];
    end
  end
  
  %Process Data
  if(flag_dataRow ==1)
    
    lineData = textscan(tline,dataFormat);
    for i=1:1:dataColumns
      if( isnan(lineData{i})==1)
        dataRow(1,i) = nanNumber;
      else
        dataRow(1,i) = lineData{i};
      end
    end
    
    if(isempty(data) ==1)
      data = dataRow;
    else
      data = [data; dataRow];
    end
        
  end
  
  %Get the next line of data
  tline  = fgetl(fid);
  
  if(flag_rowBeforeData == 1)
    flag_dataRow = 1;
  end
    
    
end


fclose(fid);

