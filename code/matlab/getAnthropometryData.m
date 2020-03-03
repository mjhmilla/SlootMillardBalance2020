function [anthroData, anthroColNames] = ...
    getAnthropometryData( dataFolderRaw,dataFolderMat,anthroFileName,...
                          headerRows,textInRowBeforeData, ...%nanNumberCode,...
                          flag_loadMatFileData, flag_verbose)
   
   


%%
% Get the subject's height and weight
%%
anthroData      = [];
anthroColNames  = [];




if(flag_loadMatFileData == 0)
  if(flag_verbose)
    disp('Anthropometry Data');
  end

  [anthroData, anthroColNames] = ...
    getFileAndColumnNames(dataFolderRaw, anthroFileName, ...
                          textInRowBeforeData, headerRows, ...%nanNumberCode,...
                          flag_verbose);

  i=strfind(anthroFileName,'.');
  fname = [anthroFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'anthroData','anthroColNames');

else
  if(flag_verbose)
    disp('Anthropometry Data: reading in saved mat structures');
  end
  
  i=strfind(anthroFileName,'.');
  fname = [anthroFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  anthroData = data.anthroData;
  anthroColNames = data.anthroColNames;
  
end

      