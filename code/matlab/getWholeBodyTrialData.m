function [wholeBodyData, wholeBodyColNames] = ...
    getWholeBodyTrialData(  dataFolderRaw,dataFolderMat,wholeBodyFileName,...
                            headerRows,textInRowBeforeData,...% nanNumberCode,...
                            flag_loadMatFileData, flag_verbose)
   
                  
%%
% Get the whole-body data 
%%
wholeBodyData = [];
wholeBodyColNames = [];



if(flag_loadMatFileData == 0)

  if(flag_verbose)
    disp('Wholebody Data');
  end  
  

  [wholeBodyData, wholeBodyColNames] = ...
    getFileAndColumnNames(dataFolderRaw, wholeBodyFileName,...
                          textInRowBeforeData, headerRows, ... %nanNumberCode,...
                          flag_verbose);

  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'wholeBodyData','wholeBodyColNames');
  
 
else
  
  if(flag_verbose)
    disp('Wholebody Data: reading in saved mat structures');
  end 
  
  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  wholeBodyData = data.wholeBodyData;
  wholeBodyColNames = data.wholeBodyColNames;
    
end                