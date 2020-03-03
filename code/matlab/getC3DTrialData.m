function [c3dTime, c3dMarkers,c3dMarkerNames,c3dMarkerUnits,...
          c3dForcePlates, c3dForcePlateInfo, c3dGrf] = ...
    getC3DTrialData(  dataFolderRaw,dataFolderMat,...
                      c3dFileName, flag_loadMatFileData,...
                      flag_MetersRadians, flag_verbose)
   
   
if(flag_loadMatFileData == 0)
  [c3dTime, c3dMarkers,c3dMarkerNames,c3dMarkerUnits,...
   c3dForcePlates, c3dForcePlateInfo, c3dGrf] = ...
    preprocessC3DData( [dataFolderRaw,c3dFileName], flag_MetersRadians,flag_verbose );
  
  idx = strfind(c3dFileName,'.');
  c3dMatFileName = c3dFileName(1:1:(idx-1));    
  save( [dataFolderMat,c3dMatFileName,'.mat'],...
      'c3dTime', 'c3dMarkers','c3dMarkerNames',...
      'c3dForcePlates','c3dForcePlateInfo','c3dGrf','c3dMarkerUnits');
else
  idx = strfind(c3dFileName,'.');
  c3dMatFileName = c3dFileName(1:1:(idx-1));    
  load( [dataFolderMat,c3dMatFileName,'.mat']);
  
  if(strcmp(c3dMarkerUnits.marker,'mm')==1)
    mm2m = 0.001;
    c3dMarkers    = c3dMarkers.*mm2m;
    c3dGrf.cop    = c3dGrf.cop.*mm2m;
    c3dGrf.moment = c3dGrf.moment.*mm2m;
  end
  
end

             