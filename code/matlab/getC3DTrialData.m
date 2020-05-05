function [c3dTime, c3dMarkers,c3dMarkerNames,c3dMarkerUnits,...
          c3dForcePlates, c3dForcePlateInfo, c3dGrf,c3dGrfDataAvailable] = ...
    getC3DTrialData(  dataFolderRaw,dataFolderMat,...
                      c3dFileName, c3dPlanarSettings, ...
                      flag_loadMatFileData, flag_MetersRadians, ...
                      flag_grfDataRecorded, flag_verbose)
   
   
c3dGrfDataAvailable = NaN;
c3dPlanar = [];

if(flag_loadMatFileData == 0)
  [c3dTime, c3dMarkers,c3dMarkerNames,c3dMarkerUnits,...
   c3dForcePlates, c3dForcePlateInfo, c3dGrf,c3dGrfDataAvailable,...
   c3dPlanar] = ...
    preprocessC3DData( [dataFolderRaw,c3dFileName], c3dPlanarSettings,...
                       flag_MetersRadians,flag_grfDataRecorded,...
                       flag_verbose );
  
  idx = strfind(c3dFileName,'.');
  c3dMatFileName = c3dFileName(1:1:(idx-1));    
  save( [dataFolderMat,c3dMatFileName,'.mat'],...
      'c3dTime', 'c3dMarkers','c3dMarkerNames','c3dMarkerUnits',...
      'c3dForcePlates','c3dForcePlateInfo','c3dGrf',...
      'c3dGrfDataAvailable');
  if(isempty(c3dPlanar)==0)
    if(c3dPlanarSettings.write == 1)
      c3dPlanarFileName = [c3dFileName(1:1:(idx-1)),'_2D.c3d'];      
      btkWriteAcquisition(c3dPlanar,[dataFolderMat,c3dPlanarFileName]);
    end
  end
else
  idx = strfind(c3dFileName,'.');
  c3dMatFileName = c3dFileName(1:1:(idx-1));    
  load( [dataFolderMat,c3dMatFileName,'.mat']);
  
  %Make sure this variable gets overwritten.
  assert(isnan(c3dGrfDataAvailable)==0)
  
  if(strcmp(c3dMarkerUnits.marker,'mm')==1)
    mm2m = 0.001;
    c3dMarkers    = c3dMarkers.*mm2m;
    c3dGrf.cop    = c3dGrf.cop.*mm2m;
    c3dGrf.moment = c3dGrf.moment.*mm2m;
  end
  
end

             