clc;
close all;
clear all;

%List of the subjects to process

flag_NoE06 = 0;

subjectsToProcess = ...
 {'configE01','configE02','configE06','configE07','configE09',...
'configH01','configH02','configH06'}; 
E06Ending = '';

if(flag_NoE06 ==1)
  subjectsToProcess = ...
   {'configE01','configE02','configE07','configE09',...
    'configH01','configH02','configH06'}; 
  E06Ending = '_NoE06';
end

pathToBTK='/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk/';
addpath(pathToBTK);

%%
% Processing Flags
%%
flag_createC3DFilesForRBDL = 0;
c3dPlanarProjection = struct('normal',[0,1,0]);

flag_loadC3DMatFileData = 0;

flag_useMetersRadiansInC3DData = 1;

flag_verbose = 0;

markerVelocityUpperBound = 1;

%%
% Setup the input/output directory structure
%%
inputDirRelative = '../../inputData/QuietStanding_SubAnalysis';
outputDirRelative = '../../outputData/QuietStanding_SubAnalysis';

codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);

gravityVector     = [0;0;-9.81];
gravityVectorDir  = gravityVector./norm(gravityVector);

contactPlanes = [0,0,0]; 

%%
%Data Records
%%
dataRecord        = zeros(6,length(subjectsToProcess));

dataStruct = struct('comSpeed',dataRecord,...
                    'angSpeed',dataRecord,...
                    'angZSpeed',dataRecord,...
                    'comCopDist',dataRecord,...
                    'comBosDist',dataRecord,...
                    'fpeCopInTDist',dataRecord,...
                    'fpeCopInSDist',dataRecord);


summaryRecord = zeros(4, 5);

subjectLabels = cell(size(subjectsToProcess));
%%
%
%%
for indexSubject = 1:1:length(subjectsToProcess)
  
  %%
  % 1. Configure the list of input/output files for each trial
  %%
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});
  cd(codeDir);
  
  %Now we have:
  %   subjectAge                  
  %   subjectGender1Male0Female   
  % 
  %   forcePlateDataRecorded 
  % 
  %   subjectId 
  %   inputFolder  
  %   inputC3DFiles        
  %   inputWholeBodyFiles  
  %   inputAnthroFiles     
  % 
  %   outputFolder = 
  %   
  
  %Set up the input folders

  inputDataFolder = [inputPath,'/',inputFolder];
  outputDataFolder= [outputPath,'/',outputFolder];  
  
  disp(['Processing Quiet Standing: ', subjectId]);
  
  subjectLabels{1,indexSubject} = subjectId;
  
  [c3dTime, ...
  c3dMarkers, ... 
  c3dMarkerNames,...
  c3dMarkerUnits,...
  c3dForcePlates, ... 
  c3dForcePlateInfo, ...
  c3dGrf,...
  c3dGrfDataAvailable] = ...
   getC3DTrialData( inputDataFolder, ...
                    outputDataFolder,...
                    inputC3DFiles,...
                    flag_createC3DFilesForRBDL,...
                    c3dPlanarProjection,...                          
                    flag_loadC3DMatFileData, ...
                    flag_useMetersRadiansInC3DData, ...
                    forcePlateDataRecorded,...
                    flag_verbose);

  indicesOfMarkerJumps = findMarkerJumps(c3dTime,c3dMarkers,...
                                         markerVelocityUpperBound);
  if(isempty(indicesOfMarkerJumps) == 0)
    fprintf('  :%i jumps (%1.3f-%1.3fs is valid)\n',...
      length(indicesOfMarkerJumps), ...
      c3dTime(1,1), ...
      c3dTime(indicesOfMarkerJumps(1,1),1));
  end
      
  idxStart = 1;
  idxEnd = length(c3dTime);
  if(isempty(indicesOfMarkerJumps)==0)
    idxEnd = indicesOfMarkerJumps(1,1);  
  end
  
                  
  headerRows          = 4; 
  textInRowBeforeData = 'ITEM';                  
                  
  [anthroData, anthroColNames] = ...
    getAnthropometryData( inputDataFolder,...
                          outputDataFolder,...
                          inputAnthroFiles,...
                          headerRows,...
                          textInRowBeforeData,... %nanNumberCode,...
                          flag_loadC3DMatFileData, ...
                          flag_verbose);                         

  %Anthropometry
  colMass     = getColumnIndex({'MASS';'METRIC';'PROCESSED';'X'},...
                              headerRows,anthroColNames);
  mass        = anthroData(1,colMass);

  colHeight   = getColumnIndex({'HEIGHT';'METRIC';'PROCESSED';'X'},...
                            headerRows,anthroColNames);
  height      = anthroData(1,colHeight);                              
                              
  [ wholeBodyData, ...
    wholeBodyColNames] = ...
      getWholeBodyTrialData( inputDataFolder,...
                             outputDataFolder,...
                             inputWholeBodyFiles,...
                             headerRows,...
                             textInRowBeforeData, ...%nanNumberCode,...
                             flag_loadC3DMatFileData,  ...
                             flag_verbose);

  assert( size(wholeBodyData,1) == size(c3dTime,1),...
  ['Error: number of items in c3d marker data',...
  ' and whole body data should match']);    

  %Whole body quantities
  colItem    = getColumnIndex({[];[];[];'ITEM'},headerRows,wholeBodyColNames);

  colComPos = zeros(1,3);
  colComPos(1,1)  = getColumnIndex(...
                {'LBody_CoM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                headerRows,wholeBodyColNames);
  colComPos(1,2) = colComPos(1,1)+1;
  colComPos(1,3) = colComPos(1,2)+1;

  colComVel = zeros(1,3);
  colComVel(1,1)  = getColumnIndex(...
                {'LBody_CoM_vel';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                headerRows,wholeBodyColNames);
  colComVel(1,2) = colComVel(1,1)+1;
  colComVel(1,3) = colComVel(1,2)+1;

  colHo = zeros(1,3);
  colHo(1,1)  = getColumnIndex(...
                {'LBody_ANGMOM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
                headerRows,wholeBodyColNames);
  colHo(1,2) = colHo(1,1)+1;
  colHo(1,3) = colHo(1,2)+1;

  colJo       = zeros(1,3);
  colJo(1,1)  = getColumnIndex(...
              {'LBody_MOMINERT';'LINK_MODEL_BASED';'PROCESSED_MATT';'0'},...
              headerRows,wholeBodyColNames);
  for i=2:1:9
    colJo(1,i) = colJo(1,i-1)+1;
  end  
  
%%
%
% Com Velocity
%
%%
  comVelocity = wholeBodyData(idxStart:idxEnd,colComVel(1,:));
  comSpeed = (comVelocity(:,1).^2 ...
             +comVelocity(:,2).^2 ...
             +comVelocity(:,3).^2).^0.5; 

  dataStruct.comSpeed(1,indexSubject) = mean(comSpeed,'omitnan').*100;           
  dataStruct.comSpeed(2,indexSubject) = min(comSpeed).*100;           
  dataStruct.comSpeed(3,indexSubject) = max(comSpeed).*100;           
  dataStruct.comSpeed(4,indexSubject) = std(comSpeed,'omitnan').*100;           
  
  comSpeedSorted = sort(comSpeed.*100);
  idx25 = round(0.25*length(comSpeedSorted));
  idx75 = round(0.75*length(comSpeedSorted));
  dataStruct.comSpeed(5,indexSubject) = comSpeedSorted(idx25,1);
  dataStruct.comSpeed(6,indexSubject) = comSpeedSorted(idx75,1);
  
  fprintf('COM Velocity (cm/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.comSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.comSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.comSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.comSpeed(4,indexSubject));
  fprintf('\t%1.3f\t25p\n', dataStruct.comSpeed(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.comSpeed(6,indexSubject));


  
%%
%
% Average Angular Velocity 
%
%%
  [ JcmEigenValuesDiag, ...
    JcmEigenVectorsRowWise] = ...
    decomposeInertiaMatrices(wholeBodyData(idxStart:idxEnd,colJo(1,:)));  

  %For plotting/debugging purposes compute the whole-body angular velocity
  %and the radius of gyration

  JcmRadiusOfGyration = sqrt(JcmEigenValuesDiag./mass);
  WcmAngularVelocity = zeros(size(JcmEigenValuesDiag,1),3);
  for i=1:1:size(JcmEigenValuesDiag,1)
    JC0 = [wholeBodyData(i,colJo(1,1:3));...
         wholeBodyData(i,colJo(1,4:6));...
         wholeBodyData(i,colJo(1,7:9))];
    HC0 = [wholeBodyData(i,colHo(1,1:3))]';

    WcmAngularVelocity(i,:) = zeros(1,3).*NaN;

    if( sum(sum(isnan(JC0))) == 0 && sum(isnan(HC0))==0)
      WcmAngularVelocity(i,:) = (JC0\HC0)';
    end
  end
  
  angularSpeed = (WcmAngularVelocity(:,1).^2 ...
                +WcmAngularVelocity(:,2).^2 ...
                +WcmAngularVelocity(:,3).^2).^0.5; 

  dataStruct.angSpeed(1,indexSubject) = mean(angularSpeed,'omitnan').*(180/pi);
  dataStruct.angSpeed(2,indexSubject) = min(angularSpeed).*(180/pi);
  dataStruct.angSpeed(3,indexSubject) = max(angularSpeed).*(180/pi);
  dataStruct.angSpeed(4,indexSubject) = std(angularSpeed,'omitnan').*(180/pi);
  
  angSpeedSorted = sort(angularSpeed.*(180/pi));
  idx25 = round(0.25*length(angSpeedSorted));
  idx75 = round(0.75*length(angSpeedSorted));
  dataStruct.angSpeed(5,indexSubject) = angSpeedSorted(idx25,1);
  dataStruct.angSpeed(6,indexSubject) = angSpeedSorted(idx75,1);
              

  fprintf('Angular Speed (deg/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.angSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.angSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.angSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.angSpeed(4,indexSubject));  
  fprintf('\t%1.3f\t25p\n', dataStruct.angSpeed(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.angSpeed(6,indexSubject));  

  
  dataStruct.angZSpeed(1,indexSubject) = mean(abs(WcmAngularVelocity(:,3)),'omitnan').*(180/pi);
  dataStruct.angZSpeed(2,indexSubject) = min(abs(WcmAngularVelocity(:,3))).*(180/pi);
  dataStruct.angZSpeed(3,indexSubject) = max(abs(WcmAngularVelocity(:,3))).*(180/pi);
  dataStruct.angZSpeed(4,indexSubject) = std(WcmAngularVelocity(:,3),'omitnan').*(180/pi);

  angZSpeedSorted = sort(WcmAngularVelocity(:,3).*(180/pi));
  idx25 = round(0.25*length(angZSpeedSorted));
  idx75 = round(0.75*length(angZSpeedSorted));
  dataStruct.angZSpeed(5,indexSubject) = angZSpeedSorted(idx25,1);
  dataStruct.angZSpeed(6,indexSubject) = angZSpeedSorted(idx75,1);
  
  
  
  fprintf('Angular Speed Z (deg/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.angZSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.angZSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.angZSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.angZSpeed(4,indexSubject));  
  fprintf('\t%1.3f\t25p\n', dataStruct.angZSpeed(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.angZSpeed(6,indexSubject));  

  here=1; 
  
%%
%
% COM-COP Alignment 
%
%%  
  idxFP = NaN;
  massErrSmall = Inf;
  flag_GrfIsUsable = 0;
  for i=1:1:length(c3dGrf)
    fmax    = max(c3dGrf(i).force(:,3));
    fpMass  = (fmax/max(abs(gravityVector)));    
    massErr = abs(fpMass-mass)/abs(0.5*(fpMass+mass)); 
    if( massErr < 0.3 )
      idxFP = i;
      flag_GrfIsUsable = 1;
    end
    if(massErr < massErrSmall)
      massErrSmall = massErr;
    end
  end
  
  if(flag_GrfIsUsable == 0)
    fprintf('Force plate data not used\n');
    fprintf('\t%1.3f\n',massErr);
    
  else
    comGroundProj = wholeBodyData(idxStart:idxEnd,colComPos(1,:));
    comGroundProj = comGroundProj ...
                  + comGroundProj.*(gravityVectorDir');
    comCopVec = comGroundProj - c3dGrf(idxFP).cop(idxStart:idxEnd,:);
    comCopDist = sqrt(comCopVec(:,1).^2 ...
                     +comCopVec(:,2).^2 ...
                     +comCopVec(:,3).^2);

  

    dataStruct.comCopDist(1,indexSubject) = mean(comCopDist,'omitnan').*(100);
    dataStruct.comCopDist(2,indexSubject) = min(comCopDist).*(100);
    dataStruct.comCopDist(3,indexSubject) = max(comCopDist).*(100);
    dataStruct.comCopDist(4,indexSubject) = std(comCopDist,'omitnan').*(100);                     

    comCopDistSorted = sort(comCopDist.*(100));
    idx25 = round(0.25*length(comCopDistSorted));
    idx75 = round(0.75*length(comCopDistSorted));
    dataStruct.comCopDist(5,indexSubject) = comCopDistSorted(idx25,1);
    dataStruct.comCopDist(6,indexSubject) = comCopDistSorted(idx75,1);
        
    fprintf('Com-Cop Alignment (cm)\n');
    fprintf('\t%1.3f\tmean\n', dataStruct.comCopDist(1,indexSubject));
    fprintf('\t%1.3f\tmin\n',  dataStruct.comCopDist(2,indexSubject));
    fprintf('\t%1.3f\tmax\n',  dataStruct.comCopDist(3,indexSubject));
    fprintf('\t%1.3f\tstd\n',  dataStruct.comCopDist(4,indexSubject));      
    fprintf('\t%1.3f\t25p\n',  dataStruct.comCopDist(5,indexSubject));
    fprintf('\t%1.3f\t75p\n',  dataStruct.comCopDist(6,indexSubject));      
    
  end

%%
%
% Distance to the BOS
%
%%

  c3dFootMarkerRightNames =...
    {'R_FM1','R_FM2','R_FM5','R_FAL','R_FCC','R_TAM'};
  c3dFootMarkerLeftNames =...
    {'L_FM1','L_FM2','L_FM5','L_FAL','L_FCC','L_TAM'};

  c3dFootMarkerNames = {c3dFootMarkerRightNames{:},...
                        c3dFootMarkerLeftNames{:}};  
                      
  comgp2FootConvexHullDist = processDistanceToConvexHull(...
                                        wholeBodyData(:,colComPos),...
                                        c3dMarkers,...
                                        c3dFootMarkerNames);    

                                      
  dataStruct.comBosDist(1,indexSubject) = mean(comgp2FootConvexHullDist.distance.*(-1),'omitnan').*(100);
  dataStruct.comBosDist(2,indexSubject) = min(comgp2FootConvexHullDist.distance.*(-1)).*(100);
  dataStruct.comBosDist(3,indexSubject) = max(comgp2FootConvexHullDist.distance.*(-1)).*(100);
  dataStruct.comBosDist(4,indexSubject) = std(comgp2FootConvexHullDist.distance.*(-1),'omitnan').*(100);                                       
                            
  comBosDistSorted = sort(comgp2FootConvexHullDist.distance.*(-100));
  idx25 = round(0.25*length(comBosDistSorted));
  idx75 = round(0.75*length(comBosDistSorted));
  dataStruct.comBosDist(5,indexSubject) = comBosDistSorted(idx25,1);
  dataStruct.comBosDist(6,indexSubject) = comBosDistSorted(idx75,1);  
  
  fprintf('Com-Bos (cm)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.comBosDist(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.comBosDist(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.comBosDist(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.comBosDist(4,indexSubject));                                         
  fprintf('\t%1.3f\t25p\n', dataStruct.comBosDist(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.comBosDist(6,indexSubject));                                         

  %%
  %
  % Evaluate the FPE / Retreive Data
  %
  %%

  %Planes to evaluate the FPE  
  [valMaxRFax, idxMaxRFax] = max(c3dMarkers.('R_FAX')(:,3));
  [valMaxLFax, idxMaxLFax] = max(c3dMarkers.('L_FAX')(:,3));
  idxSeat = 1;
  idxFloor= 2;
  
  flag_loadFpeDataFromFile=0;
  omegaSmall = 0.01;

  %Numerical tolerances on the solution               
  tol     = 1e-9;
  iterMax = 100;                

  %File to save/laod the file
  flag_fpeVerbose = 0;
  flag_fpeEvaluateDerivatives = 0;
  fpeData = process3DFootPlacementEstimator(...
              mass,...
              wholeBodyData,...
              colComPos,...
              colComVel,...
              colJo,...
              colHo,...
              gravityVector,...
              contactPlanes,...
              omegaSmall,...
              tol,...
              iterMax,...
              flag_fpeEvaluateDerivatives,...
              flag_fpeVerbose,...
              [outputPath,'/',subjectLabels{1,indexSubject},'_fpe.mat'],...
              flag_loadFpeDataFromFile)  ;
            
            
  
  rFC0t=zeros(size(fpeData.f)).*NaN;
  rFC0s=zeros(size(fpeData.f)).*NaN;
  
  for z=idxStart:1:idxEnd
    rFC0 = c3dGrf(idxFP).cop(z,:)-fpeData.r0F0(z,:);    
    rFC0t(z,1) = sum(rFC0.*fpeData.n(z,:));
    rFC0s(z,1) = sum(rFC0.*fpeData.u(z,:));
  end
            
  dataStruct.fpeCopInT(1,indexSubject) = mean(rFC0t,'omitnan').*(100);
  dataStruct.fpeCopInT(2,indexSubject) = min(rFC0t).*(100);
  dataStruct.fpeCopInT(3,indexSubject) = max(rFC0t).*(100);
  dataStruct.fpeCopInT(4,indexSubject) = std(rFC0t,'omitnan').*(100);                                       

  rFC0tSorted = sort(rFC0t.*(100));
  idx25 = round(0.25*length(rFC0tSorted));
  idx75 = round(0.75*length(rFC0tSorted));
  dataStruct.fpeCopInT(5,indexSubject) = rFC0tSorted(idx25,1);
  dataStruct.fpeCopInT(6,indexSubject) = rFC0tSorted(idx75,1);   
  
  fprintf('FPE-COP in t (cm)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.fpeCopInT(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.fpeCopInT(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.fpeCopInT(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.fpeCopInT(4,indexSubject));   
  fprintf('\t%1.3f\t25p\n', dataStruct.fpeCopInT(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.fpeCopInT(6,indexSubject));   
  
  
  dataStruct.fpeCopInS(1,indexSubject) = mean(rFC0s,'omitnan').*(100);
  dataStruct.fpeCopInS(2,indexSubject) = min(rFC0s).*(100);
  dataStruct.fpeCopInS(3,indexSubject) = max(rFC0s).*(100);
  dataStruct.fpeCopInS(4,indexSubject) = std(rFC0s,'omitnan').*(100);                                       

  rFC0sSorted = sort(rFC0s.*(100));
  idx25 = round(0.25*length(rFC0sSorted));
  idx75 = round(0.75*length(rFC0sSorted));
  dataStruct.fpeCopInS(5,indexSubject) = rFC0sSorted(idx25,1);
  dataStruct.fpeCopInS(6,indexSubject) = rFC0sSorted(idx75,1);   
  
  
  fprintf('FPE-COP in s (cm)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.fpeCopInS(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.fpeCopInS(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.fpeCopInS(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.fpeCopInS(4,indexSubject));    
  fprintf('\t%1.3f\t25p\n', dataStruct.fpeCopInS(5,indexSubject));
  fprintf('\t%1.3f\t75p\n', dataStruct.fpeCopInS(6,indexSubject));    
  
  here=1;
end        

save([outputPath,'/quietStanding',E06Ending,'.mat'],'dataStruct');

summaryRecord = zeros(4, 6);

dataFields = fieldnames(dataStruct);




for i=1:1:length(dataFields)
  
  summaryRecord(1,i) = mean( dataStruct.(dataFields{i,1})(1,:) );
  summaryRecord(2,i) = min(  dataStruct.(dataFields{i,1})(2,:) );
  summaryRecord(3,i) = max(  dataStruct.(dataFields{i,1})(3,:) );
  summaryRecord(4,i) = mean( dataStruct.(dataFields{i,1})(4,:) );
  summaryRecord(5,i) = mean( dataStruct.(dataFields{i,1})(5,:) );
  summaryRecord(6,i) = mean( dataStruct.(dataFields{i,1})(6,:) );  
  
  disp(dataFields{i,1});
  fprintf('\t%1.2f\t mean\n', summaryRecord(1,i));
  fprintf('\t%1.2f\t min\n' , summaryRecord(2,i));  
  fprintf('\t%1.2f\t max\n' , summaryRecord(3,i));
  fprintf('\t%1.2f\t std\n' , summaryRecord(4,i));  
  fprintf('\t%1.2f\t 25p\n' , summaryRecord(5,i));
  fprintf('\t%1.2f\t 75p\n' , summaryRecord(6,i));
  
end

rowHeader = {'mean','min','max','std','25p','75p'};

fid = fopen([outputPath,'/summary',E06Ending,'.csv'],'w');
for i=1:1:length(dataFields)
  fprintf(fid,',%s',dataFields{i,1});
end
fprintf(fid,'\n');

for j=1:1:4
  for i=1:1:length(dataFields)  
    if(i==1)
      fprintf(fid,'%s',rowHeader{1,j});
    end    
    fprintf(fid,',%1.2f',summaryRecord(j,i));
  end
  fprintf(fid,'\n');
end
fclose(fid);



for idxData=1:1:length(dataFields)
  fid=fopen([outputPath,'/',dataFields{idxData,1},E06Ending,'.csv'],'w');
  for i=1:1:length(subjectLabels)
    fprintf(fid,',%s',subjectLabels{1,i});
  end
  fprintf(fid,'\n');
  for i=1:1:size(dataStruct.(dataFields{idxData,1}),1)
    fprintf(fid,'%s',rowHeader{i});    
    for j=1:1:size(dataStruct.(dataFields{idxData,1}),2)
      fprintf(fid,',%1.2f',dataStruct.(dataFields{idxData,1})(i,j));
    end
    fprintf(fid,'\n');
  end  
  fclose(fid);
end

