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

%%
%Data Records
%%
dataRecord        = zeros(4,length(subjectsToProcess));

dataStruct = struct('comSpeed',dataRecord,...
                    'angSpeed',dataRecord,...
                    'angZSpeed',dataRecord,...
                    'comCopDist',dataRecord,...
                    'comBosDist',dataRecord);


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
  
  fprintf('COM Velocity (cm/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.comSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.comSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.comSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.comSpeed(4,indexSubject));


  
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
  
              

  fprintf('Angular Speed (deg/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.angSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.angSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.angSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.angSpeed(4,indexSubject));  

  dataStruct.angZSpeed(1,indexSubject) = mean(abs(WcmAngularVelocity(:,3)),'omitnan').*(180/pi);
  dataStruct.angZSpeed(2,indexSubject) = min(abs(WcmAngularVelocity(:,3))).*(180/pi);
  dataStruct.angZSpeed(3,indexSubject) = max(abs(WcmAngularVelocity(:,3))).*(180/pi);
  dataStruct.angZSpeed(4,indexSubject) = std(WcmAngularVelocity(:,3),'omitnan').*(180/pi);
  
  
  fprintf('Angular Speed Z (deg/s)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.angZSpeed(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.angZSpeed(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.angZSpeed(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.angZSpeed(4,indexSubject));  

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
                   
    fprintf('Com-Cop Alignment (cm)\n');
    fprintf('\t%1.3f\tmean\n', dataStruct.comCopDist(1,indexSubject));
    fprintf('\t%1.3f\tmin\n',  dataStruct.comCopDist(2,indexSubject));
    fprintf('\t%1.3f\tmax\n',  dataStruct.comCopDist(3,indexSubject));
    fprintf('\t%1.3f\tstd\n',  dataStruct.comCopDist(4,indexSubject));      
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
                            
  
  
  fprintf('Com-Bos (cm)\n');
  fprintf('\t%1.3f\tmean\n',dataStruct.comBosDist(1,indexSubject));
  fprintf('\t%1.3f\tmin\n', dataStruct.comBosDist(2,indexSubject));
  fprintf('\t%1.3f\tmax\n', dataStruct.comBosDist(3,indexSubject));
  fprintf('\t%1.3f\tstd\n', dataStruct.comBosDist(4,indexSubject));                                         
      
end        

summaryRecord = zeros(4, 5);

dataFields = fieldnames(dataStruct);




for i=1:1:length(dataFields)
  
  summaryRecord(1,i) = mean( dataStruct.(dataFields{i,1})(1,:) );
  summaryRecord(2,i) = min(  dataStruct.(dataFields{i,1})(2,:) );
  summaryRecord(3,i) = max(  dataStruct.(dataFields{i,1})(3,:) );
  summaryRecord(4,i) = mean( dataStruct.(dataFields{i,1})(4,:) );

  disp(dataFields{i,1});
  fprintf('\t%1.2f\t mean\n', summaryRecord(1,i));
  fprintf('\t%1.2f\t min\n' , summaryRecord(2,i));  
  fprintf('\t%1.2f\t max\n' , summaryRecord(3,i));
  fprintf('\t%1.2f\t std\n' , summaryRecord(4,i));  

  
end

rowHeader = {'mean','min','max','std'};

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

