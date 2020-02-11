clc;
close all;
clear all;

addpath('/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk');

flag_readData = 1;

dataFolderRaw              = '../data/raw/';
dataFolderMat              = '../data/mat/';

c3dFileName             = 'sts_0001_Side.c3d';
wholeBodyFileName       = 'DATA.txt';
anthroFileName          = 'Metrics.txt';

delimiter           = '\t';
textInRowBeforeData = 'ITEM';
nanNumber           = 1234567890;
headerRows          = 4;


%%
% Get the subject's height and weight
%%
anthroData      = [];
anthroColNames  = [];

if(flag_readData == 1)
  [anthroData, anthroColNames] = ...
    getFileAndColumnNames(dataFolderRaw, anthroFileName, delimiter, ...
                          textInRowBeforeData, headerRows, nanNumber);

  i=strfind(anthroFileName,'.');
  fname = [anthroFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'anthroData','anthroColNames');
else
  i=strfind(subjectAnthroFileName,'.');
  fname = [subjectAnthroFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  anthroData = data.anthroData;
  anthroColNames = data.anthroColNames;
end

colHeight = getColumnIndex(...
              {'HEIGHT';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
colMass = getColumnIndex(...
              {'MASS';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
colID = getColumnIndex(...
              {'SUBJECT_ID';'METRIC';'PROCESSED';'X'},...
              headerRows,anthroColNames);
    
height = anthroData(1,colHeight);
mass   = anthroData(1,colMass);
id     = anthroData(1,colID);            

            
            
%%
% Get the whole-body data 
%%

wholeBodyData = [];
wholeBodyColNames = [];

if(flag_readData == 1)
  [wholeBodyData, wholeBodyColNames] = ...
    getFileAndColumnNames(dataFolderRaw, wholeBodyFileName, delimiter, ...
                          textInRowBeforeData, headerRows, nanNumber);

  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  save([dataFolderMat,fname],'wholeBodyData','wholeBodyColNames');
else
  i=strfind(wholeBodyFileName,'.');
  fname = [wholeBodyFileName(1:1:i),'mat'];
  data = load([dataFolderMat,fname]);
  wholeBodyData = data.wholeBodyData;
  wholeBodyColNames = data.wholeBodyColNames;
end
                      
colItem    = getColumnIndex({[];[];[];'ITEM'},headerRows,wholeBodyColNames);
colComPosX  = getColumnIndex(...
              {'LBody_CoM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
              headerRows,wholeBodyColNames);
colComVelX  = getColumnIndex(...
              {'LBody_CoM_vel';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
              headerRows,wholeBodyColNames);
colHoX  = getColumnIndex(...
              {'LBody_ANGMOM';'LINK_MODEL_BASED';'PROCESSED_MATT';'X'},...
              headerRows,wholeBodyColNames);
colJoX  = getColumnIndex(...
            {'LBody_MOMINERT';'LINK_MODEL_BASED';'PROCESSED_MATT';'0'},...
            headerRows,wholeBodyColNames);

            
            
rmpath('/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk');