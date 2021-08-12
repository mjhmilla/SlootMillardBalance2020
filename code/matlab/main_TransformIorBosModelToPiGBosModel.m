clc;
close all;
clear all;


%%
%Global variables
%%
footType = {'leftFoot','rightFoot'};
footwareType = {'shod','bare'};

%%
%Directories
%%
inputDirRelative = '../../inputData/TrueBOS_SubAnalysis';
outputDirRelative = '../../outputData/TrueBOS_SubAnalysis';


%%
%Force plate COP offset correction
%%
tmp = load(['../../inputData/COP_SubAnalysis/',...
                    'errorForcePlateAtIndex1.mat']);
fpAtIndex1ErrorStruct = tmp.errorStruct;                  

tmp = load(['../../inputData/COP_SubAnalysis/',...
                    'errorForcePlateAtIndex2.mat']);
fpAtIndex2ErrorStruct = tmp.errorStruct; 

%%
%Plot settings
%%
normFootAxis = [-0.65,0.65,-0.4,0.9];
normXTicks = [ normFootAxis(1,1):0.1:normFootAxis(1,2)];
normYTicks = [ normFootAxis(1,3):0.1:normFootAxis(1,4)];

lineWidth   = 0.75;
boxWidth    = 0.33;
panelWidth  = 7;
panelHeight = panelWidth*2;


numberOfFiguresPerPage        = 4;
numberOfVerticalPlotRows      = 2;
numberOfHorizontalPlotColumns = numberOfFiguresPerPage/numberOfVerticalPlotRows;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;           
pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;           
plotHeight  = panelHeight;
plotWidth   = panelWidth ;

plotConfigGeneric;



%%
% Set up the input/output directory structure
%%
codeDir = pwd;
  cd(inputDirRelative);
  inputPath = pwd;
  cd(codeDir);
  cd(outputDirRelative);
  outputPath = pwd;  
cd(codeDir);

%%
% Load normalized Bos model
%%

tmp=load([outputPath,'/normBosModel.mat']);
footData=tmp.footData;


footDataPig = transformIorBosToPigBos(footData,footwareType);

figBos = figure;
figBos = plotFootBosModel(footData,figBos,subPlotPanel,subPlotPanelIndex,...
                          normFootAxis,normXTicks,normYTicks,0);
figBos = plotFootBosModel(footDataPig,figBos,subPlotPanel,subPlotPanelIndex,...
                          normFootAxis,normXTicks,normYTicks,2);

figure(figBos);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBosIorPig.pdf']);  


save([outputPath,'/normBosModelPiG.mat'],'footDataPig');



