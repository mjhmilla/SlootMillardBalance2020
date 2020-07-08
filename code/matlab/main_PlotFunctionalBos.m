clc;
close all;
clear all;


subjectsToProcess = {'configPilot1','configPilot2'};

forceLowerBound           = 0.5; %As in half body weight

barefootPadCompressionMax = 0.5*(25*0.5/1000 + 15*0.5/1000);  
shodPadCompressionMax     = 2*barefootPadCompressionMax;  



%%
%Directories
%%
inputDirRelative = '../../inputData/TrueBOS_SubAnalysis';
outputDirRelative = '../../outputData/TrueBOS_SubAnalysis';

pathToBTK='/home/mjhmilla/dev/BTKRoot/BTKCore-install/share/btk-0.3dev/Wrapping/Matlab/btk/';
addpath(pathToBTK);

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
lineWidth  = 0.75;
boxWidth   = 0.33;
panelHeight = 20;
panelWidth  = 10;

numberOfFiguresPerPage        = 2;
numberOfVerticalPlotRows      = 1;
numberOfHorizontalPlotColumns = numberOfFiguresPerPage;    
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
%Output figures
%%

figureBare = figure;
figureShod = figure;
normFootAxis = [-0.8,0.8,-0.4,0.9];

indexRightFoot = 2;
footType = {'leftFoot','rightFoot'};
footwareType = {'shod','bare'};

indexBare = 2;

dataStruct = struct('convhull',[],'convhullNorm',[],'markers',[],'markersNorm',[],'x0',[],'y0',[]);

footStruct = struct('leftFoot',dataStruct,'rightFoot',dataStruct);

for indexFoot=1:1:length(footType)

  footStruct.(footType{indexFoot}) =dataStruct;
end

subjectData(length(subjectsToProcess)) = struct('shod',footStruct, ... 
                                                'bare',footStruct);

groupData = struct('shod',footStruct,...
                   'bare',footStruct);

footData = struct('shod',dataStruct,...
                  'bare',dataStruct);
                 
                                              
for indexSubject=1:1:length(subjectsToProcess)
  subjectData(indexSubject).shod = footStruct;
  subjectData(indexSubject).bare = footStruct;  
end
           
groupData.shod = footStruct;
groupData.bare = footStruct;
footData.shod = dataStruct;
footData.bare = dataStruct;

gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;              
  
subjectColor = [gorgeousGreen(1,:);...
                bellaBlue(1,:)];
footTypeLineType = {'-','--'};
footwareLineWidth = [0.5,0.5];              

flag_markersLabelled = [0,0;0,0];

for indexSubject = 1:1:length(subjectsToProcess)
  cd(inputDirRelative);
   eval(subjectsToProcess{indexSubject});   
  cd(codeDir);
  
  subjectDataFile =[outputPath,'/',outputFolder,'/data',subjectId,'.mat'];
  load(subjectDataFile);
  
  footLength = 0;
  footWidth  = 0;
    
  
  for indexTrial = 1:1:length(trialData)
    if(trialData(indexTrial).isValid==1)
      
      for indexFoot = 1:1:length(footType)
        footware = 'bare';
        if(withShoes(indexTrial,1)==1)
          footware = 'shod';
        end
        
        fprintf('%i,%i\n',indexTrial,indexFoot);
        flag_outlier = 0;
        for z=1:1:length(processedOutliers)
          if(processedOutliers(1,z) == indexTrial)
            flag_outlier=1;
          end
        end
        
        
        if(isempty(trialData(indexTrial).(footType{indexFoot}))==0 ...
            && flag_outlier == 0 )
          if(trialData(indexTrial).(footType{indexFoot}).isValid == 1)
            idxBos = find(trialData(indexTrial).(footType{indexFoot}).f0(:,3)...
                           > forceLowerBound*bodyWeight ...
                      & abs(trialData(indexTrial).(footType{indexFoot}).ea321(:,3)) ...
                           < trialData(indexTrial).angleXUB ...
                      & abs(trialData(indexTrial).(footType{indexFoot}).ea321(:,2)) ...
                          < trialData(indexTrial).angleYUB);          
            if(length(idxBos)>0)

              footLength = trialData(indexTrial).footLength;
              footWidth  = trialData(indexTrial).footWidth;

              %Update the convex hull of the functional BOS
              dataX = [];
              dataY = [];
              dataXNorm = [];
              dataYNorm = [];
              trialX = [];
              trialY = [];
              trialXNorm = [];
              trialYNorm = [];
              trialXConvHull = [];
              trialYConvHull = [];
              trialXNormConvHull = [];
              trialYNormConvHull = [];
              
              trialX = [trialData(indexTrial).(footType{indexFoot}).rFCF(idxBos,1)];
              trialY = [trialData(indexTrial).(footType{indexFoot}).rFCF(idxBos,2)];
              
              idxCH = convhull(trialX,trialY);
                trialXConvHull = trialX(idxCH,:);
                trialYConvHull = trialY(idxCH,:);
        
              trialXNorm = [trialData(indexTrial).(footType{indexFoot}).rFCF(idxBos,1)./footWidth];
              trialYNorm = [trialData(indexTrial).(footType{indexFoot}).rFCF(idxBos,2)./footLength];
                
              idxCh = convhull(trialXNorm,trialYNorm);
                trialXNormConvHull = trialXNorm(idxCH,:);
                trialYNormConvHull = trialYNorm(idxCH,:);
                
              if(isempty(subjectData(indexSubject).(footware).(footType{indexFoot}).convhull))
                dataX = trialX;
                dataY = trialY;

                dataXNorm = trialXNorm;
                dataYNorm = trialYNorm;

              else
                dataX = [subjectData(indexSubject).(footware).(footType{indexFoot}).convhull(:,1);...      
                                        trialX];
                dataY = [subjectData(indexSubject).(footware).(footType{indexFoot}).convhull(:,2);... 
                                        trialY];              

                dataXNorm = [subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,1);...      
                                            trialXNorm];
                dataYNorm = [subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,2);... 
                                            trialYNorm];
              end



              %Update the convex hull of the norm. functional BOS         
              idxCH = convhull(dataX,dataY);
              subjectData(indexSubject).(footware).(footType{indexFoot}).convhull = ...
                [dataX(idxCH,1),dataY(idxCH,1)];

              idxCHNorm = convhull(dataXNorm,dataYNorm);
              subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm = ...
                [dataXNorm(idxCHNorm,1),dataYNorm(idxCHNorm,1)];

              %Append the marker locations            
              markerNames =  fields(trialData(indexTrial).(footType{indexFoot}).markers);

              for indexMarker=1:1:length(markerNames)

                mkrPos = trialData(indexTrial).(footType{indexFoot}).markers.(markerNames{indexMarker})(idxBos,:);
                mkrPosNorm = mkrPos.*[(1/footWidth),(1/footLength),1];

                if(isfield(subjectData(indexSubject).(footware).(footType{indexFoot}).markers,markerNames{indexMarker}))
                  subjectData(indexSubject).(footware).(footType{indexFoot}).markers.(markerNames{indexMarker}) = ...
                    [subjectData(indexSubject).(footware).(footType{indexFoot}).markers.(markerNames{indexMarker});...
                    mkrPos];

                  subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker}) = ...
                    [subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker});...
                    mkrPosNorm];  
                else
                  subjectData(indexSubject).(footware).(footType{indexFoot}).markers.(markerNames{indexMarker}) = mkrPos;
                  subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker})=mkrPosNorm;                
                end

              end
              
              indexFootware=1;
              if(withShoes(indexTrial,1)==1)
                indexFootware=2;
              end
              
              if(withShoes(indexTrial,1)==1)
                figure(figureShod);
              else
                figure(figureBare);
              end
              
              [row,col] = find(subPlotPanelIndex==indexFoot);          
              subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
              subplot('Position',subPlotVec); 

              
              lineColor = subjectColor(indexSubject,:).*0.5 ...
                         +[1,1,1].*0.5;
              lineType  = '-';%footTypeLineType{1,indexFoot};
              lineWidth = footwareLineWidth(1,indexFootware);

              plot( trialXNormConvHull,...
                    trialYNormConvHull,...
                    lineType,'Color',lineColor,'LineWidth',lineWidth);
              hold on;
              
              markerColor = subjectColor(indexSubject,:).*0.5 + [1,1,1].*0.5;
             
              markerNames = fields(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm);

              for indexMarker=1:1:length(markerNames)
                mkrPos = trialData(indexTrial).(footType{indexFoot}).markers.(markerNames{indexMarker})(idxBos,:);
                mkrPosNorm = mkrPos.*[(1/footWidth),(1/footLength),1];                
                
                x0   = mean(mkrPosNorm(:,1));
                xStd = std(mkrPosNorm(:,1));
                y0   = mean(mkrPosNorm(:,2));
                yStd = std(mkrPosNorm(:,2));
                
                mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);
                
                fill( mkrEllipse(:,1),...
                      mkrEllipse(:,2),...
                      markerColor,'EdgeColor',subjectColor(indexSubject,:));
                hold on;
                if(flag_markersLabelled(indexFoot,indexFootware) == 0)
                  markerLabel = markerNames{indexMarker};
                  idx = strfind(markerLabel,'_');
                  markerLabel(1,idx) = ' ';
                  vAlign = 'bottom';
                  hAlign = 'center';
                  dy = 0.05;   
                  if(y0<0)
                    dy = dy*-1.;
                    vAlign = 'top';
                  end
                  text(x0,y0+dy,markerLabel,...
                       'VerticalAlignment',vAlign,...
                       'HorizontalAlignment',hAlign);
                  hold on;
                end
              end
              flag_markersLabelled(indexFoot,indexFootware)=1;             

              %axis equal;
              box off;
              xlabel('Norm. Width');
              ylabel('Norm. Width');
              if(indexFoot==indexRightFoot)
                title(['Right Foot: ', footware ]);
              else
                title(['Left Foot: ', footware]);                
              end
              axis(normFootAxis);
            end
          end
        end
      end
    end    
  end


 %Update the figure
 for indexFoot = 1:1:length(footType)    
    for indexFootware = 1:1:length(footwareType)  
      if(indexFootware==1)        
        figure(figureShod);       
      else
        figure(figureBare); 
      end
      footware = footwareType{indexFootware};

      

      [row,col] = find(subPlotPanelIndex==indexFoot);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
      subplot('Position',subPlotVec);

      lineColor = subjectColor(indexSubject,:);
      lineType  = '-';%footTypeLineType{1,indexFoot};
      lineWidth = 3*footwareLineWidth(1,indexFootware);

      
      
      plot(subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,1),...
        subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,2),...
        lineType,'Color',[1,1,1],'LineWidth',lineWidth*3);
      hold on;
      
      plot(subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,1),...
        subjectData(indexSubject).(footware).(footType{indexFoot}).convhullNorm(:,2),...
        lineType,'Color',lineColor,'LineWidth',lineWidth);
      hold on;
      axis(normFootAxis);
    end    
 end

%Get the average convex hull across all subjects
%  Do this by interpolating using polar coordinates and averaging
%  across the r's



     

%  for indexFoot = 1:1:length(footType)

%     
%     for indexFootware = 1:1:length(footwareType)
%       
%       xSign=1;
%       if(indexRightFoot==indexFoot)
%         xSign=-1;
%       end
%       
%       [row,col] = find(subPlotPanelIndex==indexFootware);          
%           subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
%           subplot('Position',subPlotVec);      
%           
%       lineColor = subjectColor(indexSubject,:);
%       lineType  = footTypeLineType{1,indexFoot};
%       lineWidth = footwareLineWidth(1,indexFootware);
%       
%       plot( xSign.*subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1),...
%             subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
%             lineType,'Color',lineColor,'LineWidth',lineWidth);
%       hold on;
%       
%       markerColor = lineColor;
%       if(indexRightFoot==indexFoot)
%         markerColor = markerColor.*0.5 + [1,1,1].*0.5;
%       end
%       
%       markerNames = fields(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm);
% 
% %       for indexMarker=1:1:length(markerNames)
% %         plot( subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker})(:,1),...
% %                 subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker})(:,2),...
% %                 lineType,'Color',markerColor','LineWidth',lineWidth);
% %         hold on;
% %         if(flag_markersLabelled(1,indexFootware) == 0)
% %           xloc = mean(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker})(:,1));
% %           yloc = mean(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerNames{indexMarker})(:,2));
% %           text(xloc,yloc,markerNames{indexMarker});
% %         end
% %       end
%       flag_markersLabelled(1,indexFootware)=1;
%                  
%     end
%   end
  
end


flag_debugGroup=0;
figDebugGroup=[];
if(flag_debugGroup==1)
  figDebugGroup = figure;
end


for indexFoot = 1:1:length(footType)    
  for indexFootware = 1:1:length(footwareType) 
    angles = [-0.5:0.01:0.5]'.*(2*pi);
    nAngles=length(angles);
    xPts = zeros(nAngles,length(subjectsToProcess));
    yPts = zeros(nAngles,length(subjectsToProcess));
    radius= zeros(nAngles,length(subjectsToProcess));
    
    bosCenterX = zeros(1,length(subjectsToProcess));
    bosCenterY = zeros(1,length(subjectsToProcess));

    markersNorm = subjectData(1).(footware).(footType{indexFoot}).markersNorm;
    markerFields = fields(markersNorm);
    for z=1:1:length(markerFields)
      markersNorm.(markerFields{z}) = zeros(length(subjectsToProcess),4);
    end
    
    for indexSubject = 1:1:length(subjectsToProcess)
      
      for z=1:1:length(markerFields)
        x    = mean(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerFields{z})(:,1));
        xStd = std(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerFields{z})(:,1));        
        y    = mean(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerFields{z})(:,2));
        yStd = std(subjectData(indexSubject).(footware).(footType{indexFoot}).markersNorm.(markerFields{z})(:,2));
        
        markersNorm.(markerFields{z})(indexSubject,:) = [x,y,xStd,yStd];
      end
      
      x0 = 0;
      y0 = 0;
      if(indexFoot == indexRightFoot)
        if(indexFootware == indexBare)
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.0;
        else
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.15;
        end
        groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0 = 0.3;
      else
        if(indexFootware == indexBare)
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.0;
        else
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = -0.15;
        end
        groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0 = 0.3;        
      end
      
      x0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0;
      y0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0;
      
      bosXY = subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm;
      bosXY = bosXY-[x0,y0];
      angleBos = atan2(bosXY(:,2),bosXY(:,1));
      [angleBosSorted,idxSorted]=sort(angleBos);
      bosXYSorted = bosXY(idxSorted,:);
      radiusBosSorted = (bosXYSorted(:,1).^2+bosXYSorted(:,2).^2).^0.5;
      for j=1:1:nAngles                
        pt = calcRayToConvexHullIntersectionPoint(angles(j,1),angleBosSorted,radiusBosSorted);
        xPts(j,indexSubject) = pt(1,1);
        yPts(j,indexSubject) = pt(1,2);           
        radius(j,indexSubject) = sqrt(pt(1,1)*pt(1,1) + pt(1,2)*pt(1,2));
      end
      if(flag_debugGroup==1)
        figure(figDebugGroup)
        clf(figDebugGroup);  

        plot(bosXY(:,1),bosXY(:,2),'k');
        hold on;
        plot(xPts(:,indexSubject),yPts(:,indexSubject),'xm');
        hold on;
        xlabel('X');
        ylabel('Y');
        here=1;
      end
    end

    radiusMean = mean(radius,2);
    
    x0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0;
    y0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0;
    
    xPtsMean = radiusMean.*cos(angles)+x0;
    yPtsMean = radiusMean.*sin(angles)+y0;    
    
    idxCH = convhull(xPtsMean,yPtsMean);
    
    groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm = ...
        [xPtsMean(idxCH,1),yPtsMean(idxCH,1)];
      
    groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm = ...
      markersNorm;
    
    for z=1:1:length(markerFields)
      groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z}) = ...
        mean(markersNorm.(markerFields{z}),1) ;
    end    
    
    markersNorm.(markerFields{z})(indexSubject,:)
    here=1;
  end
end
   
figFootware = figure;
   
for indexFootware = 1:1:length(footwareType)
  
  angles = [-0.5:0.01:0.5]'.*(2*pi);
  nAngles=length(angles);
  xPts = zeros(nAngles,2);
  yPts = zeros(nAngles,2);
  radius= zeros(nAngles,2);

  bosCenterX = 0;
  bosCenterY = 0.2;
  xSign = [1,1];
  xSign(1,indexRightFoot)=-1;
  
  markerSummary = struct('FAL',[],'TAM',[],'FCC',[],'FM1',[],'FM2',[],'FM5',[]);
  
  
  for indexFoot = 1:1:length(footType)    
       
    bosXY      = groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm;
    bosXY(:,1) = bosXY(:,1).*xSign(1,indexFoot);
    
    bosXY    = bosXY-[bosCenterX,bosCenterY];
    angleBos = atan2(bosXY(:,2),bosXY(:,1));
    [angleBosSorted,idxSorted]=sort(angleBos);
    bosXYSorted = bosXY(idxSorted,:);
    radiusBosSorted = (bosXYSorted(:,1).^2+bosXYSorted(:,2).^2).^0.5;
    for j=1:1:nAngles                
      pt = calcRayToConvexHullIntersectionPoint(angles(j,1),angleBosSorted,radiusBosSorted);
      xPts(j,indexFoot) = pt(1,1);
      yPts(j,indexFoot) = pt(1,2);           
      radius(j,indexFoot) = sqrt(pt(1,1)*pt(1,1) + pt(1,2)*pt(1,2));
    end    

    
    [row,col] = find(subPlotPanelIndex==indexFootware);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
      subplot('Position',subPlotVec);

    plot(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1).*xSign(1,indexFoot),...
         groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
         '-','Color',[1,1,1].*0.75);
    hold on;    
    
    markerFields = fields(groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm);
    
    for z=1:1:length(markerFields)
      
      
      
      markerInfo = groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z});
      x0 = markerInfo(1,1).*xSign(1,indexFoot);
      y0 = markerInfo(1,2);
      xStd = markerInfo(1,3);
      yStd = markerInfo(1,4);
      
      mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);
                
      fill( mkrEllipse(:,1),...
            mkrEllipse(:,2),...
            [1,1,1].*0.75,'EdgeColor',[1,1,1].*0.5);
      hold on;
      
    end  
  end
  radiusMean = mean(radius,2);

  xPtsMean = radiusMean.*cos(angles)+bosCenterX;
  yPtsMean = radiusMean.*sin(angles)+bosCenterY;    

  idxCH = convhull(xPtsMean,yPtsMean);  
  
  footData.(footwareType{indexFootware}).convhullNorm = [xPtsMean(idxCH,1),yPtsMean(idxCH,1)];
  
  [row,col] = find(subPlotPanelIndex==indexFootware);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);

    
    
  plot(footData.(footwareType{indexFootware}).convhullNorm(:,1),...
       footData.(footwareType{indexFootware}).convhullNorm(:,2),...
       '-','Color',[0,0,0],'LineWidth',2);
  hold on;
  
  xlabel('Norm. Width');
  ylabel('Norm. Length');
  title(['Norm Bos Template: ',footwareType{indexFootware}]);
    
  box off;
  axis(normFootAxis);
end

    

for indexFoot = 1:1:length(footType)    
  for indexFootware = 1:1:length(footwareType) 
    if(indexFootware==1)        
      figure(figureShod);       
    else
      figure(figureBare); 
    end

    [row,col] = find(subPlotPanelIndex==indexFoot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);

    plot(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1),...
         groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
         '-','Color',[1,1,1],'LineWidth',4);
    hold on;
    plot(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1),...
         groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
         '-','Color',[1,1,1].*0,'LineWidth',2);
    hold on;
    axis(normFootAxis);
  end  
end



save( [outputPath,'/subjectData.mat'] ,'subjectData');
save( [outputPath,'/averageSubjectData.mat'] ,'groupData');



% for indexFootware=1:1:length(footwareType)
%   [row,col] = find(subPlotPanelIndex==indexFootware);          
%       subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
%       subplot('Position',subPlotVec); 
%       
%   xlabel('Norm. Width');
%   ylabel('Norm. Length');
%   title(['Norm. Functional Bos: ',footwareType{indexFootware},' foot']);
% end
figure(figureBare);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBosBarefootLeftRight.pdf']);  
figure(figureShod);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBosShodLeftRight.pdf']);  
figure(figFootware);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBos.pdf']);  
