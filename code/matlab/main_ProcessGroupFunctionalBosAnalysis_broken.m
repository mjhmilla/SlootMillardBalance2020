clc;
close all;
clear all;


subjectsToProcess = { 'configPilot1',...
                      'configPilot2',...
                      'configP01',...
                      'configP02',...
                      'configP03',...
                      'configP04',...
                      'configP05',...
                      'configP06',...
                      'configP07',...
                      'configP08',...                      
                      'configP10'};
%'configP09',... Left foot looks very offset

forceLowerBound           = 0.5; %As in half body weight

barefootPadCompressionMax = (2.0*max(25*0.5/1000,15*0.5/1000)); %A high limit to detect clear lift off
shodPadCompressionMax     = barefootPadCompressionMax+5/1000;  

minimumNumberOfPoints = 150;


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
lineWidth   = 0.75;
boxWidth    = 0.33;
panelWidth  = 7;
panelHeight = panelWidth*2*(16/14.5);


numberOfFiguresPerPage        = 3;
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

figureBareFoot = figure;
figureRunningShoe = figure;
figureHikingShoe = figure;
figFootware = figure;

normFootAxis = [-0.7,0.7,-0.4,0.9];
normXTicks = [ normFootAxis(1,1):0.1:normFootAxis(1,2)];
normYTicks = [ normFootAxis(1,3):0.1:normFootAxis(1,4)];

indexRightFoot = 2;
footType = {'leftFoot','rightFoot'};
footwareType = {'run','hike','bare'};

indexBare = 3;


dataStruct = struct('convhull',[],'convhullNorm',[],'markers',[],'markersNorm',[],'x0',[],'y0',[],'isValid',0);
dataStructNorm=struct('convhullNorm',[],'markersNorm',[],...
                      'x',[],'y',[],'xErr',[],'yErr',[],...
                      'xHeel',[],'yHeel',[],'xToe',[],'yToe',[],...
                      'xLat',[],'yLat',[],'xMed',[],'yMed',[]);

footStruct = struct('leftFoot',dataStruct,'rightFoot',dataStruct);

for indexFoot=1:1:length(footType)

  footStruct.(footType{indexFoot}) =dataStruct;
end

subjectData(length(subjectsToProcess)) = struct('run',footStruct, ... 
                                                'hike',footStruct, ... 
                                                'bare',footStruct);

groupData = struct('run',footStruct,...
                   'hike',footStruct,...
                   'bare',footStruct);

footData = struct('run',dataStructNorm,...
                  'hike',dataStructNorm,...
                  'bare',dataStructNorm);
                 
                                              
for indexSubject=1:1:length(subjectsToProcess)
  subjectData(indexSubject).run = footStruct;
  subjectData(indexSubject).hike = footStruct;  
  subjectData(indexSubject).bare = footStruct;  
end

groupData.run = footStruct;
groupData.hike = footStruct;
groupData.bare = footStruct;

footData.run = dataStructNorm;
footData.hike = dataStructNorm;
footData.bare = dataStructNorm;

gorgeousGreen       = [102 204 0]./255; 
bellaBlue           = [51 153 255]./255; 
ostentatiousOrange  = [255 128 0]./255;              
red = [1,0,0];

grey = [1,1,1].*0.5;
n1 = 0.0;
n2 = 0.25;
n3 = 0.50;
n4 = 0.75;

subjectColor = [(gorgeousGreen.*(1-n1)+bellaBlue.*(n1));...
                (gorgeousGreen.*(1-n2)+bellaBlue.*(n2));...
                (gorgeousGreen.*(1-n3)+bellaBlue.*(n3));...
                (gorgeousGreen.*(1-n4)+bellaBlue.*(n4));...
                (bellaBlue.*(1-n1)+red.*(n1));...
                (bellaBlue.*(1-n2)+red.*(n2));...
                (bellaBlue.*(1-n3)+red.*(n3));...
                (bellaBlue.*(1-n4)+red.*(n4));...
                (red.*(1-n1)+ostentatiousOrange.*n1);...
                (red.*(1-n2)+ostentatiousOrange.*n2);...
                (red.*(1-n3)+ostentatiousOrange.*n3);...
                (red.*(1-n4)+ostentatiousOrange.*n4)];
footTypeLineType = {'-','-','--'};
footwareLineWidth = [0.5,0.5,0.5];              

flag_markersLabelled = [0,0;0,0;0,0];

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
      footware = '';
      for indexFoot = 1:1:length(footType)
        
        if(indexTrial > length(withShoes))
          here=1;
        end
        
        switch withShoes(indexTrial,1)
          case 0
            footware = 'bare';
          case 1
            footware = 'run';            
          case 2
            footware = 'hike';
          otherwise
            assert(0)
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
            if(length(idxBos)>minimumNumberOfPoints)

              footLength = trialData(indexTrial).footLength;
              footWidth  = trialData(indexTrial).footWidth;
              fprintf('%1.3f by %1.3f : %i Subject, %i Trial\n',footLength,footWidth,indexSubject,indexTrial);
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
              

              
              switch withShoes(indexTrial,1)
                case 0
                  indexFootware=1;
                  figure(figureBareFoot);
                case 1
                  indexFootware=2;
                  figure(figureRunningShoe);
                case 2
                  indexFootware=2;
                  figure(figureHikingShoe);
                otherwise
                  assert(0);
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
 footware = '';
 for indexFoot = 1:1:length(footType)    
    for indexFootware = 1:1:length(footwareType)
      

      switch indexFootware
        case 1
          figure(figureRunningShoe); 
        case 2
          figure(figureHikingShoe); 
        case 3
          figure(figureBareFoot);
        otherwise
          assert(0);
      end
      
      footware = footwareType{indexFootware};

      

      [row,col] = find(subPlotPanelIndex==indexFoot);          
      subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
      subplot('Position',subPlotVec);

      lineColor = subjectColor(indexSubject,:);
      lineType  = '-';%footTypeLineType{1,indexFoot};
      lineWidth = 3*footwareLineWidth(1,indexFootware);


      
      if(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).isValid==1)
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
 end

  
end


flag_debugGroup=0;
figDebugGroup=[];
if(flag_debugGroup==1)
  figDebugGroup = figure;
end

footware = '';
for indexFoot = 1:1:length(footType)    
  for indexFootware = 1:1:length(footwareType) 
    angles = [-0.5:0.01:0.5]'.*(2*pi);
    nAngles=length(angles);
    xPts = zeros(nAngles,length(subjectsToProcess));
    yPts = zeros(nAngles,length(subjectsToProcess));
    radius= zeros(nAngles,length(subjectsToProcess));
    
    bosCenterX = zeros(1,length(subjectsToProcess));
    bosCenterY = zeros(1,length(subjectsToProcess));

    markersNorm = subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).markersNorm;
    markerFields = fields(markersNorm);
    for z=1:1:length(markerFields)
      markersNorm.(markerFields{z}) = zeros(length(subjectsToProcess),4);
    end
    
    indicesOfValidSubjectData = [];
    
    for indexSubject = 1:1:length(subjectsToProcess)
      if(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).isValid==1)
        for z=1:1:length(markerFields)
          x    = mean(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z})(:,1));
          xStd = std(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z})(:,1));        
          y    = mean(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z})(:,2));
          yStd = std(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{z})(:,2));

          markersNorm.(markerFields{z})(indexSubject,:) = [x,y,xStd,yStd];
        end

        x0 = 0;
        y0 = 0;
        if(indexFoot == indexRightFoot)
          switch indexFootware
            case 1
              %run
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.15;
            case 2
              %hike
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.15;
            case 3
              %bare
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.0;
            otherwise
              assert(0);
          end
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0 = 0.3;

        else
          switch indexFootware
            case 1
              %run
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = -0.15;
            case 2
              %hike
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = -0.15;
            case 3
              %bare
              groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0 = 0.0;
            otherwise
              assert(0);
          end
          groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0 = 0.3;

        end

        x0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).x0;
        y0 = groupData.(footwareType{indexFootware}).(footType{indexFoot}).y0;

        bosXY = subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm;
        bosXY = bosXY-[x0,y0];

        flag_bosValid = 0;
        minX = min(bosXY(:,1));
        maxX = max(bosXY(:,1));
        minY = min(bosXY(:,2));
        maxY = max(bosXY(:,2));

        if(minX < 0 && minY < 0 && maxX > 0 && maxY > 0)
          flag_bosValid = 1;
          indicesOfValidSubjectData = [indicesOfValidSubjectData;...
                                       indexSubject];
        else
          fprintf('Rejected: %i\tSubject\t%i\tTrial\n',indexSubject,indexTrial);
        end

        subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).isValid = ...
          flag_bosValid;

        angleBos = atan2(bosXY(:,2),bosXY(:,1));
        [angleBosSorted,idxSorted]=sort(angleBos);
        bosXYSorted = bosXY(idxSorted,:);
        radiusBosSorted = (bosXYSorted(:,1).^2+bosXYSorted(:,2).^2).^0.5;
        for j=1:1:nAngles      
          if(flag_bosValid==1)
            pt = calcRayToConvexHullIntersectionPoint(angles(j,1),angleBosSorted,radiusBosSorted);
            xPts(j,indexSubject) = pt(1,1);
            yPts(j,indexSubject) = pt(1,2);           
            radius(j,indexSubject) = sqrt(pt(1,1)*pt(1,1) + pt(1,2)*pt(1,2));
          else
            xPts(j,indexSubject) = NaN;
            yPts(j,indexSubject) = NaN;           
            radius(j,indexSubject) = NaN;          
          end
        end
        if(flag_debugGroup==1)
          figure(figDebugGroup)
          clf(figDebugGroup);  

          if(flag_boxValid==1)
            plot(bosXY(:,1),bosXY(:,2),'k');
            hold on;
            plot(xPts(:,indexSubject),yPts(:,indexSubject),'xm');
            hold on;
            xlabel('X');
            ylabel('Y');
            here=1;
          end
        end
      end

    end

    if(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).isValid==1)
      radiusMean = mean(radius(:,indicesOfValidSubjectData),2);

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

    end
    here=1;
  end
end
   
figure(figFootware);   

xSign = [1,1];
xSign(1,indexRightFoot)=-1;

for indexSubject = 1:1:length(subjectsToProcess)
  for indexFoot = 1:1:length(footType)
    for indexFootware = 1:1:length(footwareType)
          [row,col] = find(subPlotPanelIndex==indexFootware);          
            subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
            subplot('Position',subPlotVec);
            
          if(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).isValid==1)
            plot(subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1).*xSign(1,indexFoot),...
                 subjectData(indexSubject).(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
                 '-','Color',[1,1,1].*0.75);
            hold on;    
          end
    end
  end
end

for indexFootware = 1:1:length(footwareType)
  
  angles = [-0.5:0.01:0.5]'.*(2*pi);
    
  nAngles=length(angles);
  xPts = zeros(nAngles,2);
  yPts = zeros(nAngles,2);
  radius= zeros(nAngles,2);

  bosCenterX = 0;
  bosCenterY = 0.2;
  markerSummary = struct('FAL',[],'TAM',[],'FCC',[],'FM1',[],'FM2',[],'FM5',[]);
  markerSummaryFields = fields(markerSummary);
  
  isRadiusValid = zeros(length(footType),1);
  
  for indexFoot = 1:1:length(footType)    
       
    bosXY      = groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm;
    
    if(isempty(bosXY)==0)
      isRadiusValid(indexFoot,1) = 1;
      
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
    end
    
    
    

    
       
    if(isfield(groupData.(footwareType{indexFootware}).(footType{indexFoot}),'markersNorm'))
      if(isempty(groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm)==0)
        markerFields = fields(groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm);
        for idxMarker=1:1:length(markerFields)           

          markerInfo = groupData.(footwareType{indexFootware}).(footType{indexFoot}).markersNorm.(markerFields{idxMarker});
          x0 = markerInfo(1,1).*xSign(1,indexFoot);
          y0 = markerInfo(1,2);
          xStd = markerInfo(1,3);
          yStd = markerInfo(1,4);

          %mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);

          %fill( mkrEllipse(:,1),...
          %      mkrEllipse(:,2),...
          %      [1,1,1].*0.75,'EdgeColor',[1,1,1].*0.5);
          %hold on;

          found = 0;
          for z=1:1:length(markerSummaryFields)
            if(contains(markerFields{idxMarker},markerSummaryFields{z}))
              assert(found==0);
              found = 1;
              if(isempty(markerSummary.(markerSummaryFields{z}))==1)
                markerSummary.(markerSummaryFields{z}) = [x0,y0,xStd,yStd];
              else
                markerSummary.(markerSummaryFields{z}) = ...
                  [markerSummary.(markerSummaryFields{z})(:,:);...
                    x0,y0,xStd,yStd];            
              end
            end
          end

        end 
      end
    end
  end
  
  countValidRadius=0;
  radiusMean=zeros(size(radius,1),1);
  for z=1:1:length(isRadiusValid)
    if(isRadiusValid(z,1)==1)
      radiusMean = radiusMean + radius(:,z);
      countValidRadius=countValidRadius+1;
    end
  end
  radiusMean = radiusMean./countValidRadius;
  
 

  xPtsMean = radiusMean.*cos(angles)+bosCenterX;
  yPtsMean = radiusMean.*sin(angles)+bosCenterY;  

  footData.(footwareType{indexFootware}).x = zeros(size(radius,1),size(radius,2));
  footData.(footwareType{indexFootware}).y = zeros(size(radius,1),size(radius,2));  
  
  footData.(footwareType{indexFootware}).xHeel = zeros(1,size(radius,2));
  footData.(footwareType{indexFootware}).yHeel = zeros(1,size(radius,2));
  
  footData.(footwareType{indexFootware}).xToe = zeros(1,size(radius,2));
  footData.(footwareType{indexFootware}).yToe = zeros(1,size(radius,2));

  footData.(footwareType{indexFootware}).xMed = zeros(1,size(radius,2));
  footData.(footwareType{indexFootware}).yMed = zeros(1,size(radius,2));

  footData.(footwareType{indexFootware}).xLat = zeros(1,size(radius,2));
  footData.(footwareType{indexFootware}).yLat = zeros(1,size(radius,2));
  
  for z=1:1:size(radius,2)
    
    if(isRadiusValid(z,1)==1)
        
      footData.(footwareType{indexFootware}).x(:,z) = ...
        radius(:,z).*cos(angles)+bosCenterX;    
      footData.(footwareType{indexFootware}).y(:,z) = ...
        radius(:,z).*sin(angles)+bosCenterY;    

      footData.(footwareType{indexFootware}).xErr(:,z) = ...
        footData.(footwareType{indexFootware}).x(:,z) ...
        - xPtsMean;    
      footData.(footwareType{indexFootware}).yErr(:,z) = ...
        footData.(footwareType{indexFootware}).y(:,z) ...
        - yPtsMean;    

      %'xHeel',[],'yHeel',[],'xToe',[],'yToe',[],...
      %                      'xLat',[],'yLat',[],'xMed',[],'yMed',[]    

      [val,idxHeel] = min(footData.(footwareType{indexFootware}).y(:,z));

      footData.(footwareType{indexFootware}).xHeel(1,z) = ...
        footData.(footwareType{indexFootware}).x(idxHeel,z);
      footData.(footwareType{indexFootware}).yHeel(1,z) = ...
        footData.(footwareType{indexFootware}).y(idxHeel,z);    

      [val,idxToe] = max(footData.(footwareType{indexFootware}).y(:,z));

      footData.(footwareType{indexFootware}).xToe(1,z) = ...
        footData.(footwareType{indexFootware}).x(idxToe,z);
      footData.(footwareType{indexFootware}).yToe(1,z) = ...
        footData.(footwareType{indexFootware}).y(idxToe,z);    

      [val,idxLat] = min(footData.(footwareType{indexFootware}).x(:,z));

      footData.(footwareType{indexFootware}).xLat(1,z) = ...
        footData.(footwareType{indexFootware}).x(idxLat,z);
      footData.(footwareType{indexFootware}).yLat(1,z) = ...
        footData.(footwareType{indexFootware}).y(idxLat,z);    

      [val,idxMed] = max(footData.(footwareType{indexFootware}).x(:,z));

      footData.(footwareType{indexFootware}).xMed(1,z) = ...
        footData.(footwareType{indexFootware}).x(idxMed,z);
      footData.(footwareType{indexFootware}).yMed(1,z) = ...
        footData.(footwareType{indexFootware}).y(idxMed,z);   
    end
    
  end
  
  disp(['Norm. Foot Template Heel Stats: ',footwareType{indexFootware}]);
  fprintf('  %1.3f: err mean\n',mean(footData.(footwareType{indexFootware}).yHeel));
  fprintf('  %1.3f: err min\n',  min(footData.(footwareType{indexFootware}).yHeel));
  fprintf('  %1.3f: err max\n',  max(footData.(footwareType{indexFootware}).yHeel));  
  fprintf('All: \n');
  disp(footData.(footwareType{indexFootware}).yHeel);
  
  if(sum(isRadiusValid) > 0)
    idxCH = convhull(xPtsMean,yPtsMean);    
    footData.(footwareType{indexFootware}).convhullNorm = [xPtsMean(idxCH,1),yPtsMean(idxCH,1)];
  end
  
  [row,col] = find(subPlotPanelIndex==indexFootware);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);

    footData.(footwareType{indexFootware}).markersNorm=markerSummary;
    
  for z=1:1:length(markerSummaryFields)
    
    if(isempty(markerSummary.(markerSummaryFields{z}))==0)
      
      x0    = mean(markerSummary.(markerSummaryFields{z})(:,1));
      y0    = mean(markerSummary.(markerSummaryFields{z})(:,2));    
      xStd  = mean(markerSummary.(markerSummaryFields{z})(:,3));
      yStd  = mean(markerSummary.(markerSummaryFields{z})(:,4));    

      footData.(footwareType{indexFootware}).markersNorm.(markerSummaryFields{z}) = [x0,y0,xStd,yStd];

      mkrEllipse = getEllipse([x0,y0],[xStd,yStd],20);

      fill( mkrEllipse(:,1),...
            mkrEllipse(:,2),...
            [1,1,1].*0.75,'EdgeColor',[1,1,1].*0.5);
      hold on;

      vAlign = 'bottom';
      hAlign = 'center';
      dy = 0.05;   
      if(y0<0)
        dy = dy*-1.;
        vAlign = 'top';
      end

      text(x0,y0+dy,markerSummaryFields{z},...
            'VerticalAlignment',vAlign,...
            'HorizontalAlignment',hAlign);
      hold on;
    end
    
  end
    
  if(isempty(footData.(footwareType{indexFootware}).convhullNorm)==0)
    plot(footData.(footwareType{indexFootware}).convhullNorm(:,1),...
         footData.(footwareType{indexFootware}).convhullNorm(:,2),...
         '-','Color',[0,0,0],'LineWidth',2);
    hold on;
  end
  
  plot(0,0,'ok','MarkerFaceColor','k');
  hold on;
  plot([0;0.2],[0;0],'-k','MarkerFaceColor','k');
  hold on;
  plot([0;0],[0;0.2],'-k','MarkerFaceColor','k');
  hold on;
  
  if(isempty(footData.(footwareType{indexFootware}).convhullNorm)==0)
    xMinCH = min(footData.(footwareType{indexFootware}).convhullNorm(:,1));
    xMaxCH = max(footData.(footwareType{indexFootware}).convhullNorm(:,1));
    xZeroCH = 0;
    xTicksCH = sort([xMinCH,xZeroCH,xMaxCH]);
    yMinCH = min(footData.(footwareType{indexFootware}).convhullNorm(:,2));
    yMaxCH = max(footData.(footwareType{indexFootware}).convhullNorm(:,2));
    yZeroCH = 0;
    yTicksCH = sort([yMinCH,yZeroCH,yMaxCH]);
  end
  
  xlabel('Norm. Width');
  ylabel('Norm. Length');
  title(['Functional BOS : ',footwareType{indexFootware}]);
    
  box off;
  if(isempty(footwareType{indexFootware}.convhullNorm)==0)  
    axis(normFootAxis);
    xticks(round(xTicksCH,2));
    yticks(round(yTicksCH,2));
  end
  grid on;
end

save([outputPath,'/normBosModel.mat'],'footData');
    

for indexFoot = 1:1:length(footType)    
  for indexFootware = 1:1:length(footwareType) 
    if(indexFootware==1)        
      figure(figureRunningShoe);       
    else
      figure(figureBareFoot); 
    end

    [row,col] = find(subPlotPanelIndex==indexFoot);          
    subPlotVec = reshape(subPlotPanel(row,col,:),1,4);         
    subplot('Position',subPlotVec);

    if(isempty(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm)==0)
      plot(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1),...
           groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
           '-','Color',[1,1,1],'LineWidth',4);
      hold on;
      plot(groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,1),...
           groupData.(footwareType{indexFootware}).(footType{indexFoot}).convhullNorm(:,2),...
           '-','Color',[1,1,1].*0,'LineWidth',2);
      hold on;
      axis(normFootAxis);
      xticks(normXTicks);
      yticks(normYTicks);
    end
    grid on;    
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
figure(figureBareFoot);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBosBarefootLeftRight.pdf']);  
figure(figureRunningShoe);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBosShodLeftRight.pdf']);  
figure(figFootware);
configPlotExporter;
print('-dpdf',[outputPath,'/normFunctionalBos.pdf']);  
