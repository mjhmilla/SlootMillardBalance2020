clc;
close all;
clear all;
 
filesToProcess = {'test_Lframe_0002'};

flag_plotScene = 0;
flag_printError = 0;

 
data =load([filesToProcess{1},'.mat']);
data =data.(filesToProcess{1});


assert(data.Trajectories.Labeled.Count==0)

figTest=[];
dTime = 1;
if(flag_plotScene==1)
  figTest=figure;
end

%Plot configuration
lineWidth  = 0.75;
boxWidth   = 0.33;
panelWidth  = 25;
panelHeight = (400/1200)*panelWidth;
numberOfFiguresPerPage        = 3;
numberOfVerticalPlotRows      = numberOfFiguresPerPage;
numberOfHorizontalPlotColumns = 1;    
assert(numberOfVerticalPlotRows*numberOfHorizontalPlotColumns ...
         >= numberOfFiguresPerPage);

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;           
pageHeight  = numberOfVerticalPlotRows*(panelHeight) + 2;
pageWidth   = numberOfHorizontalPlotColumns*(panelWidth) + 2;          
plotHeight  = panelHeight;
plotWidth   = panelWidth ;

plotConfigGeneric;


forceThreshold = 100;
copErrorRecord(2) = struct('opticalCop',  zeros(data.Frames,3),...
                           'forcePlateCop', zeros(data.Frames,3),...
                           'errorCop', zeros(data.Frames,3),...
                           'errorMagCop', zeros(data.Frames,1),...                           
                           'force', zeros(data.Frames,3),...
                           'isValid',zeros(data.Frames,1),...
                           'indexTime',zeros(data.Frames,1));
for k=1:1:2
  copErrorRecord(k).opticalCop    = zeros(data.Frames,3);
  copErrorRecord(k).forcePlateCop = zeros(data.Frames,3);
  copErrorRecord(k).errorCop      = zeros(data.Frames,3);
  copErrorRecord(k).errorMagCop   = zeros(data.Frames,3);
  copErrorRecord(k).force         = zeros(data.Frames,3);
  copErrorRecord(k).isValid       = zeros(data.Frames,3);
  copErrorRecord(k).indexTime     = zeros(data.Frames,3);  
end
                         

idxActiveForcePlate = 0;                         

idxForcePlateZero = 0;

modelMarkers      = 0;
idxMarkerForce    = 0;

flag_configuringTrial=1;

%Trial specific script configuration
if(strfind(filesToProcess{1},'test_Lframe_0001') == 1)
  idxForcePlateZero = 2500*data.Force(1).SamplingFactor;
  modelMarkers = [1,3,5,2];
  idxMarkerForce=1;
  idxActiveForcePlate = 2;
  sf = data.Force(1).SamplingFactor;
  for k=1:1:length(data.Force)
    assert(sf == data.Force(k).SamplingFactor);
  end
  
  modelGeometry = zeros(length(modelMarkers),3);
  markerCount = 1;
  for idxMarker=1:1:data.Trajectories.Unidentified.Count
    flag_isMarker=0;
      for k=1:1:length(modelMarkers)
        if(idxMarker == modelMarkers(1,k))
          flag_isMarker=1;
        end
      end
      if(flag_isMarker==1)
        modelGeometry(markerCount,:) = ...
          [data.Trajectories.Unidentified.Data(idxMarker,1,1),...
           data.Trajectories.Unidentified.Data(idxMarker,2,1),...
           data.Trajectories.Unidentified.Data(idxMarker,3,1)];
        markerCount=markerCount+1;
      end
  end
  save('modelLBracket.mat','modelGeometry');
  flag_configuringTrial=0;
end

if(strfind(filesToProcess{1},'test_Lframe_0002') == 1)
  idxForcePlateZero = 3691*data.Force(1).SamplingFactor;
  modelMarkers = [3,4,6,2];
  idxMarkerForce = 1;
  idxActiveForcePlate = 1;
  sf = data.Force(1).SamplingFactor;
  for k=1:1:length(data.Force)
    assert(sf == data.Force(k).SamplingFactor);
  end
  flag_configuringTrial=0;
end


if(strfind(filesToProcess{1},'test_Lframe_0003') == 1)
  assert(0,'There are no markers visible in this trial');
end

if(strfind(filesToProcess{1},'test_Lframe_20001') == 1)
  idxForcePlateZero = 1*data.Force(1).SamplingFactor;
  modelMarkers = [1,4,3,2];
  idxMarkerForce = 1;
  idxActiveForcePlate = 2;
  sf = data.Force(1).SamplingFactor;
  for k=1:1:length(data.Force)
    assert(sf == data.Force(k).SamplingFactor);
  end
  
  flag_configuringTrial=0;
end

if(strfind(filesToProcess{1},'test_Lframe_20002') == 1)
  idxForcePlateZero = 1*data.Force(1).SamplingFactor;
  modelMarkers     = [2,1,3,5];
  
  idxMarkerForce = 1;
  idxActiveForcePlate = 1;
  sf = data.Force(1).SamplingFactor;
  for k=1:1:length(data.Force)
    assert(sf == data.Force(k).SamplingFactor);
  end
  
  flag_configuringTrial=0; %Not yet setup
end

modelMarkerOrderMap = [1,2,3,4];

if(flag_plotScene==1)
  dTime = round(data.FrameRate*0.2);
end

for idxTime=1:dTime:size(data.Trajectories.Unidentified.Data,3)

  
  
  if(flag_plotScene==1)
    clf(figTest);  
    plot3([0;1].*100,[0;0],[0;0],'r');
    hold on;
    plot3([0;0],[0;1].*100,[0;0],'g');
    hold on;
    plot3([0;0],[0;0],[0;1].*100,'b');
    hold on;
  end

  if(flag_printError==1)
    disp(idxTime);
  end
  %idxTime = modelTimeIndex;
  
  if(idxTime==5207)
    here=1;
  end
  
  posMarkers = zeros(length(modelMarkers),3).*NaN;
  idxValidMarker=1;
  for idxMarker=1:1:data.Trajectories.Unidentified.Count

    mkrColor = 'b';
    flag_isMarker=0;
    idxMarkerInsertion=0;
    for k=1:1:length(modelMarkers)
      if(idxMarker == modelMarkers(1,k))
        flag_isMarker=1;
        idxMarkerInsertion = modelMarkerOrderMap(1,k);
      end
    end

    if(flag_isMarker==1)
      
      if(flag_plotScene==1)
        plot3(data.Trajectories.Unidentified.Data(idxMarker,1,idxTime),... 
              data.Trajectories.Unidentified.Data(idxMarker,2,idxTime),...
              data.Trajectories.Unidentified.Data(idxMarker,3,idxTime),['o',mkrColor]);
        hold on;
        text(data.Trajectories.Unidentified.Data(idxMarker,1,idxTime),... 
              data.Trajectories.Unidentified.Data(idxMarker,2,idxTime),...
              data.Trajectories.Unidentified.Data(idxMarker,3,idxTime),...
              num2str(idxMarker));
        hold on;
        axis equal;
      end
      posMarkers(idxMarkerInsertion,:) ...
        = [data.Trajectories.Unidentified.Data(idxMarker,1,idxTime),...
           data.Trajectories.Unidentified.Data(idxMarker,2,idxTime),...
           data.Trajectories.Unidentified.Data(idxMarker,3,idxTime)];
      idxValidMarker=idxValidMarker+1;
    end
  end
  idxValidMarker = idxValidMarker-1;

  


  if(idxValidMarker >= 4 && sum(sum(isnan(posMarkers))) == 0)
    
    %Get the cop location using the optical measurement
    vec13 = posMarkers(3,:)-posMarkers(1,:);
    vec13 = vec13 ./ norm(vec13);
    vec34= posMarkers(4,:)-posMarkers(3,:);
    vec34 = vec34 ./ norm(vec34);
    
 
    offsetAxis = cross(-vec13,vec34);  
    offsetAxis = offsetAxis./norm(offsetAxis);
    if(offsetAxis(1,3) > 0)
      offsetAxis = offsetAxis.*-1;
    end

    
    copOptical = posMarkers(idxMarkerForce,:) + offsetAxis.*50;
    
    if(flag_plotScene==1)
      plot3(copOptical(1,1),copOptical(1,2),copOptical(1,3),'x')
      hold on;
    end
    
    %Get the cop location from the force plate
    
    for idxFP=1:1:length(data.Force)
      sq = [data.Force(idxFP).ForcePlateLocation;...
            data.Force(idxFP).ForcePlateLocation(1,:)];

      if(flag_plotScene==1)
        fpColor = ['m','g'];
        for idxEdge = 1:1:(size(sq,1)-1)
          plot([sq(idxEdge,1),sq(idxEdge+1,1)],...
               [sq(idxEdge,2),sq(idxEdge+1,2)],...
               fpColor(idxFP))
          hold on;
        end
      end

      idxFpSample = (idxTime)*data.Force(idxFP).SamplingFactor;
      force =  data.Force(idxFP).Force(:,idxFpSample)';%
      forceZero = [0,0,0];
      if(isnan(idxForcePlateZero)==0)
        forceZero = data.Force(idxFP).Force(:,idxForcePlateZero)';
      end
      force = force-forceZero;

      forceLine = [data.Force(idxFP).COP(:,idxFpSample)';...
                  (data.Force(idxFP).COP(:,idxFpSample)'...
                  + force)]; 

      if(flag_plotScene==1)
        lineColor = 'k';
        if(abs(force(1,3)) > forceThreshold)
          lineColor = fpColor(idxFP);
        end
        
        plot3(forceLine(:,1),forceLine(:,2),forceLine(:,3),fpColor(idxFP));
        hold on;
        plot3(forceLine(1,1),forceLine(1,2),forceLine(1,3),'ok');
        hold on;
      end

      %if(flag_printError==1)
      %  fprintf('\t%1.3f N\t%i\n',data.Force(idxFP).Force(3,idxFpSample),idxFP);
      %end

      %Get the error
      copErr = copOptical - data.Force(idxFP).COP(:,idxFpSample)';

      %if(flag_printError==1)
      %  fprintf('\t%1.3f mm\tCop Error\n',norm(copErr));
      %end

      if(abs(force(1,3)) > forceThreshold || flag_configuringTrial==1)

        if(flag_configuringTrial==0)
          if(~(offsetAxis(1,3) < -0.9))
            here=1;
          end
          assert(offsetAxis(1,3) < -0.9); 
        end
        %The L bracket is on the ground: the offset axis should point 
        %nearly perfectly down.
        
        idxFpSample = (idxTime)*data.Force(idxFP).SamplingFactor;

        copErrorRecord(idxFP).isValid(idxTime,1) = 1;
        copErrorRecord(idxFP).opticalCop(idxTime,:) = copOptical;
        copErrorRecord(idxFP).forceCop(idxTime,:)   = data.Force(idxFP).COP(:,idxFpSample)';
        copErrorRecord(idxFP).errorCop(idxTime,:)   = copErr;
        copErrorRecord(idxFP).errorMagCop(idxTime,:)   = norm(copErr);      
        copErrorRecord(idxFP).force(idxTime,:) = force;        
        copErrorRecord(idxFP).indexTime(idxTime,1)=idxTime;
      end


    end    
    
    
  end
  
  if(flag_plotScene==1)
    axis([0,1200,0,400,-200,200]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    pause(0.1);  
  end

end

%Get rid of small blocks
minBlockLength = 10;
for idxFP = 1:1:length(data.Force)
  blockCount = 0;
  for idxTime = 1:1:(length(copErrorRecord(idxFP).isValid))
    if( copErrorRecord(idxFP).isValid(idxTime,1)==1)
      blockCount = blockCount+1;
    else
      %We've gotten to the end of a block: zero it out if its less than
      %the min length
      if(blockCount < minBlockLength && blockCount > 0)
        idxEnd = idxTime-1;
        idxStart= max([1,(idxTime-blockCount-1)]);
        copErrorRecord(idxFP).isValid([idxStart:1:idxEnd],1) = 0;
      end
      
      blockCount = 0;
    end
    
    
  end
end


timeVec = [1:1:data.Frames];%.*(1./data.FrameRate);

figCopError = figure;
figure(figCopError);

[row,col] = find(subPlotPanelIndex==1);          
subPlotVec = reshape(subPlotPanel(row,col,:),1,4);   
subplot('Position',subPlotVec);
  lineErrX = plot( timeVec,...
                    copErrorRecord(idxActiveForcePlate).errorCop(:,1),'r');
  hold on;
  lineErrY = plot( timeVec,...
                    copErrorRecord(idxActiveForcePlate).errorCop(:,2),'g');
  hold on;
  lineErrZ = plot( timeVec,...
                   copErrorRecord(idxActiveForcePlate).errorCop(:,3),'b');
  hold on;

  minCop = min(min(copErrorRecord(idxActiveForcePlate).errorCop(:,:)));
  maxCop = max(max(copErrorRecord(idxActiveForcePlate).errorCop(:,:)));      
  
  flag_newValidData=1;
  newValidDataCount = 1;
  
  for i=1:1:data.Frames

    if(copErrorRecord(idxActiveForcePlate).isValid(i) == 1)          
      if(flag_newValidData==1)
        plot([timeVec(i),timeVec(i)],[minCop,maxCop],'Color',[1,1,1].*0.5);
        hold on;
        text(timeVec(i),maxCop,num2str(newValidDataCount),...
             'VerticalAlignment','top');
        hold on;
        newValidDataCount = newValidDataCount+1;
        flag_newValidData = 0;
      end
    else
      flag_newValidData = 1;
    end

  end      
  
  legend( [lineErrX,lineErrY,lineErrZ],...
          ['err-x'],...
          ['err-y'],...
          ['err-z'],'Location','SouthEast');    
    
  hold on;
  grid on;
  box off;
  xlabel('Sample');
  ylabel('Distance (mm)');
  title(['COP Error: ',data.Force(idxActiveForcePlate).ForcePlateName]);


[row,col] = find(subPlotPanelIndex==2);          
subPlotVec = reshape(subPlotPanel(row,col,:),1,4);   
subplot('Position',subPlotVec);

  lineFx = plot( timeVec,...
        copErrorRecord(idxActiveForcePlate).force(:,1),'r');
  hold on;
  lineFy = plot( timeVec,...
        copErrorRecord(idxActiveForcePlate).force(:,2),'g');
  hold on;
  lineFz = plot( timeVec,...
        copErrorRecord(idxActiveForcePlate).force(:,3),'b');
  hold on;



  minForce = min(min(copErrorRecord(idxActiveForcePlate).force(:,:)));
  maxForce = max(max(copErrorRecord(idxActiveForcePlate).force(:,:)));      
  
  flag_newValidData=1;
  newValidDataCount = 1;
  
  for i=1:1:data.Frames

    if(copErrorRecord(idxActiveForcePlate).isValid(i) == 1)          
      if(flag_newValidData==1)
        plot([timeVec(i),timeVec(i)],[minForce,maxForce],'Color',[1,1,1].*0.5);
        hold on;
        text(timeVec(i),maxForce,num2str(newValidDataCount),...
             'VerticalAlignment','top');
        hold on;
        newValidDataCount = newValidDataCount+1;
        flag_newValidData = 0;
      end
    else
      flag_newValidData = 1;
    end

  end        
         
  hold on;
  grid on;
  box off;
  
  legend( [lineFx,lineFy,lineFz],...
          'Fx',...
          'Fy',...
          'Fz','Location','SouthEast'); 
 
  
  xlabel('Sample');
  ylabel('Force (N)');
  title(['Applied Force: ',data.Force(idxActiveForcePlate).ForcePlateName]);
  
[row,col] = find(subPlotPanelIndex==3);          
subPlotVec = reshape(subPlotPanel(row,col,:),1,4);   
subplot('Position',subPlotVec);  


  fpColor = ['m','g'];
  for idxFP = 1:1:length(data.Force)
    sq = [data.Force(idxFP).ForcePlateLocation;...
          data.Force(idxFP).ForcePlateLocation(1,:)];
    for idxEdge = 1:1:(size(sq,1)-1)
      plot([sq(idxEdge,1),sq(idxEdge+1,1)],...
           [sq(idxEdge,2),sq(idxEdge+1,2)],...
           fpColor(idxFP));
      hold on;
    end  
  end

  flag_newValidData=1;
  newValidDataCount = 1;
  for i=1:1:data.Frames

    if(copErrorRecord(idxActiveForcePlate).isValid(i) == 1)    
      ptA = copErrorRecord(idxActiveForcePlate).opticalCop(i,1:2);
      ptB = copErrorRecord(idxActiveForcePlate).forceCop(i,1:2);

      rAB = ptB-ptA;
      ptC = ptA + rAB.*10;
      
      plot( [ptA(1,1),ptC(1,1)],[ptA(1,2),ptC(1,2)],'Color',[1,1,1].*0.5);
      hold on;
      plot(ptA(1,1),ptA(1,2),'.','Color',[0,0,0]);
      hold on;

      if(flag_newValidData==1)
        text(ptA(1,1)+10,ptA(1,2)+10,num2str(newValidDataCount));
        newValidDataCount = newValidDataCount+1;
        flag_newValidData = 0;
      end
    else
      flag_newValidData = 1;
    end

  end

  %axis equal;
  grid on;
  box off;
  axis([0,1200,0,400]);
  xlabel('X (mm)');
  ylabel('Y (mm)');
  title('Vector from optical (.) to force-plate COP (x10 mag)');

figure(figCopError);  
configPlotExporter;



print('-dpdf',[filesToProcess{1},'.pdf']);  

