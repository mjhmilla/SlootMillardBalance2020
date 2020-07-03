clc;
close all;
clear all;
 
filesToProcess = {'test_cilinder_0002'};

 
data =load([filesToProcess{1},'.mat']);
data =data.(filesToProcess{1});


assert(data.Trajectories.Labeled.Count==0)

mkrCylinderBottom = [2,5,7,13,15,16];
mkrCylinderTop    = [4,6,10,11,12,8];
mkrTop      = [3,14,9];

figTest=figure;

for i=1:1:data.Trajectories.Unidentified.Count
  idxT =1;
  mkrColor = 'b';
  flagIsTop =0;
  flagIsCylinderBottom=0;
  flagIsCylinderTop=0;
  
  for k=1:1:length(mkrTop)
    if  mkrTop(1,k) == i
      flagIsTop=1;
      mkrColor = 'r';
    end
  end
  for k=1:1:length(mkrCylinderBottom)
    if  mkrCylinderBottom(1,k) == i
      flagIsCylinderBottom=1;
      mkrColor = 'b';
    end
  end
  for k=1:1:length(mkrCylinderTop)
    if  mkrCylinderTop(1,k) == i
      flagIsCylinderTop=1;
      mkrColor = 'g';
    end
  end
  
  if(flagIsTop == 1 || flagIsCylinderTop==1 || flagIsCylinderBottom==1)
    plot3(data.Trajectories.Unidentified.Data(i,1,idxT),... 
          data.Trajectories.Unidentified.Data(i,2,idxT),...
          data.Trajectories.Unidentified.Data(i,3,idxT),['o',mkrColor]);
    hold on;
    text(data.Trajectories.Unidentified.Data(i,1,idxT),... 
          data.Trajectories.Unidentified.Data(i,2,idxT),...
          data.Trajectories.Unidentified.Data(i,3,idxT),...
          num2str(i));
    hold on;
    axis equal
  end
end