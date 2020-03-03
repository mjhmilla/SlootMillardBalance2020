subjectAge                  = 25;
subjectGender1Male0Female   = 0;

subjectId = 'H10';
inputFolder  = [subjectId,'/'];
outputFolder = [subjectId,'/'];

dataDir = pwd;
cd(inputFolder)
%Rename these c3d files so they match the directories
  fileNamesOrig = {'sts_0005_Rob2.c3d'};
  fileNamesUpd  = {'sts_0005_Rob.c3d'};

  list = dir;
  for k=1:1:length(list)
    if list(k).isdir == 0
      for z=1:1:length(fileNamesOrig)
        if strcmp(list(k).name,fileNamesOrig{z})==1
          movefile(fileNamesOrig{z},fileNamesUpd{z});
        end
      end
    end
  end
cd(dataDir);

processConfig;
