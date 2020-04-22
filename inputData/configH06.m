subjectAge                  = 25;
subjectGender1Male0Female   = 0;

forcePlateDataRecorded = 1;

subjectId = 'H06';
inputFolder  = [subjectId,'/'];
outputFolder = [subjectId,'/'];


dataDir = pwd;
cd(inputFolder)
%Rename these c3d files so they match the directories
  fileNamesOrig = {'sts_07_Rob_impr.c3d'};
  fileNamesUpd  = {'sts_07_Rob_instr.c3d'};

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
