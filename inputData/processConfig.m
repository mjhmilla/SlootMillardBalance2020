%Standard Trials collected
index_Static      = 1;
index_Chest       = 2;
index_Conv        = 3; %What does this stand for?
index_Leg         = 4;
index_Side        = 5;
index_Rob         = 6; %Rob?

numberOfDataTrials = 5;
offsetTrialIndex = 1;

index_FeetForcePlate = 1;
index_ChairForcePlate = 2;

trialTypeNames = {'Static','Chest','Conv','Leg','Side','Rob'};
numberOfTrialTypes = length(trialTypeNames);

%Keywords that distinguish C3D data from each trial
%  and the folder that contains the corresponding data from Visual 3D
c3dFileKeyWords = {  'static', 'Che','Con', 'Leg','Sid','Rob'};
folderKeyWords  = {        '', 'Che','Conv','Leg','Sid','Rob'};

%Get the list of files and directories in the data folder
rootDir = pwd;
cd(inputFolder)
  list = dir;
cd(rootDir);

numberOfC3DFiles = 0;
for i=1:1:length(list)
  if(list(i).isdir == 0)
    if(isempty(strfind(list(i).name,'.c3d'))==0)
      numberOfC3DFiles = numberOfC3DFiles +1;
    end
  end
end

if(numberOfTrialTypes < numberOfC3DFiles)
  numberOfTrialTypes = numberOfC3DFiles;
end


inputSubjectFolder = { [inputPath,'/',inputFolder] };
outputSubjectFolder = {[outputPath,'/',outputFolder] };

if(exist(outputSubjectFolder{1},'dir')==0)
  [success,msg,msgid] = mkdir(outputPath,outputFolder);  
end
inputC3DFolders     = cell(numberOfTrialTypes,1);
inputC3DFiles       = cell(numberOfTrialTypes,1);

inputV3DFolders     = cell(numberOfTrialTypes,1);
inputWholeBodyFiles = cell(numberOfTrialTypes,1);
inputAnthroFiles    = cell(numberOfTrialTypes,1);

outputTrialFolders      = cell(numberOfTrialTypes,1);
outputMeshupGrfFiles    = cell(numberOfTrialTypes,1);
outputMeshupFpeFiles    = cell(numberOfTrialTypes,1);
outputMeshupCapFiles    = cell(numberOfTrialTypes,1);
outputMeshupComFiles    = cell(numberOfTrialTypes,1);

outputFpeFileNames = cell(numberOfTrialTypes,1);
outputCapFileNames = cell(numberOfTrialTypes,1);
outputSegmentationFileNames = cell(numberOfTrialTypes,1);
outputMovementSequenceFileNames = cell(numberOfTrialTypes,1);

outputFpeToFootHullDistanceFileNames   = cell(numberOfTrialTypes,1);
outputCapToFootHullDistanceFileNames   = cell(numberOfTrialTypes,1);
outputComGPToFootHullDistanceFileNames = cell(numberOfTrialTypes,1);
outputCopToFootHullDistanceFileNames   = cell(numberOfTrialTypes,1);

%Model factory output file settings
outputModelFactoryAnthropometryFile = ['modelFactoryAnthropometry'];
outputModelFactoryEnvironmentFile = ['modelFactoryEnvironment.env'];
flag_modelFactoryAddMarkers = 1;
modelFactoryLuaModel = [subjectId,'.lua'];



modelFactoryCommonFileDirectory = [outputPath,'/ModelFactoryModelFiles/'];
modelFactoryDescriptionFile = ['3DHumanHeiAge_Description'];

[success,message,messageId] = ...
  copyfile([modelFactoryCommonFileDirectory,modelFactoryDescriptionFile],...
            [outputSubjectFolder{1},modelFactoryDescriptionFile],'f');

modelFactoryScalingAlgorithm = 'deLeva1996_segmentedTrunk';

for i=1:1:length(list)
  
  if(list(i).isdir ==0)
    for j=1:1:length(c3dFileKeyWords)
      if(    isempty( strfind(list(i).name,c3dFileKeyWords{j}))==0 ...
          && isempty( strfind(list(i).name,'.c3d'            ))==0 )
        
        k=0;
        flag_indexFound = 0;
        if(j==index_Rob)
          for z=0:1:(numberOfC3DFiles - j)
            if(isempty(inputC3DFolders{j+z}) == 1 && flag_indexFound == 0)
              k=z;
              flag_indexFound = 1;
            end
          end
        end
        
        inputC3DFolders{j+k}  = [list(i).folder,'/'];
        inputC3DFiles{j+k}    = [list(i).name];
        
        idx = strfind(list(i).name,'.');        
        fname = list(i).name(1:1:(idx-1));               
        outputMeshupGrfFiles{j+k} = ['meshupGrf.ff'];
        outputMeshupFpeFiles{j+k} = ['meshupFpe.csv'];
        outputMeshupCapFiles{j+k} = ['meshupCap.csv'];
        outputMeshupComFiles{j+k} = ['meshupCom.csv'];
      end
    end        
  end  
end

for j=1:1:numberOfTrialTypes

  
  
  fname = inputC3DFiles{j};
  if(isempty(fname) ==0)
    if(  contains(fname,'static') == 0)
      %We have to extract the keyword snippet out of the end of the name.
      %This is actually pretty nasty because who ever named the files did
      %not choose a standard form. WTF


      idx0 = 1;  
      idx1 = max(strfind(fname,'.'));
      for k=1:1:idx1
        tmp = fname(1,k);
        for z=1:1:10
          if(strcmp(num2str(z-1),tmp))
            idx0 = k;
          end        
        end
      end

      keyword = fname((idx0+2):1:(idx1-1));

      outputTrialFolders{j} = [outputPath,'/',outputFolder,keyword,'/']; 
      inputWholeBodyFiles{j} = ['DATA.txt'];
      inputAnthroFiles{j}    = ['SUBMETRICS.txt'];         
      outputFpeFileNames{j} = 'fpe.mat'; 
      outputCapFileNames{j} = 'cap.mat';
      outputSegmentationFileNames{j} = 'motionSegmentation.mat';
      outputMovementSequenceFileNames{j} = 'motionSequence.mat';

      outputFpeToFootHullDistanceFileNames{j} = 'fpe2FootConvexHullDist.mat';
      outputCapToFootHullDistanceFileNames{j} = 'cap2FootConvexHullDist.mat';
      outputComGPToFootHullDistanceFileNames{j} = 'comgp2FootConvexHullDist.mat';
      outputCopToFootHullDistanceFileNames{j} = 'cop2FootConvexHullDist.mat';


      idxInputFolder = 0;
      bestNumberOfMatchingCharacters = 0;
      for k=1:1:length(list)
        if(list(k).isdir == 1)
          s1 = keyword;
          s2 = list(k).name;
          n = min(length(s1),length(s2));

          score = zeros(1,n);

          for z=1:1:n
            score(1,z) = strcmp(s1(1,z),s2(1,z));
          end

          if(sum(score) > bestNumberOfMatchingCharacters)
            idxInputFolder = k;
            bestNumberOfMatchingCharacters = sum(score);
          end
        end
      end
      assert(idxInputFolder ~= 0, ...
        sprintf('No folder with the exact keyword (%s) from c3d file (%s)\n',...
                keyword,fname));

      inputV3DFolders{j} = ...
         [list(idxInputFolder).folder,'/',list(idxInputFolder).name,'/'];


    else
      outputTrialFolders{j} = [outputPath,'/',outputFolder];            
    end
  end
end

