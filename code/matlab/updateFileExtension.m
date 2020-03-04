function fileNameUpd = updateFileExtension(fileName, newExt)

idx = strfind(fileName,'.');
fileNameUpd = [fileName(1,1:1:idx),newExt];