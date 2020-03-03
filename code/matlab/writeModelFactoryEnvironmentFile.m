function success  = writeModelFactoryEnvironmentFile(...
                      fileName, ...
                      flag_addMarkers, ...
                      luaFileName,...
                      modelDescriptionFile,...                      
                      scalingAlgorithm,...
                      anthropometryFileName)
                    
                    
fid = fopen(fileName,'w');                  
success = 0;

fprintf(fid,'AddMarkers\n');
  assert(flag_addMarkers == 0 || flag_addMarkers == 1);
fprintf(fid,'%d\n',flag_addMarkers);
fprintf(fid,'humanModel_Save\n');
  assert( contains(luaFileName,'.lua'));
fprintf(fid,'%s\n',luaFileName);
fprintf(fid,'humanModel_DescriptionFile\n');
fprintf(fid,'%s\n',modelDescriptionFile);
fprintf(fid,'humanModel_ScalingAlgorithmChoice\n');
fprintf(fid,'%s\n',scalingAlgorithm);
fprintf(fid,'humanModel_AnthropometryFile\n');
fprintf(fid,'%s\n',anthropometryFileName);

fclose(fid);
success = 1;
                    