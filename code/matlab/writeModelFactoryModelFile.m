function success = writeModelFactoryModelFile(...
                    fileName,...
                    ageYears,...
                    heightInMeters,...
                    massInKg, ...
                    gender1Male0Female, ...
                    pelvisASISWidthInMeters, ...
                    pelvisHipJointWidthInMeters, ...
                    shoulderCenterWidthInMeters, ...
                    heelToAnkleLengthInMeters,...
                    ankleHeightToGroundInMeters,...
                    shoulderCenterToC7VerticalDistanceInMeters,...
                    torsoDepthInMeters,...
                    footWidthInMeters)
                  
fid = fopen(fileName,'w');                  
                  
success  =0;

fprintf(fid,'# Data specifying anthropometric details of subject\n');
fprintf(fid,'# Age in years\n');
fprintf(fid,'age, %1.1f\n',ageYears);
fprintf(fid,'# Height in m\n');
fprintf(fid,'height, %1.3f\n',heightInMeters);
fprintf(fid,'# Weight in kg\n');
fprintf(fid,'weight, %1.3f\n',massInKg);
fprintf(fid,'# Gender (1 = Male, 0 = Female)\n');
fprintf(fid,'gender, %d\n',gender1Male0Female);
fprintf(fid,'# Maximum width of pelvis in m\n');
fprintf(fid,'pelvisWidth, %1.3f\n',pelvisASISWidthInMeters);
fprintf(fid,'# Distance between hip centers in m\n');
fprintf(fid,'hipCenterWidth, %1.3f\n',pelvisHipJointWidthInMeters);
fprintf(fid,'# Distance between shoulder centers in m\n');
fprintf(fid,'shoulderCenterWidth, %1.3f\n',shoulderCenterWidthInMeters);
fprintf(fid,'# Distance between heel and ankle center along foot length in m\n');
fprintf(fid,'heelAnkleXOffset, %1.3f\n',heelToAnkleLengthInMeters);
fprintf(fid,'# Vertical distance between heel and ankle center in m\n');
fprintf(fid,'heelAnkleZOffset, %1.3f\n',ankleHeightToGroundInMeters);
fprintf(fid,'# Vertical distance between shoulder centers and C7 in m\n');
fprintf(fid,'shoulderNeckZOffset, %1.3f\n',...
          shoulderCenterToC7VerticalDistanceInMeters);
fprintf(fid,'# Maximum width of foot in m\n');
fprintf(fid,'footWidth, %1.3f\n',footWidthInMeters);
fprintf(fid,'# Depth of the torso between TV7 and SXS in m\n');
fprintf(fid,'torsoDepth, %1.3f\n',torsoDepthInMeters);

fclose(fid);
     
success  =1;
                  
                    