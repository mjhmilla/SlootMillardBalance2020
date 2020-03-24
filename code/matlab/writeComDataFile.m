function success = writeComDataFile(...
                      fname, c3dTime, comPosition,gravityVector, nanVal)
                    
success = 0;

fid       = fopen(fname,'w');
nVector   = reshape(-gravityVector./norm(gravityVector),1,3);

for i=1:1:(length(c3dTime))

  if(  sum(isnan(comPosition(i,:))) == 0)
    
    comPositionGP=comPosition(i,:)-(comPosition(i,:).*nVector).*nVector;

    fprintf(fid,'%1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f\n',...
            c3dTime(i,1), ...
        comPosition(i,1), comPosition(i,2),  comPosition(i,3),...
      comPositionGP(1,1),comPositionGP(1,2),comPositionGP(1,3));

  else
    fprintf(fid,'%1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f, %1.6f\n',...
            c3dTime(i,1), ...
            nanVal, nanVal, nanVal,...
            nanVal, nanVal, nanVal);

  end
end

fclose(fid);

success = 1;                    