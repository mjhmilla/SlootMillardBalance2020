function [frameLeft, frameRight]=getFootFrames(index, mkrPos, mkrNames,...
                                        frameLeftOffset, frameRightOffset)


frameRight = struct('r',zeros(3,1),'E',zeros(3,3));
frameLeft  = struct('r',zeros(3,1),'E',zeros(3,3));

frameRight.r = (1/2).*(...%mkrPos.('R_FCC')(index,:)'...
                     + mkrPos.('R_FAL')(index,:)'...
                     + mkrPos.('R_TAM')(index,:)');
           

ey = 0.5.*(mkrPos.('R_FM1')(index,:)+mkrPos.('R_FM5')(index,:)) ...
        - mkrPos.('R_FCC')(index,:); 
ey = ey./norm(ey);

ex = mkrPos.('R_FAL')(index,:) - mkrPos.('R_TAM')(index,:);
ex = ex - sum(ey.*ex).*ey;
ex = ex./norm(ex);

ez = cross(ex,ey);

assert(ex*ey' < 1e-6);
assert(ex*ez' < 1e-6);


frameRight.E = frameRightOffset.E*[ex', ey', ez'];

frameRight.r = frameRight.r + frameRight.E*frameRightOffset.r;


frameLeft.r = (1/2).*(...%mkrPos.('L_FCC')(index,:)'...
                     + mkrPos.('L_FAL')(index,:)'...
                     + mkrPos.('L_TAM')(index,:)');



ey = 0.5.*(mkrPos.('L_FM1')(index,:)+mkrPos.('L_FM5')(index,:)) ...
        - mkrPos.('L_FCC')(index,:); 
ey = ey./norm(ey);

ex = mkrPos.('L_TAM')(index,:) - mkrPos.('L_FAL')(index,:);
ex = ex - sum(ey.*ex).*ey;
ex = ex./norm(ex);

ez = cross(ex,ey);

assert(ex*ey' < 1e-6);
assert(ex*ez' < 1e-6);

frameLeft.E = frameLeftOffset.E*[ex', ey', ez'];

frameLeft.r = frameLeft.r + frameLeft.E*frameLeftOffset.r;


