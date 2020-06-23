function [frameLeftOffset, frameRightOffset]=getFootOffsetFrames(index, mkrPos, ...
                         mkrNames)


frameRightOffset = struct('r',zeros(3,1),'E',zeros(3,3));
frameLeftOffset  = struct('r',zeros(3,1),'E',zeros(3,3));

%Right frame
r = (1/2).*(...%mkrPos.('R_FCC')(index,:)'...
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


%Right frame desired orientation
exD = ex;
exD(1,3) = 0;
exD = exD./norm(exD);

eyD = ey;
eyD(1,3) = 0;
eyD = eyD./norm(eyD);

ezD = cross(exD,eyD);

%Offset
frameRightOffset.E = [exD',eyD',ezD']*([ex',ey',ez']');
frameRightOffset.r = [0;0;-r(3,1)];

%Left frame
r = (1/2).*(...%mkrPos.('L_FCC')(index,:)'...
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


%Left frame desired orientation
exD = ex;
exD(1,3) = 0;
exD = exD./norm(exD);

eyD = ey;
eyD(1,3) = 0;
eyD = eyD./norm(eyD);

ezD = cross(exD,eyD);
frameLeftOffset.E = [exD',eyD',ezD']*([ex',ey',ez']');
frameLeftOffset.r = [0;0;-r(3,1)];



