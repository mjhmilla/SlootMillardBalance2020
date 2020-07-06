function [bosLeftSym,bosRightSym] = getFunctionalFootBaseOfSupport(index,...
                                mkrPos, mkrNames,frameLeft,frameRight)

bosLeftSym = [];
bosRightSym= [];


%Left foot: markers moved into the boundary of the COP region
tmp     = frameLeft.E'*(mkrPos.('L_FCC')(index,:)'-frameLeft.r);
lfcc    = [0,(4.5/6.5)*tmp(2,1),0];

tmp  = frameLeft.E'*(mkrPos.('L_TAM')(index,:)'-frameLeft.r);
ltam = [tmp(1,1),tmp(2,1),0].*(2.0/4.);

tmp   = frameLeft.E'*(mkrPos.('L_FAL')(index,:)'-frameLeft.r);
lfal = [tmp(1,1),tmp(2,1),0].*(2.0/4.);

lfm15  = 0.5.*(mkrPos.('L_FM1')(index,:)'+mkrPos.('L_FM5')(index,:)');


lfm1 = frameLeft.E'*((4./5.0).*(mkrPos.('L_FM1')(index,:)'-lfm15) ...
                    +0.9*(lfm15-frameLeft.r));
lfm1 = [lfm1(1,1),lfm1(2,1),0];

lfm50 = frameLeft.E'*((3.0./4).*(mkrPos.('L_FM5')(index,:)'-lfm15)...
                     +(9/11.5)*(lfm15-frameLeft.r));
lfm50 = [lfm50(1,1),lfm50(2,1),0];

lfm51 = frameLeft.E'*((2.0./4).*(mkrPos.('L_FM5')(index,:)'-lfm15)...
                     +(12/11.5)*(lfm15-frameLeft.r));
lfm51 = [lfm51(1,1),lfm51(2,1),0];

tmp     = frameLeft.E'*(mkrPos.('L_FM2')(index,:)'-frameLeft.r);
lfm20     = [tmp(1,1),tmp(2,1),0].*(15./19.);

lfm22 = lfm1;
lfm22(1,2) = lfm20(1,2);
lfm21 = (lfm22+lfm20).*0.5;

%lfm20(1,2) = lfm20(1,2).*0.95;

lfm20(1,2)=lfm20(1,2).*1.15;
lfm21(1,2)=lfm21(1,2).*1.15;
lfm22(1,2)=lfm22(1,2).*1.15;

ltam1 = 0.5.*ltam;
ltam1(1,2) = lfcc(1,2);

lfal1 = 0.5.*lfal;
lfal1(1,2) = lfcc(1,2);



bosLeft = [ltam1;ltam;lfm1;lfm22;lfm21;lfm51;lfm50;lfal;lfal1];


%Right foot: markers moved into the boundary of the COP region
tmp     = frameRight.E'*(mkrPos.('R_FCC')(index,:)'-frameRight.r);
rfcc    = [0,(4.5/6.5)*tmp(2,1),0];

tmp  = frameRight.E'*(mkrPos.('R_TAM')(index,:)'-frameRight.r);
rtam = [tmp(1,1),tmp(2,1),0].*(2.0./4.);

tmp   = frameRight.E'*(mkrPos.('R_FAL')(index,:)'-frameRight.r);
rfal = [tmp(1,1),tmp(2,1),0].*(2.0/4.);

rfm15  = (0.5.*(mkrPos.('R_FM1')(index,:)'+mkrPos.('R_FM5')(index,:)'));


rfm1 = frameRight.E'*((4./5.0).*(mkrPos.('R_FM1')(index,:)'-rfm15)...
                    +0.9*(rfm15-frameRight.r));
rfm1 = [rfm1(1,1),rfm1(2,1),0];

rfm50 = frameRight.E'*((3.0./4).*(mkrPos.('R_FM5')(index,:)'-rfm15)...
                     +(9/11.5)*(rfm15-frameRight.r));
rfm50 = [rfm50(1,1),rfm50(2,1),0];

rfm51 = frameRight.E'*((2.0./4).*(mkrPos.('R_FM5')(index,:)'-rfm15)...
                     +(12/11.5)*(rfm15-frameRight.r));
rfm51 = [rfm51(1,1),rfm51(2,1),0];


tmp     = frameRight.E'*(mkrPos.('R_FM2')(index,:)'-frameRight.r);
rfm20     = [tmp(1,1),tmp(2,1),0].*(15./19.);

rfm22 = rfm1;
rfm22(1,2) = rfm20(1,2);

rfm21 = (rfm22+rfm20).*0.5;
%rfm22(1,2) = rfm22(1,2).*0.95;

rfm20(1,2)=rfm20(1,2).*1.15;
rfm21(1,2)=rfm21(1,2).*1.15;
rfm22(1,2)=rfm22(1,2).*1.15;



rtam1 = 0.5.*rtam;
rtam1(1,2) = rfcc(1,2);

rfal1 = 0.5.*rfal;
rfal1(1,2) = rfcc(1,2);


bosRight = [rtam1;rtam;rfm1;rfm22;rfm21;rfm51;rfm50;rfal;rfal1];

bosRightFlip = [-bosRight(:,1),bosRight(:,2),bosRight(:,3)];

bosMean = 0.5.*(bosLeft+bosRightFlip);


bosLeftSym = bosMean;
bosRightSym= [-bosMean(:,1),bosMean(:,2),bosMean(:,3)];

here=1;
