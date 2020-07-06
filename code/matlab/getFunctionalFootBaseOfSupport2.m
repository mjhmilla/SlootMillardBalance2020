function [bosLeft,bosRight] = getFunctionalFootBaseOfSupport2(index,...
                                mkrPos, mkrNames,frameLeft,frameRight)

%L :center of the left foot frame
%R : center of the right foot frame
rLHL =  frameLeft.E'*(mkrPos.('L_FCC')(index,:)'-frameLeft.r);
rRHR = frameRight.E'*(mkrPos.('R_FCC')(index,:)'-frameRight.r);

lH = 0.5*(abs(rRHR(2,1))+abs(rLHL(2,1)));

rMTL =  frameLeft.E'*(mkrPos.('L_TAM')(index,:)'-mkrPos.('L_FAL')(index,:)');
rMTR = frameRight.E'*(mkrPos.('R_TAM')(index,:)'-mkrPos.('R_FAL')(index,:)');

wA = 0.5*(abs(rMTL(1,1))+abs(rMTR(1,1)));

rLFL = frameLeft.E'*( 0.5.*(mkrPos.('L_FM1')(index,:)'...
                           +mkrPos.('L_FM5')(index,:)') - frameLeft.r);
rRFR = frameRight.E'*(0.5.*(mkrPos.('R_FM1')(index,:)'...
                           +mkrPos.('R_FM5')(index,:)') - frameRight.r);

lF = 0.5*(rLFL(2,1)+rRFR(2,1));

r51L =  frameLeft.E'*(mkrPos.('L_FM5')(index,:)'-mkrPos.('L_FM1')(index,:)');
r51R = frameRight.E'*(mkrPos.('R_FM5')(index,:)'-mkrPos.('R_FM1')(index,:)');

wF = 0.5*(abs(r51L(1,1))+abs(r51R(1,1)));

rLTL =  frameLeft.E'*(mkrPos.('L_FM2')(index,:)'-frameLeft.r);
rRTR = frameRight.E'*(mkrPos.('R_FM2')(index,:)'-frameRight.r);

lT = 0.5*( rLTL(2,1) + rRTR(2,1) );

%Medial toe offset
oTm = 0.5*(abs(rLTL(1,1))+abs(rRTR(1,1)));


%Generic left foot model

%Reference foot
lHRef  =      6.078755029885653e-02;
wARef  =      9.911568159439875e-02;
lFRef  =      1.184518482349464e-01;
wFRef  =      9.911568159439875e-02;
lTRef  =      1.935116252459572e-01;
oTmRef =      4.548369119558601e-03;


delta = 0.005;

%Reference Cop Movement
lHCopRef = 0.03136+2*delta;
wHaCopRef = 0.004+delta;
wHbCopRef = 0.008;

tamCopRef = 0.01183 +delta*2;
wFCopRef  = 0.02779 + delta;
lFCopRef  = 0.1167;



w2aCopRef = 0.02705+delta*3;
l2aCopRef = 0.1492 +delta*5;

w2bCopRef = 0.01806 - 2*delta;
l2bCopRef = 0.1492  + 2*delta ;

w5aCopRef =  0.02052;
l5aCopRef =  0.111;

w5bCopRef =  0.02816+delta;
l5bCopRef =  0.09321-delta;

falCopRef = 0.01603+delta*2;



fcca= [ (wHaCopRef/wARef)*wA, -(lHCopRef/lHRef)*lH,0];
fccb= [-(wHbCopRef/wARef)*wA, -(lHCopRef/lHRef)*lH,0];

tam = [(tamCopRef/wARef)*wA,0,0];

m1  = [(wFCopRef/wFRef)*wF,(lFCopRef/lFRef)*lF,0];

m2a = [(w2aCopRef/wFRef)*wF,(l2aCopRef/lTRef)*lT, 0];
m2b = [(w2bCopRef/wFRef)*wF,(l2bCopRef/lTRef)*lT,0];

m5a = [-(w5aCopRef/wFRef)*wF,(l5aCopRef/lFRef)*lF,0];
m5b = [-(w5bCopRef/wFRef)*wF,(l5bCopRef/lFRef)*lF,0];

fal = [-(falCopRef/wARef)*wF,0,0];

bosLeft = [fcca;tam;m1;m2a;m2b;m5a;m5b;fal;fccb];

bosRight = bosLeft.*[-1,1,1];

here=1;

