clc;
close all;
clear all;

r000 = [0;0;0];

r0C0 = [0;0;1];
v0C0 = [0;0;0];

mass  =1;
JC0 = eye(3,3);
w0C0a = [0;0;0];
w0C0b = [0;1;0];

HC0a = JC0*w0C0a;   
HC0b = JC0*w0C0b;   

g0 = [0;0;-9.81];

tol     = 1e-9;
iterMax = 100;
flag_fpeEvaluateDerivatives = 0;
flag_fpeVerbose = 0;

vx = [0:0.01:1]';


phiA = zeros(length(vx),1);
phiB = zeros(length(vx),1);

fpeInfoAPrevious = [];
fpeInfoBPrevious = [];

for j=1:1:length(vx)

  v0C0(1,1) = vx(j,1);
  
  
  fpeInfoA = calc3DFootPlacementEstimatorInfo(mass,...
                                          r0C0,...
                                          v0C0,...                                                    
                                          JC0,...                                                    
                                          HC0a,...
                                          r000,...
                                          g0,...
                                          tol,...
                                          iterMax,...
                                          flag_fpeEvaluateDerivatives,...
                                          flag_fpeVerbose,...
                                          fpeInfoAPrevious);
  fpeInfoAPrevious = fpeInfoA;
  
  phiA(j,1) = fpeInfoA.phi;
  
  fpeInfoB = calc3DFootPlacementEstimatorInfo(mass,...
                                          r0C0,...
                                          v0C0,...                                                    
                                          JC0,...                                                    
                                          HC0b,...
                                          r000,...
                                          g0,...
                                          tol,...
                                          iterMax,...
                                          flag_fpeEvaluateDerivatives,...
                                          flag_fpeVerbose,...
                                          fpeInfoBPrevious);  
  fpeInfoBPrevious = fpeInfoB;
  phiB(j,1) = fpeInfoB.phi;

end

fig = figure;
plot(vx,phiA(:,1),'b');
hold on;
plot(vx,phiB(:,1),'r');
xlabel('Forward Velocity (m/s)');
ylabel('Leg angle (phi)');
legend('omega = 0', 'omega > 0');
