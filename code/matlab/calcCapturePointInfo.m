function capInfo = calcCapturePointInfo(mass, r0C0, v0C0, r0S0, g0)
%%
% Evaluates the capture point as detailed in Pratt et al., but follows the
% notation that appears in Millard et al. This implementation assumes
% that the foot is contacting a half plane that goes through the point 
% r0S0
%
% Pratt, J., Carff, J., Drakunov, S., and Goswami, A., 2006, “Capture Point: A
% Step Toward Humanoid Push Recovery,” Proceedings of the IEEE-RAS Inter-
% national Conference on Humanoid Robots, Genoa, Italy.
%
% Millard M, McPhee J, Kubica E. Foot placement and balance in 3D. Journal 
% of computational and nonlinear dynamics. 2012 Apr 1;7(2).
%
% @param mass: total mass of the human/robot
% @param r0C0: position vector from the inertial frame to the center-of-mass
%             resolved in the coordinates of the inertial frame
% @param v0C0: velocity vector " ... "
%
% @param  r0S0 (contact surface position): An 3-by-1 vector containing one 
%               position of the surface that the object will be contacting in 
%               order to balance.
%
% @param g0: a 3-by-1 vector of gravity resolved in the inertial frame
%
% @return a struct with the following fields
%         r0F0: location of the capture point
%         r0G0: location of the center-of-mass ground projection
%         lu  : length of the capture-point-step in the u direction from
%               the center-of-mass ground projection r0G0
%         u   : horizontal direction of travel
%         k   : vertical vector
%         vu  : velocity in the direction of travel
%         vk  : velocity in the vertical direction
%         h   : height of the center-of-mass above the ground plane
%         eorb: a quantity that is assumed to be conserved
%%
capInfo = struct( 'r0F0', NaN.*zeros(3,1),...
                  'r0G0', NaN.*zeros(3,1),...  
                    'lu', NaN,...
                    'vu', NaN,...
                    'vk', NaN,...
                     'u', NaN.*zeros(3,1),...
                     'n', NaN.*zeros(3,1),...
                     'k', NaN.*zeros(3,1),...
                     'h', NaN,...
                  'eorb', NaN);
flag_isValid = 1;
if(isnan(mass)==1)
  flag_isValid = 0;
end
if(sum(isnan(r0C0))~=0)
  flag_isValid = 0;
end
if(sum(isnan(v0C0))~=0)
  flag_isValid = 0;
end
if(sum(isnan(r0S0))~=0)
  flag_isValid = 0;
end   
if(sum(isnan(g0))~=0)
  flag_isValid = 0;
end

if(flag_isValid == 1)
  g          = norm(g0);
  capInfo.k  = -g0./g;

  v0C0k      = (v0C0'*capInfo.k)*capInfo.k;
  v0C0u      = v0C0 - (v0C0'*capInfo.k)*capInfo.k;
  v0C0uMag   = sqrt(v0C0u'*v0C0u);

  capInfo.u  = v0C0u./(max(v0C0uMag, 1e-9));

  capInfo.n  = getCrossProductMatrix(capInfo.k)*capInfo.u; 
  
  capInfo.vu = v0C0u'*capInfo.u;
  capInfo.vk = v0C0k'*capInfo.k;

  %Ground projection of the center of mass
  rSG0         = r0C0 - r0S0;
  capInfo.r0G0 = rSG0 - (rSG0'*capInfo.k)*capInfo.k;
  capInfo.h    = rSG0'*capInfo.k;

  capInfo.lu   = capInfo.vu*sqrt(capInfo.h/g);
  capInfo.r0F0 = capInfo.r0G0 + capInfo.u.*capInfo.lu;
  capInfo.eorb = 0.5*(capInfo.vu*capInfo.vu) ...
               - 0.5*(g/capInfo.h)*(capInfo.lu*capInfo.lu);

end           


