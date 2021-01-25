function fpeInfo = calc3DFootPlacementEstimatorInfo(m,...
                                          r0C0,...
                                          v0C0,...                                                    
                                          JC0,...                                                    
                                          HC0,...
                                          r0S0,...
                                          g0,...
                                          omegaSmall,...
                                          numericTolerance,...
                                          maximumIterations,...
                                          flag_evaluateDerivatives,...
                                          flag_verbose,...
                                          fpeInfoGuess)
%%
% This function computes the location of a balance-restoring step on a horizontal
% surface using the 3DFPE algorithm described by Millard et al. [1]. All 
% quantities are expressed in MKS units. 
%
%
% [1] Millard, M., McPhee, J., & Kubica, E. (2012). Foot placement and balance 
%  in 3D. Journal of Computational and Nonlinear Dynamics, 7(2), 021015.
%
% A small apology for the terse variable names. These terse 
% names are needed to do write these functions without making obvious errors.
% Nice human-readable variable names could have been used, however, time and 
% memory would have to have been spent to copy these over. In case you need to 
% modify the guts of this code a description of the variable convention appears
% below.
%
% The variables follow a notation convention that is similar to that used in 
% multibody-dynamics. 
%                                              
%  Global quantities           : [ Descriptor ] [Resolved in frame]
%   e.g. g0 : gravity vector in K0, the inertial frame
%
%  Single-point/body quantites : [ Descriptor ] [About Point      ] [Resolved in frame] [*optional* operator]
%  Two-point quantites         : [ Descriptor ] [From Point       ] [To Point         ] [Resolved in frame] [*optional* operator]
%
% Global descriptors
%   g : gravity
%
% Single Point/Body Descriptors 

%   e  : unit vector
%   rm : rotation matrix
%   K  : a frame: consists of a point in space and a rotation matrix
%   J  : inertia
%   H  : angular momentum
%
%   (specific to this code)
%   u : vector in the direction of a balancing step (See Fig. 6 of Millard et al.)
%   v : vertical vector
%   n : vector normal to the the plane
%   P : the projection frame: centered at the CoM ground projection and 
%           with a u, k direction vectors that are perpendicular to the 
%           horizontal component of the angular momentum vector (when taken)
%           about the CoM ground projection.
%
% Two-Point Descriptors:
%   r: position
%   v: velocity
%   w: angular velocity
%
% Point
%   C: the whole-body-center of mass
%   0: the origin of the lab/inertial frame
%   S: the contact surface frame origin
%   P: the projection plane origin
%
% Frames
%   0: the lab/inertial frame
%   S: the contact surface frame origin
%
%
% (optional) Operator
%   x: Cross-product matrix. For example r0C0 is the vector, r0C0x is the
%      cross-product matrix of r0C0.
%   p: projection
%
% Examples:
%   r0C0 : (r) position vector 
%          (0) from the lab origin 
%          (C) to the center-of-mass (CoM) 
%          (0) resolved in the coordinates of the lab frame
%
%   rmS0 : (rm) rotation matrix
%          (S) from the surface frame 
%          (0) to the coordinates of the root frame

%   JC0 :  (J) inertia matrix
%          (C) about the center-of-mass (CoM) 
%          (0) resolved in the coordinates of the lab frame
%
%
% @param mass: A scalar of the mass of the whole body
%
% @param r0C0: (com position) An 3-by-1 vector of center of mass CoM location
%
% @param v0C0: (com velocity) An 3-by-1 vector of CoM velocity 
%
% @params JC0 (whole body inertia): A 3-by-3 matrix containing the inertia 
%           matrix. The inertia matrix is taken about the current CoM location 
%           and resolved in the coordinates inertial (or lab) frame.
%
% @params HC0 (whole body angular momentum): An 3-by-1 vector containing the 
%           whole-body-angular momentum about the CoM location resolved in
%           the coordinates of the inertial (or lab) frame
% 
% @param  r0S0 (contact surface position): An 3-by-1 vector containing one 
%               position of the surface that the object will be contacting in 
%               order to balance.
%
% @param g0: a 3-by-1 vector of gravity resolved in the inertial frame
%
% @param numericTolerance: a small but non-zero number used to define the
%                          tolerance to which Eqn. 45 of Millard et al. should
%                          be satisfied.
%
% @param maximumIterations: the maximum number of iterations the root-solving
%                           routine is permitted to take to solve for the
%                           3DFPE location within the desired numeric tolerance.
%
% @param flag_verbose: setting this to 1 will print a great deal of debugging
%                      information to the terminal.
%
% @returns fpeInfo a structure contain data on the location of the 3DFPE, its
%             partial derivatives, and the error (both numerical and physical)
%             in its calculation. These quantites are quite specific, and so
%             where possible the descriptions below will mention the equation
%             number in Millard et al. that correspond to the quantity.
%
%  .fpe               : The location of the balancing point on the contact surface (See Fig. 6)
%  .phi               : The foot contact angle: the angle between the gravity vector and the vector between the CoM and the contact location (See Fig. 6)
%
%  .f                 : Numerical 3DFPE constraint error (Eqn. 45)
%  .projectionError   : Percentage of the body's momentum that is not in the projection plane (Eqn. 44)
%  .n                 : Vector that is normal to the projection plane (Eqn. 43)
%  .r0G0              : CoM ground projection location
%  .Df_Dphi           : Partial derivative of f w.r.t the contact angle
%  .Df_Dw0C0n         : " " of f w.r.t. the whole body angular velocity
%  .Df_Dh             : " " of f w.r.t. the contact height (See Fig. 6)
%  .Df_Dv0C0u            : " " of f w.r.t. the velocity of the CoM parallel to the surface
%  .Df_Dv0C0k            : " " of f w.r.t. the velocity of the CoM normal to the surface
%  .Df_DJ             : " " of f w.r.t. the inertia of the body
%  .Df_Dm             : " " of f w.r.t. the mass of the body
%  .Df_Dg             : " " of f w.r.t. gravity
%
% ** FEATURE NOT YET IMPLEMENTED **
% It is possible to redefine this method to find steps on an angled surface, 
% or even a curved surface. This feature has not been implemented yet as it
% requires re-deriving the constraint equation (Eqn. 45 in Millard et al.) 
% that must be satisfied. To implement this feature in
% the most generic manner:
%
%   0. Define two constraint equations in two unknow0C0ns: L and phi.
%   1. Eqn 1: Change Eqn. 45 so that L is left as L: do not substitute 
%      L/cos(phi)
%   2. Eqn 2: The implicit surface function evaluated at the location of the
%             foot.
%   3. Solve this system using a Newton method. Be careful to limit numerical
%      blunders by limiting the norm of the largest step size and rescale things
%      appropriately for the surface.
%
% @param  rmS0: ** FEATURE NOT YET IMPLEMENTED **
%               3-by-3 rotation matrix that rotates vectors from the S frame
%               into the inertial frame. This variable is required when the 
%               contacting surface is not perpendicular to the gravity vector
%
% @param enS0: ** FEATURE NOT YET IMPLEMENTED **
%             3-by-1 direction vector that defines the surface normal.
%             This variable is required when the contacting surface is not 
%             perpendicular to the gravity vector
%%%


fpeInfo = struct('f'               ,  zeros(1,1).*NaN,...
                 'phi'             ,  zeros(1,1).*NaN,...
                 'r0F0'            ,  zeros(3,1).*NaN,...
                 'v0F0'            ,  zeros(3,1).*NaN,... 
                 'projectionError' ,  zeros(1,1).*NaN,...
                 'n'               ,  zeros(3,1).*NaN,...
                 'u'               ,  zeros(3,1).*NaN,...
                 'k'               ,  zeros(3,1).*NaN,...                 
                 'r0G0'            ,  zeros(3,1).*NaN,...
                 'w0C0n'           ,  zeros(1,1).*NaN,...
                 'w0C0'            ,  zeros(3,1).*NaN,...
                 'w0G0'            ,  zeros(3,1).*NaN,...                 
                 'h'               ,  zeros(1,1).*NaN,...
                 'v0C0u'           ,  zeros(1,1).*NaN,...
                 'v0C0k'           ,  zeros(1,1).*NaN,...
                 'v0G0'            ,  zeros(3,1).*NaN,...
                 'J'               ,  zeros(1,1).*NaN,...                 
                 'Df_Dphi'         ,  zeros(1,1).*NaN,...
                 'Df_Dw0C0n'       ,  zeros(1,1).*NaN,...
                 'Df_Dh'           ,  zeros(1,1).*NaN,...
                 'Df_Dv0C0u'       ,  zeros(1,1).*NaN,...
                 'Df_Dv0C0k'       ,  zeros(1,1).*NaN,...
                 'Df_DJ'           ,  zeros(1,1).*NaN,...
                 'Df_Dm'           ,  zeros(1,1).*NaN,...
                 'Df_Dg'           ,  zeros(1,1).*NaN,...
                   'Ds_Dl'           ,  zeros(1,1).*NaN,...
                   'Ds_DJ'           ,  zeros(1,1).*NaN,...
                   'Ds_DE'           ,  zeros(1,1).*NaN,... 
                   'Ds_Dv0C0u'       ,  zeros(1,1).*NaN,...
                   'Ds_Dv0C0k'       ,  zeros(1,1).*NaN,...
                   'Ds_Dw0C0n'       ,  zeros(1,1).*NaN,...   
                 'Dphi_Dl'           ,  zeros(1,1).*NaN,...
                 'Dphi_DJ'           ,  zeros(1,1).*NaN,...
                 'Dphi_DE'           ,  zeros(1,1).*NaN,... 
                 'Dphi_Dv0C0u'       ,  zeros(1,1).*NaN,...
                 'Dphi_Dv0C0k'       ,  zeros(1,1).*NaN,...
                 'Dphi_Dw0C0n'       ,  zeros(1,1).*NaN,...                      
                   'l'               ,  zeros(1,1).*NaN,...
                   'E'               ,  zeros(1,1).*NaN);


%===============================================================================
%
% Check if the inputs are valid
%
%===============================================================================
flag_validInputs=1;
if(sum(isnan(m)) > 0 || m <= 0)
  flag_validInputs =0;
end
if(sum(isnan(r0C0)) > 0)
  flag_validInputs =0;
end
if(sum(isnan(v0C0)) > 0)
  flag_validInputs =0;
end
if(sum(sum(isnan(JC0))) > 0)
  flag_validInputs =0;
end
if(sum(isnan(HC0)) > 0)
  flag_validInputs =0;
end
if(sum(isnan(r0S0)) > 0)
  flag_validInputs =0;
end
if(sum(isnan(g0)) > 0)
  flag_validInputs =0;
end
if(sum(isnan(numericTolerance)) > 0 || numericTolerance <= 0)
  flag_validInputs =0;
end
if(sum(isnan(maximumIterations)) > 0 || maximumIterations <= 0)
  flag_validInputs =0;
end

if(flag_evaluateDerivatives ~= 0 && flag_evaluateDerivatives ~= 1)
  flag_validInputs =0;
end
if(flag_verbose ~= 0 && flag_verbose ~= 1)
  flag_validInputs =0;
end


if(flag_validInputs==1)               
  %===============================================================================
  %1. Compute the projection plane
  %   a. Resolve the whole body angular momentum about the CoM ground projection
  %   b. Direction vectors ux and uz
  %   c. epsilon - the percentage of the body's Hg
  %===============================================================================


  %   a. Resolve the whole body angular momentum about the CoM ground projection

  %Gravity normal vector in the lab frame
  gNorm= norm(g0);
  eg0 = g0./gNorm;

  %Compute the whole-body average angular velocity
  w0C0 = JC0 \ HC0;

  fpeInfo.w0C0 = w0C0;
  %Get the CoM projection onto the contacting surface
  %rSCS =  rmS0' * (r0C0-r0S0);
  %egS  = (rmS0' * g0) ./ norm(g0);
  %rSGS =  rSCS - (rSCS' * egS).*egS;
  %rSG0 =  rmS0 * rSGS;
  %r0G0 =  r0S0 + rSG0;

  rSC0  = (r0C0-r0S0);
  rGC0  = (rSC0'*eg0)*eg0;
  r0G0  = r0C0 - rGC0;
  fpeInfo.r0G0 = r0G0;

  %Compute the whole body angular momentum vector about the COM ground projection
  rGC0x = getCrossProductMatrix(rGC0);
  JG0   =  JC0 - m.*rGC0x*rGC0x;
  
  %Checking the generalization of the parallel axis theorem -rGC0x*rGC0x
  %rr = (rGC0'*rGC0)*eye(3,3);
  %rrt= rGC0*rGC0';
  %pa = rr - rrt;
  
  HG0   = HC0 + rGC0x * (m.*v0C0);

  w0G0 = JG0 \ HG0;
  fpeInfo.w0G0 = w0G0;  
  
  fpeInfo.HG0 = HG0;


  %   b. Direction vectors ux and uz

  %The projection plane has a positive direction k that opposes gravity
  ev0   = -eg0;
  fpeInfo.k = ev0;
  %The projection plane has a normal vector in the direction of the horizontal
  %component of HG0.
  HG0k  = HG0'*eg0;
  HG0p  = HG0 - (HG0k).*eg0;    
  en0   = HG0p ./ ( max(numericTolerance,norm(HG0p)) );

  %The direction of the step lies along the plane found using the cross product
  %of the plane's vertical and normal vectors. Note that the stabilizing
  %step can be in the positive/negative direction of this direction vector.
  eu0   = getCrossProductMatrix(en0)*ev0;
  fpeInfo.u = eu0;

  omegaEps = [omegaSmall;omegaSmall;omegaSmall];
  HG0eps = JG0*omegaEps;
  
  fpeInfo.projectionError = norm(HG0k)/max( norm(HG0), norm(HG0eps));
  fpeInfo.n               = en0;


  %==========================================================================
  %2. Project the state of the body and its inertia onto
  %   the projection plane.
  %==========================================================================

  JC0n    = en0' * JC0 * en0;
  w0C0n  = en0' * w0C0;
  v0C0n  = [(eu0'*v0C0);(ev0'*v0C0)];




  %hV is going to have to re-calculated at every new phi to allow for
  %contact surfaces that are not normal to the gravity vector.

  %==========================================================================
  %3. Solve for phi the acute angle between the gravity direction vector and 
  %   the  diretion vector from the CoM to the foot placement location.
  %==========================================================================

  %The surface is not necessarily horizontal. To calculate the leg length
  %that will reach the surface we need to know the angle alpha the surface makes
  %with the horizontal vector eu0 in the direction of travel

  %etS0 : unit tangent direction vector on the surface S resolved in K0
  %etS0 = getCrossProductMatrix(en0)*enS0; 
  %etS0 = etS0 ./ (max(numericTolerance,norm(etS0)));
  %alpha: the angle between etS0 and the horizontal in the direction of travel
  %alpha = asin( etS0'*g0 );


  g       = gNorm;
  phi     = pi/4;  
  maxStep = pi/8;


  %m : already defined
  %g : gravity magnitude - already defined
  J     = JC0n;
  h     = sqrt(rGC0'*rGC0);
  v0C0u    = v0C0n(1,1);
  v0C0k    = v0C0n(2,1);


  fpeInfo.J = J;
  fpeInfo.h = h;
  fpeInfo.v0C0u = v0C0u;
  fpeInfo.v0C0k = v0C0k;
  fpeInfo.w0C0n = w0C0n;

  f     = numericTolerance*10; 
  iter  = 1;
 
  if(isempty(fpeInfoGuess)==1)
    cosphi    = 0;
    cos2phi   = 0;
    sinphi    = 0;
    h2        = 0;   

    %Get close to the root using the bisection method

    %Initial solution
    cosphi    = cos(phi);
    cos2phi   = cosphi*cosphi;
    sinphi    = sin(phi);
    h2        = h*h;      
    f      = (cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))^2/(cos2phi*J+h2*m)+2*(cosphi-1)*cosphi*g*h*m;

    phiBest = phi;
    errBest = abs(f);
    fBest = f;

    delta = phi*0.5;

    for i=1:1:5
      phiL = phiBest-delta;
      cosphi    = cos(phiL);
      cos2phi   = cosphi*cosphi;
      sinphi    = sin(phiL);
      h2        = h*h;      
      fL      = (cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))^2/(cos2phi*J+h2*m)+2*(cosphi-1)*cosphi*g*h*m;
      errL   = abs(fL);

      phiR = phiBest+delta;
      cosphi    = cos(phiR);
      cos2phi   = cosphi*cosphi;
      sinphi    = sin(phiR);
      h2        = h*h;      
      fR      = (cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))^2/(cos2phi*J+h2*m)+2*(cosphi-1)*cosphi*g*h*m;
      errR   = abs(fR);

      if(errL < errBest && errL <= errR)
        phiBest = phiL;
        errBest = errL;
        fBest = fL;
      end

      if(errR < errBest && errR < errL)
        phiBest = phiR;
        errBest = errR;    
        fBest = fR;
      end

      delta = delta*0.5;

    end

    f = fBest;
    phi=phiBest;
  else
    
    phi = fpeInfoGuess.phi;
    
  end
  
  

  %Polish off using Newton's method
  while( abs(f) > numericTolerance && iter < maximumIterations )

    cosphi    = cos(phi);
    cos2phi   = cosphi*cosphi;
    sinphi    = sin(phi);
    h2        = h*h;      

    f      = (cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))^2/(cos2phi*J+h2*m)+2*(cosphi-1)*cosphi*g*h*m;
    DfDphi = (2*cosphi*sinphi*J*(cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))^2)...
                /(cos2phi*J+h2*m)^2+(2*(cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))*((-2*cosphi*w0C0n*sinphi*J)-h*m*sinphi*(sinphi*v0C0k+cosphi*v0C0u)+cosphi*h*m*(cosphi*v0C0k-sinphi*v0C0u)))/(cos2phi*J+h2*m)-2*cosphi*g*h*m*sinphi-2*(cosphi-1)*g*h*m*sinphi;


    deltaPhi = -f/DfDphi;

    if(abs( deltaPhi ) > maxStep)
       deltaPhi = maxStep*sign(deltaPhi); 
    end
    if(isnan(deltaPhi) == 0 && abs(f) > numericTolerance)
      phi = phi + deltaPhi;
    end

   iter = iter+1;      
  end
  
  if(abs(f) >= numericTolerance || phi <= -numericTolerance)
    here=1;
  end
  assert(abs(f) <= numericTolerance);  
  assert(phi >= -numericTolerance);   %Else we found the (non-physical negative root)
  
  tanphi = tan(phi);
  rGF0    = h*tanphi*eu0;
  r0F0    = r0G0 + rGF0;


  
  fpeInfo.f  = f;
  fpeInfo.phi= phi;
  fpeInfo.r0F0 = r0F0;

  if(flag_evaluateDerivatives==1)
    fpeInfo.Df_Dphi    = (2*J*cosphi*sinphi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^2+(2*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)*((-h*m*sinphi*(sinphi*v0C0k+cosphi*v0C0u))+cosphi*h*m*(cosphi*v0C0k-sinphi*v0C0u)-2*J*cosphi*w0C0n*sinphi))/(h2*m+J*cos2phi)-2*cosphi*g*h*m*sinphi-2*(cosphi-1)*g*h*m*sinphi;
    fpeInfo.Df_Dw0C0n  = (2*J*cos2phi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi);
    fpeInfo.Df_Dh      = (-(2*h*m*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^2)+(2*cosphi*m*(sinphi*v0C0k+cosphi*v0C0u)*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi)+2*(cosphi-1)*cosphi*g*m;
                                
    fpeInfo.Df_Dv0C0u     = (2*cos2phi*h*m*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi);
    fpeInfo.Df_Dv0C0k     = (2*cosphi*h*m*sinphi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi);
    fpeInfo.Df_DJ      = (2*cos2phi*w0C0n*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi)-(cos2phi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^2;
    fpeInfo.Df_Dm      = (-(h2*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^2)+(2*cosphi*h*(sinphi*v0C0k+cosphi*v0C0u)*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi)+2*(cosphi-1)*cosphi*g*h;
    fpeInfo.Df_Dg      = 2*(cosphi-1)*cosphi*h*m;               
    
    Ds_Dphi = h*(1+tanphi*tanphi);
    
    
    
    %l = h/cos(phi)
    %h = l*cos(phi)
    %dh/dl = cos(phi)
    dh_dl = cosphi;
    % f = fo + (Df/Dphi)Dphi + (Df/DJ)*dJ. 
    % f = fo = 0
    % Dphi/DJ = -(Df/DJ)/(Df/Dphi)
    
    fpeInfo.Dphi_Dl = (-1/fpeInfo.Df_Dphi)*(fpeInfo.Df_Dh)*dh_dl;
    fpeInfo.Ds_Dl = (Ds_Dphi)*fpeInfo.Dphi_Dl;
    
    fpeInfo.Dphi_DJ = (-1/fpeInfo.Df_Dphi)*(fpeInfo.Df_DJ);
    fpeInfo.Ds_DJ = (Ds_Dphi)*fpeInfo.Dphi_DJ;

    omegaPlusTest = (cos2phi*w0C0n*J+cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u))/(cos2phi*J+h2*m);
    omegaPlus=(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)/(h2*m+J*cos2phi);
    DomegaPlusDphi=((-h*m*sinphi*(sinphi*v0C0k+cosphi*v0C0u))+cosphi*h*m*(cosphi*v0C0k-sinphi*v0C0u)-2*J*cosphi*w0C0n*sinphi)/(h2*m+J*cos2phi)+(2*J*cosphi*sinphi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n))/(h2*m+J*cos2phi)^2;
    tv=((1/2)*((h2*m)/cos2phi+J)*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^2+g*h*m;
    DtvDphi=(2*(1/2)*h2*m*sinphi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(cosphi^3*(h2*m+J*cos2phi)^2)+(4*(1/2)*J*cosphi*((h2*m)/cos2phi+J)*sinphi*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)^2)/(h2*m+J*cos2phi)^3+(2*(1/2)*((h2*m)/cos2phi+J)*(cosphi*h*m*(sinphi*v0C0k+cosphi*v0C0u)+J*cos2phi*w0C0n)*((-h*m*sinphi*(sinphi*v0C0k+cosphi*v0C0u))+cosphi*h*m*(cosphi*v0C0k-sinphi*v0C0u)-2*J*cosphi*w0C0n*sinphi))/(h2*m+J*cos2phi)^2;    
        
    
    v0G0 = omegaPlus*(h/cosphi).*(eu0+ev0+en0);
    
    Ea = m*g*h/cosphi;
    Eb = 0.5*(J + m*(h/cosphi)^2)*omegaPlus^2 +m*g*h;
    if(abs(Ea-Eb) > 1e-3)
      fprintf('%1.6f\n',abs(Ea-Eb));
    end
    
    %  l = h/cosphi;
    %
    Dl_Dphi = h*sinphi/(cos2phi);
    Dphi_Dl = (1/Dl_Dphi);
    %cosphi=h/l
    
    %E = mgl
    Dl_DE = 1/(m*g);    
    fpeInfo.Dphi_DE = (Dphi_Dl*Dl_DE);
    fpeInfo.Ds_DE = (Ds_Dphi)*(Dphi_Dl*Dl_DE);
    fpeInfo.E = Eb;
    
    %Angular momentum = r x mv + Jw
    %The sensitivity of the foot placement to the input momentum under
    %th assumption that m, r, and J are constant is just the sensitivity
    %of the foot placement location to the input velocities v0C0u, vz and w
    %
    fpeInfo.Dphi_Dv0C0u = (-1/fpeInfo.Df_Dphi)*fpeInfo.Df_Dv0C0u;
    fpeInfo.Ds_Dv0C0u = (Ds_Dphi)*fpeInfo.Dphi_Dv0C0u;
    
    fpeInfo.Dphi_Dv0C0k = (-1/fpeInfo.Df_Dphi)*fpeInfo.Df_Dv0C0k;
    fpeInfo.Ds_Dv0C0k = (Ds_Dphi)*fpeInfo.Dphi_Dv0C0k;
    
    fpeInfo.Dphi_Dw0C0n = (-1/fpeInfo.Df_Dphi)*fpeInfo.Df_Dw0C0n;
    fpeInfo.Ds_Dw0C0n = (Ds_Dphi)*fpeInfo.Dphi_Dw0C0n;
    here=1;
  end
end

