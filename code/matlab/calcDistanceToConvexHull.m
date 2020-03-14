function [distanceToConvexHull, normalConvexHull, closestPointOnConvexHull]...
        = calcDistanceToConvexHull(r0P0, r0Q0V)
%%
% 
% @param r0P0 : 1 x 3 vector of the point of interest
% @param r0Q0V: n x 3 matrix of the points that define the closed convex
%               hull (beginning and end points must be the same, and the
%               points must be ordered so that the line between the
%               indices (i+1) and i form a line segment that is a boundary
%               of the convex hull.
% @return distanceToConvexHul
%         the distance between r0P0 and the nearest edge. Positive numbers
%         indicates that r0P0 is outside of the convex hull while
%         negative numbers mean r0P0 is inside the convex hull.
%
%         normalConvexHull 
%         the direction of the convex hull line segment normal that
%         intersects the point r0P0
%
%         closestPointOnConvexHull
%         the point on the convex hull that is closest to the point
%         given
%%
distanceToConvexHull      = NaN;
normalConvexHull          = zeros(size(r0P0)).*NaN; 
closestPointOnConvexHull  = zeros(size(r0P0)).*NaN;

if((sum(isnan(r0P0))==0) && (sum(sum(isnan(r0Q0V)))==0) )
  assert(length(r0P0) ==2);
  assert(size(r0Q0V,2)==2);
  
  %Center the points
  r0C0   = [mean(r0Q0V(:,1)),mean(r0Q0V(:,2))];
  rCP0   = r0P0  - r0C0;  
  distanceToConvexHull = Inf;
 
  rC10 = zeros(size(rCP0));
  rC20 = zeros(size(rCP0));
  r120 = zeros(size(rCP0));
  e12  = zeros(size(rCP0));
  n12  = zeros(size(rCP0));
  rCX0 = zeros(size(rCP0));
  rXP0 = zeros(size(rCP0));
  
  for i=2:1:size(r0Q0V,1)
  
    rC10 = r0Q0V(i-1,:) - r0C0;
    rC20 = r0Q0V(i,:)   - r0C0;
    r1P0 = rCP0-rC10;
    
    %Direction vector along the line segment
    r120 = rC20-rC10;
    e12 = r120 ./ norm(r120);
   
    %Direction vector pointing outside of the convex hull
    n12  = rC10 - (rC10*e12').*e12;
    n12  = n12 ./ norm(n12);
   
    A = [e12(1,1), n12(1,1);...
         e12(1,2), n12(1,2)];

    b = [r1P0(1,1);...
         r1P0(1,2)];
       
    x = A\b;
    
    %Check if the point of closest approach is within the line segment
    if( x(1) <= norm(r120) && x(1) >= 0)
      if( abs(x(2)) < abs(distanceToConvexHull))
        distanceToConvexHull      = x(2,1);
        normalConvexHull          = n12;        
        closestPointOnConvexHull  = r0C0 + rC10 + e12.*x(1,1);
        here=1;
      end
    end
    
  end
  

  
  flag_debug = 0;
  if(flag_debug ==1)
%     if(exist('figDebugClosestDistanceConvexHull')==0)
%       figDebugClosestDistanceConvexHull = figure;
%     else
%       clf(figDebugClosestDistanceConvexHull);
%     end
    plot(r0Q0V(:,1), r0Q0V(:,2),'r');
    hold on;
    plot(r0P0(1,1), r0P0(1,2),'x','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(closestPointOnConvexHull(1,1), closestPointOnConvexHull(1,2),...
          'ob','MarkerSize',10);
    hold on;
    
    plot([closestPointOnConvexHull(1,1); ...
            closestPointOnConvexHull(1,1)+normalConvexHull(1,1).*0.1],...
         [closestPointOnConvexHull(1,2); ...
            closestPointOnConvexHull(1,2)+normalConvexHull(1,2).*0.1],'--b')
    hold on;
    axis equal;
    here=1;    
  end
  
end