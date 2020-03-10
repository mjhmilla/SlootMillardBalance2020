function distanceToConvexHull = calcDistanceToConvexHull(r0P0, convexHull)

distanceToConvexHull = NaN;

if(sum(isnan(r0P0))==0)
  convexHullCenter = mean(convexHull);
  assert( abs(convexHullCenter(1,1)) < 1e-6 ...
       && abs(convexHullCenter(1,2)) < 1e-6);
     
  eP0   = r0P0 ./ norm(r0P0);
  
  i = 2;
  flag_found = 0;

  r010 = zeros(size(r0P0));
  r020 = zeros(size(r0P0));
  e10 = zeros(size(r0P0));
  e20 = zeros(size(r0P0));
  d1P = 0;
  d2P = 0;
  d12 = 0;
  
  while i <= size(convexHull,1) && flag_found == 0

    r010 =   convexHull(i,:);
    e10 = r010./norm(r010);
    
    r020 = convexHull(i-1,:);
    e20 = r020./norm(r020);
    
    d1P = eP0*e10';
    d2P = eP0*e20';
    d12 = e20*e10';
    
    if(d1P >= d12 && d2P >= d12 && ( (d1P > 0) || (d2P > 0)))
        flag_found = 1;
    end
      
    i=i+1;
  end
  if(flag_found ==0)
    here=1;
  end
  
  assert(flag_found ~= 0);

  %w1 = d1P./(d1P + d2P);
  %w2 = d2P./(d1P + d2P);
  %r0Q0 = w1.*r010 + w2.*r020;

  r120 = r020-r010;
  e120 = r120./norm(r120);
  
  A = [eP0(1,1), -e120(1,1);...
       eP0(1,2), -e120(1,2)];
  
  b = [r010(1,1);...
       r010(1,2)];
     
  x = A\b;
  r0Q0 = eP0.*x(1);
     
  distanceToConvexHull = (r0P0-r0Q0)*eP0';
  
  flag_debug = 0;
  if(flag_debug ==1)
    fig = figure;
    plot(convexHull(:,1), convexHull(:,2),'r');
    hold on;
    plot([0;r0P0(1,1)], [0;r0P0(1,2)],'k','LineWidth',2);
    hold on;
    plot(r0P0(1,1), r0P0(1,2),'x','LineWidth',2,'MarkerSize',10);
    hold on;
    plot([0;r0Q0(1,1)], [0;r0Q0(1,2)],'g');
    hold on;
    plot(r0Q0(1,1), r0Q0(1,2),'go','MarkerSize',10);
    hold on;
    axis equal;
    
    
    
  end
  
end