function pt = calcRayToConvexHullIntersectionPoint(angle, angleConvHull,radiusConvHull)

assert(angle >= -pi);
assert(angle <=  pi);

assert( min(angleConvHull) > -pi);
assert( max(angleConvHull) <  pi);


pt = [0,0];

angle1 = 0;
angle2 = 0;
angle3 = 0;
l1 = 0;
l2 = 0;
l3 = 0;

x1 = 0;
y1 = 0;
x2 = 0;
y2 = 0;

if(angle <angleConvHull(1,1) && angle >= -pi)
  angle1 = angleConvHull(end,1)-2*pi;
  angle2 = angleConvHull(1,1);
  l1 = radiusConvHull(end,1);
  l2 = radiusConvHull(1,1);
  
  x1 = l1*cos(angle1);
  y1 = l1*sin(angle1);
  x2 = l2*cos(angle2);
  y2 = l2*sin(angle2);  
  
  l3 = sqrt((x2-x1)^2 +(y2-y1)^2);
  angle3 = atan2(y2-y1,x2-x1);
  
elseif(angle > angleConvHull(end,1) && angle <= pi)

  angle1 = angleConvHull(end,1);
  angle2 = angleConvHull(1,1)+2*pi;
  l1 = radiusConvHull(end,1);
  l2 = radiusConvHull(1,1);
  
  x1 = l1*cos(angle1);
  y1 = l1*sin(angle1);
  x2 = l2*cos(angle2);
  y2 = l2*sin(angle2);  
  
  l3 = sqrt((x2-x1)^2 +(y2-y1)^2);
  angle3 = atan2(y2-y1,x2-x1);
    
else
  flag_found = 0;
  idx=2;
  while(flag_found==0 && idx <= size(angleConvHull,1))

    angle2 = angleConvHull(idx);% atan2(convHullXY(idx,2) , convHullXY(idx,1));
    angle1 = angleConvHull(idx-1);%atan2(convHullXY(idx-1,2),convHullXY(idx-1,1));

    if( angle >= angle1 && angle <= angle2 )

      l1 = radiusConvHull(idx-1,1);
      l2 = radiusConvHull(idx,1);
      
      x1 = l1*cos(angle1);
      y1 = l1*sin(angle1);
      x2 = l2*cos(angle2);
      y2 = l2*sin(angle2);       
      
      dx = x2-x1;
      dy = y2-y1;    
      l3 = sqrt(dx*dx+dy*dy);
      
      angle3 = atan2(dy,dx);

      flag_found = 1;
    end
    idx=idx+1;  
  end
  assert(flag_found==1);
end


da1 = angle-angle1;
da2 = angle2-angle;
dx = x2-x1;
dy = y2-y1;

if(da1>da2)
  h = l1*sin(angle-angle1);
  da3 = angle3-angle;
  d = h/sin(da3);
  normDist = d/l3;
  pt = [x1,y1] + [dx,dy].*(normDist);
else
  h = l2*sin(angle2-angle);
  da3 = angle3-angle;
  d = h/sin(da3);
  if(d<=0)
    here=1;
  end
  assert(d>0);
  normDist = d/l3;
  pt = [x2,y2] - [dx,dy].*(normDist);
end



