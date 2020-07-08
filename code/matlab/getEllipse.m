function ellipseXY = getEllipse(centerXY, radiiXY, npts)

ellipseXY = zeros(npts,2);

for i=1:1:npts
  th = 2*pi*((i-1)/(npts-1));
  ellipseXY(i,:) = centerXY+[radiiXY(1,1)*cos(th),radiiXY(1,2)*sin(th)];
end


