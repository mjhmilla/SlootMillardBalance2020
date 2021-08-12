function rPXP = transformXY(r0X0, r0P0, EP0)
  rPXP = zeros(size(r0X0));
  for i=1:1:size(r0X0,1)
    rPXP(i,:) = r0X0(i,:)-r0P0;
    rPXP(i,:) = ( EP0 * rPXP(i,:)' )';
  end