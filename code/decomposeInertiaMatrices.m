function [eigenValues, eigenVectors] = ...
  decomposeInertiaMatrices( vectorOfRowWiseInertiaMatrices)


n = size(vectorOfRowWiseInertiaMatrices,1);
eigenValues         = zeros(n,3);
eigenVectors        = zeros(n,9);

for i=1:1:n

  jrow = vectorOfRowWiseInertiaMatrices(i,:);
  Jcmi = [jrow(1,1),jrow(1,2),jrow(1,3);...
          jrow(1,4),jrow(1,5),jrow(1,6);...
          jrow(1,7),jrow(1,8),jrow(1,9)];  

  [v,d] = eig(Jcmi);
 
  eigenValues(i,:)        = [d(1,1),d(2,2),d(3,3)];
  eigenVectors(i,:)       = [v(1,:),v(2,:),v(3,:)];

  assert( d(1,1) > 0 && d(2,2) > 0 && d(3,3) > 0,...  
          ['Error: whole body inertia matrix has a negative eigen value,'...
          'which is physically impossible. ']);
end