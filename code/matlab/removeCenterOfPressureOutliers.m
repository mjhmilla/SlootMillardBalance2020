function [idxMainCluster, idxExcluded] = removeCenterOfPressureOutliers(cop,acceptFraction)


minDistance = 0.001; %So there's no divide by zero

%Get the distance between all pairs
d=pdist2(cop,cop,'squaredeuclidean');

%Calculate a density for every point
rho = zeros(size(cop,1),1);
for i=1:1:size(cop,1)
  rho(i,1) = sum((d(i,:)+minDistance).^(-1));
end

[rhoSorted,indexRhoSorted]=sort(rho,'descend');

idxMax = round(max( (size(cop,1)*acceptFraction), 1 ));


idxMainCluster = sort(indexRhoSorted(1:idxMax));
idxExcluded = sort(indexRhoSorted((idxMax+1):end,1));