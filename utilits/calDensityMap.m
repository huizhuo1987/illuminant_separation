function [qq_density, density]= calDensityMap(qq, numBins)


xbins = numBins;
ybins = numBins;

qx = (qq(1, :) +1)./2;
qy = (qq(2, :) +1)./2;
% qx = (qq(1, :) - min(qq(1, :)))./(max(qq(1, :)) - min(qq(1, :)));
% qy = (qq(2, :) - min(qq(2, :)))./(max(qq(2, :)) - min(qq(2, :)));

xgrid = round(1 + qx * (xbins - 1));
ygrid = round(1 + qy * (ybins - 1));

idx = sub2ind([ybins xbins], ygrid, xgrid);

density = zeros(ybins, xbins);
idU = unique(idx);
qq_density = zeros(1, size(qq, 2));

for kk = 1:length(idU)
    tempIDD = find(idx == idU(kk));
    density(idU(kk)) = length(tempIDD);
    qq_density(:, tempIDD) = length(tempIDD);
end
