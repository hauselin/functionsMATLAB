function [cl] = cbsymmetric(x)
colorbar
if nargin == 1
    caxis([-x x]);
else
    cl = [-max(abs(caxis)) max(abs(caxis))];
    caxis(cl);
end
end

