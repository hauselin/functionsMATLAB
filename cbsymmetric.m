function [cl] = cbsymmetric(x)  % symmetric colorbar
colorbar
if nargin == 1
    caxis([-x x]);
else
    cl = [-max(abs(caxis)) max(abs(caxis))];
    caxis(cl);
end
end

