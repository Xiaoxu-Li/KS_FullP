% radial basis functions \chi_i
function y = r_basis(i,r,R)
y = r^(i-1) *i/(R^i);    % normalization
% y = r^(i-1) ;
return;