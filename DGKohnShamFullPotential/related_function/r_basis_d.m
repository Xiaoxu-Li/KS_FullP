% derivative of radial basis functions \chi'_i
function y = r_basis_d(i,r,R)
if i == 1
    y = 0;
else
    y = (i-1)*r^(i-2) *i/(R^i);  % normalization
%    y = (i-1)*r^(i-2) ;
end
return;