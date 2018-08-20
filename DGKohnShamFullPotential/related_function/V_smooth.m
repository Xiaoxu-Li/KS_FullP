function y = V_smooth(x,R)
% return the function y = -2/R            r>=R
%                         -1/r-1/(2R-r)            r<R
% which is a smooth function away from the singularity

if abs(x) >= R
    y = -2.0/R;    
else
    y = -1.0/abs(x) - 1.0/(2*R-abs(x));
end

return;