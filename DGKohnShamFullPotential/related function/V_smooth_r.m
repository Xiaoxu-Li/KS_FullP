function y = V_smooth_r(r,R)
% return the function y = -2*r/R            r<R
%                         -1-r/(2R-r)   r>R

if r >= R
    y = -2.0*r/R;    
else
    y = -1.0 - r/(2*R-r);
end

return;