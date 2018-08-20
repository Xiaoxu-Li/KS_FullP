function y = int_vext_sink(k,R,n_integral)
% integration y=1/k*\int_0^R r*v(r)*sin(kr) dr
% using the potential V_smooth

N = n_integral;   % integration discretization
h = R/N;
s = 0;
if k == 0
    for i=1:N
        xl = (i-1)*h;
        xr = i*h;
        s = s + 0.5 * h * (xl*V_smooth_r(xl,4.5) + xr*V_smooth_r(xr,4.5));
    end
else
    for i=1:N
        xl = (i-1)*h;
        xr = i*h;
        s = s + 0.5 * h * (V_smooth_r(xl,4.5)*sin(k*xl) + V_smooth_r(xr,4.5)*sin(k*xr));
    end
    s = s/k;
end
y=s;
return;