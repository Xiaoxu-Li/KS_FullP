% integration y=1/k^3*\int_0^R r*sin(kr)*chi_i(r)*chi_j(r) dr
function y=int_vext_per_sin(K,i,j,R,n_integral)
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (xl*sin(K*xl)*r_basis(i,xl,R)*r_basis(j,xl,R) + xr*sin(K*xr)*r_basis(i,xr,R)*r_basis(j,xr,R));
end
y=s/K^3 ;
return;  