% integration y=\int_0^R chi_i(r)*chi_j(r)*chi_n(r) dr
function y = int_overlap_3(i,j,n,R,n_integral)
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (r_basis(i,xl,R)*r_basis(j,xl,R) *r_basis(n,xl,R)+ r_basis(i,xr,R)*r_basis(j,xr,R)*r_basis(n,xr,R) );
end
y=s;
return;
