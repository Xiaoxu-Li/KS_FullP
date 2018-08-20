% integration y=\int_0^R r^2*chi_i(r)*chi_j(r) dr
function y = int_overlap_radius_slater(i,j,R,n_integral)
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (xl^2*r_basis_slater(i,xl,R)*r_basis_slater(j,xl,R) +xr^2* r_basis_slater(i,xr,R)*r_basis_slater(j,xr,R));
end
y=s;