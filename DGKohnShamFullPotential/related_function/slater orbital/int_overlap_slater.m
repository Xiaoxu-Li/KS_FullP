% integration y=\int_0^R chi_i(r)*chi_j(r) dr
function y = int_overlap_slater(i,j,R,n_integral)
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (r_basis_slater(i,xl,R)*r_basis_slater(j,xl,R) + r_basis_slater(i,xr,R)*r_basis_slater(j,xr,R));
end
y=s;
return;
