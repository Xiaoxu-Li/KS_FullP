function y = r_basis_d_slater(i,r,R)
if i == 1
    y = -2^i*sqrt(2/factorial(2*i)) * exp(-r);
else
    y = 2^i*sqrt(2/factorial(2*i)) * ((i-1)*r^(i-2) - r^(i-1)) * exp(-r);
end
return;