function y = r_basis_slater(i,r,R)
    y = 2^i*sqrt(2/factorial(2*i)) * r^(i-1) * exp(-r);
return;
