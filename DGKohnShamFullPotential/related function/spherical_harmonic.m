function [real, imag] = spherical_harmonic(l,m,theta,phi)
% return Y_lm(\theta,\phi)
% a=[ 1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800 87178291200 1307674368000 20922789888000 355687428096000];

 pre_part = sqrt(factorial(l-m)/factorial(l+m));
%pre_part = sqrt(a(l-m+1)/a(l+m+1));
if m >= 0
    legendre_part = legendre(l,cos(theta));
    real = pre_part * legendre_part(m+1) * cos(m*phi);
    imag = pre_part * legendre_part(m+1) * sin(m*phi);
else
     legendre_part = (-1)^m * factorial(l+m)/factorial(l-m) * legendre(l,cos(theta));
    %legendre_part = (-1)^m * a(l+m+1)/a(l-m+1) * legendre(l,cos(theta));
    real = pre_part * legendre_part(-m+1) * cos(m*phi);
    imag = pre_part * legendre_part(-m+1) * sin(m*phi);
end

 real = sqrt((2*l+1)/(4*pi)) * real;
 imag = sqrt((2*l+1)/(4*pi)) * imag;