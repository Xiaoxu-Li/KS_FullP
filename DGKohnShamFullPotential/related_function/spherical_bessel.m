% spherical Bellsel functions
function y = spherical_bessel(n,z)
y = besselj(n+1/2,z)*sqrt(pi/(2*z));
return;
