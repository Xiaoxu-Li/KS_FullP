% Numerical errors for different types of radial basis functions
% Apr 7th, 2018

[lambda,~] = eigen_dg_per_new(10,3,5,6,5,10000);
stand = lambda(1);
for i = 1:4
    x(i) = i;
    [lambda,~] =eigen_dg_per_new(10,3,5,i,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
semilogy(x,error_dg,'ro-');

hold on
clear x
x=0:5/100:5;
y=2*exp(-1.5*x);
semilogy(x,y,'k--');

hold on
clear x
[lambda,~] = eigen_dg_per_slater_new(10,3,5,6,5,10000);
stand = lambda(1);
for i = 1:4
    x(i) = i;
    [lambda,~] =eigen_dg_per_slater_new(10,3,5,i,5,10000);
    error_slater(i) = abs(lambda(1)-stand);
end
semilogy(x,error_slater,'bd-');

xlabel('degrees of radial basis');
ylabel('numerical errors of eigenvalue');
legend('polynomials','y=2e^{-1.5x}','slater type orbitals')
