% Convergence rates of DG approximations with respect to the number of
% planewaves (different R)
% Apr 7th, 2018

[lambda,~] = eigen_dg_per_new(10,3,5,5,5,10000);
stand = lambda(1);
for i = 1:8
    x(i) =0.5* i;
    [lambda,~] =eigen_dg_per_new(10,3,0.5*i,5,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
semilogy(x,error_dg,'ro-');

hold on
clear x
[lambda,~] = eigen_dg_per_new(10,1,5,5,5,10000);
stand = lambda(1);
for i = 1:8
    x(i) =0.5* i;
    [lambda,~] =eigen_dg_per_new(10,1,0.5*i,5,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
semilogy(x,error_dg,'gs-');

hold on
clear x
x=0:8/100:8;
y=exp(-x);
semilogy(x,y,'b:');

xlabel('truncation of plane wave K');
ylabel('numerical errors of eigenvalue');
legend('R=3','R=1','y=e^{-x}')