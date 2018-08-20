% Convergence rates of DG approximations with respect to the number of
% planewaves for different eigenvalues
% Apr 19th, 2018

[lambda,~] = eigen_dg_2atom_multieigenvalues(10,1,5,5,5,5,5,10000);
stand = lambda;
for i = 1:8
    x(i) =0.5* i;
    [lambda,~] =eigen_dg_2atom_multieigenvalues(10,1,0.5*i,5,5,5,5,10000);
    error_dg1(i) = abs(lambda(1,1)-stand(1,1));
    error_dg2(i) = abs(lambda(2,2)-stand(2,2));
    error_dg3(i) = abs(lambda(3,3)-stand(3,3));
end
semilogy(x,error_dg1,'ro-',x,error_dg2,'gs-',x,error_dg3,'bd-');

hold on
clear x
x=0:5/100:5;
y=2*exp(-x);
semilogy(x,y,'k--');

xlabel('truncation of plane wave K');
ylabel('numerical errors of eigenvalues');
axis([0.3,4.5,10^(-3),10]);
legend('1st eigenvalue','2nd smallest eigenvalue','3rd smallest eigenvalue','y=2e^{-x}')
