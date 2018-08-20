% Compare numerical errors of planewaves methods and DG methods using
% periodic potential in two-atom system
% Apr 16th, 2018

% [lambda,~] = eigen_dg_2atom_new(10,1,5,5,5,5,5,10000);
% stand = lambda(1);
stand=-0.2721;          % eigen_dg_2atom_new(10,1,5,5,5,5,5,10000);
for i = 1:9
    x(i) =0.5* i;
    [lambda,~] =eigen_dg_2atom_new(10,1,0.5*i,5,5,5,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
loglog(x,error_dg,'ro-');

hold on
clear x
[lambda,~] = eigen_solve_2atom(10, 8);
stand = lambda(1);
for i = 1:14
    x(i) = 0.5*i;
    [lambda,~] = eigen_solve_2atom(10, i*0.5);
    error_pw(i) = abs(lambda(1)-stand);
end
loglog(x,error_pw,'gx-');

hold on
clear x
x=0:8/100:8;
y=3*x.^(-3);
loglog(x,y,'b:');

xlabel('truncation of plane wave K');
ylabel('numerical errors of eigenvalue');
legend('2 atom (R=1)','planewave method','y=3x^{-3}')
axis([0.3,10,10^(-4),10^2])