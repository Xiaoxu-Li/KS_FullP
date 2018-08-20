% Compare numerical errors of planewaves methods and DG methods using
% periodic potential
% Dec 6th, 2017

[lambda,~] = eigen_dg_per_new(10,3,5,5,5,10000);
% eval(['save C:\Users\lixiaoxu\Documents\MINE\MATLAB\data','\M',num2str(j),' lambda'])
stand = lambda(1);
for i = 1:8
    x(i) =0.5* i;
    [lambda,~] =eigen_dg_per_new(10,3,0.5*i,5,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
loglog(x,error_dg,'ro-');

hold on
clear x
[lambda,phi] = eigen_solve_per(10, 8);
stand = lambda(1);
for i = 1:15
    x(i) = 0.5*i;
    [lambda,phi] = eigen_solve_per(10, i*0.5);
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
legend('DG method(R=3)','planewave method','y=3x^{-3}')
axis([0.3,10,10^(-4),10^2])