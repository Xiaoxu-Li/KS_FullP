% Compare numerical errors of planewaves methods and DG methods
% using V_smooth
% Dec 14th, 2017

[lambda,~] = eigen_dg(10,1,5,5,5,10000);
stand = lambda(1);
for i = 1:8
    x(i) =0.5* i;
    [lambda,~] =eigen_dg(10,1,0.5*i,5,5,10000);
    error_dg(i) = abs(lambda(1)-stand);
end
loglog(x,error_dg,'ro-');

hold on
clear x
[lambda,phi] = eigen_solve(10, 8);
stand = lambda(1);
for i = 1:15
    x(i) = 0.5*i;
    [lambda,phi] = eigen_solve(10, i*0.5);
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
legend('DG method(R=1)','planewave method','y=3x^{-3}')
axis([0.3,10,10^(-4),10^2])
% title(' Compare numerical errors of planewaves methods and DG methods')