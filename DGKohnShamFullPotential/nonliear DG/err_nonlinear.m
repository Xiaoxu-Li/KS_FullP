% compute numerical errors of energy with respect to plane wave K 
% Aug 20th, 2018

% eval(['save C:\Users\lixiaoxu\Documents\MINE\MATLAB\data','\M',num2str(j),' lambda'])

[energy,~,~]=eigen_dg_nonlinear(10,3,5,4,3,10000);
stand=energy;
for i = 1:5
    x(i) =0.5* (i+3);
    [energy,~,~] =eigen_dg_nonlinear(10,1,0.5*(i+3),4,3,10000);
    error_energy(i) = abs(energy-stand);
end
plot(x,error_energy,'ro-');

% hold on
% clear x
% [lambda,phi] = eigen_solve_nonlinear(10, 8);
% stand = lambda(1);
% for i = 1:15
%     x(i) = 0.5*i;
%     [lambda,phi] = eigen_solve_per(10, i*0.5);
%     error_pw(i) = abs(lambda(1)-stand);
% end
% loglog(x,error_pw,'gx-');
% 
% hold on
% clear x
% x=0:8/100:8;
% y=3*x.^(-3);
% loglog(x,y,'b:');

xlabel('truncation of plane wave K');
ylabel('numerical errors of energy');
% legend('DG method(R=3)','planewave method','y=3x^{-3}')
% axis([0.3,10,10^(-4),10^2])