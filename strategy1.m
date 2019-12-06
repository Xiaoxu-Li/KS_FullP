% strategy 1
% predict E_opt by linear regression
function [Eopt1] = strategy1(Ecv,indiv,eps)
% Ecv Ecut vector;
% indiv indicator vector;
b = polyfit(sqrt(Ecv),log(indiv),1);
slope = b(1);
cept = b(2);
syms x
Eopt1 = eval(solve(eps-exp(slope*sqrt(x)+cept),x));
return


