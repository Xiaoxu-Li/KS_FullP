% strategy A
% predict E_opt by linear regression
function [Eopt1] = strategyA(Ecv,indiv,eps)
% Ecv Ecut vector;
% indiv indicator vector;
if indiv(end) > indiv(end-1) || indiv(end) == indiv(end-1)
   Eopt1 = Ecv(end) + 0.5;
   return
end
b = polyfit(sqrt(Ecv),log(indiv),1);
slope = b(1);
cept = b(2);
syms x
Eopt1 = eval(solve(eps-exp(slope*sqrt(x)+cept),x));
return