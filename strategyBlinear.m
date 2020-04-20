% strategy B
% predict E_opt by bi-section method
function [Eopt2] = strategyBlinear(L,Ec,Neig,l,phi,eps,bieps)
Eg = 4*Ec;
% Ecv Ecut vector;
El = Ec; 
Er = Eg;
Em = (Er + El)/2;
[err_post1,~] = residualinv(L, Ec, Neig, Eg, phi);
[err_post2,~] = residualinv(L, Ec, Neig, Eg, phi);
%[err_post1,~,~, ~,~,~,~,~] = PMresidual(L, Ec, Neig, Eg, phi);
%[err_post2,~,~, ~,~,~,~,~] = PMresidual(L, Ec, Neig, Em, phi);
f = err_post1(l) - err_post2(l) - eps;
 fprintf('Em = %f\n',Em);
delta = abs(Er - El);
while delta > bieps
    if   f > 0
        El = Em;
    else
        Er = Em;
    end
    delta = abs(Er - El);
    Em = (Er + El)/2;
[err_post1,~] = residualinv(L, Ec, Neig, Eg, phi);
[err_post2,~] = residualinv(L, Ec, Neig, Eg, phi);
%[err_post1,~,~, ~,~,~,~,~] = PMresidual(L, Ec, Neig, Eg, phi);
%[err_post2,~,~, ~,~,~,~,~] = PMresidual(L, Ec, Neig, Em, phi);
f = err_post1(l) - err_post2(l) - eps;
    fprintf('Em = %f\n',Em);
end
fprintf('maxE = %f\n',Em);
Eopt2 = Em;
return