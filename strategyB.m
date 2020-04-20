% strategy B
% predict E_opt by bi-section method
function [Eopt2] = strategyB(L,Ec,Neig,l,phi,eps,bieps,beta,rho)
% Ecv Ecut vector;
Eg = 4*Ec;
El = Ec;
Er = Eg;
Em = (Er + El)/2;
[err_post1, ~, ~, ~] = EcutestiGP(L, Ec, Neig, Eg, phi, beta, rho);
[err_post2, ~, ~, ~] = EcutestiGP(L, Ec, Neig, Em, phi, beta, rho);
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
[err_post1, ~, ~, ~] = EcutestiGP(L, Ec, Neig, Eg, phi, beta, rho);
[err_post2, ~, ~, ~] = EcutestiGP(L, Ec, Neig, Em, phi, beta, rho);
f = err_post1(l) - err_post2(l) - eps;
    fprintf('Em = %f\n',Em);
end
fprintf('maxE = %f\n',Em);
Eopt2 = Em;
return