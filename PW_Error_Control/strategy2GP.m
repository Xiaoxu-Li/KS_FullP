% strategy 2
% predict E_opt by bi-section method
function [Eopt2] = strategy2GP(L,Ec,Eg,Neig,l,phi,eps,bieps,beta,rho)
% Ecv Ecut vector;
El = Ec;
Er = Eg;
Em = (Er + El)/2;
err_post = residualcrrGP(L, Ec, Neig, Eg, phi, beta, rho);
err_post_ = residualcrrGP(L, Ec, Neig, Em, phi, beta, rho);
f = err_post(l) - err_post_(l) - eps;
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
err_post = residualcrrGP(L, Ec, Neig, Eg, phi, beta, rho);
err_post_ = residualcrrGP(L, Ec, Neig, Em, phi, beta, rho);
f = err_post(l) - err_post_(l) - eps;
    fprintf('Em = %f\n',Em);
end
fprintf('maxE = %f\n',Em);
Eopt2 = Em;
return