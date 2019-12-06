function [lambda, phi, delta, err_post] = algorithm2(L, Ec0, Neig, l, Eg, tol, bieps, maxscfk, maxa1k, beta)
Ng = floor(sqrt(Eg));
rhon = rand(2*Ng+1, 2*Ng+1);
rho = rhon./norm(rhon);
scfk = 0;
% 1. solve
[lambda,phi] = solve_eigenGP(L, Ec0, Neig, Eg, beta, rho);
% 2. rhonew
[u] = backreal(L, Ec0, Eg, phi);
rho = 0.8*rho+0.2*abs(u).^2;
% scf
[lambdanew,phi] = solve_eigenGP(L, Ec0, Neig, Eg, beta, rho);
[u] = backreal(L, Ec0, Eg, phi);
scfk = scfk + 1;
% 3.estimate
delta = abs(lambdanew - lambda);
% 4.compare
while delta >= tol
    [err_post] = residualcrrGP(L, Ec0, Neig, Eg, phi, beta, rho)
    if err_post > delta
    [lambda, phi, err_post, Ec0] = algorithm1_2DGP(L, Ec0, Neig, l, Eg, delta, bieps, maxa1k, beta, rho);
    else
    %scf
    rho = 0.2*rho + 0.8*abs(u).^2;
    [lambdanew,phi] = solve_eigenGP(L, Ec0, Neig, Eg, beta, rho);
    [u] = backreal(L, Ec0, Eg, phi);
    delta = abs(lambdanew - lambda)
    lambda = lambdanew;
    scfk = scfk + 1
       if scfk > maxscfk
       fprintf('Exit since the scf iteration steps has been reached!\n');
        return
       end
    end
end
return