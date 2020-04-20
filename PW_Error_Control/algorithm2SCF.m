function [lambda, phi, errorSCF, errorEc] = algorithm2SCF(L, Ec, tol)
% solve the nonlinear GP eigenvalue problem :
% (£­Delta+V_{ext}+\beta\rho)\phi=\lambda \phi with periodic boundary condition on [-L:L]
% using SCF iteration balance iteration and discretization error.
% Input:
% L: width of the domain
% Ec: inital energy cutoff
% tol : tolrence
% Output:
% lambda: eigenvalues
% phi: eigenfunctions
% error_SCF: iteration error
% error_Ec: discretization error
beta = 1; %nonlinear term coefficient
Neig = 1; l = 1; %ground state
maxscfk = 100; %max of scf step
maxk = 10; %max of algorithm 1
%% Initiation
Ec0 = Ec;
scfk = 0;
Eg = 1000;
% use linear solution as initiation
[~,phi,~,~] = solve_eigen(L, Ec0, Neig, Eg);
[u] = backrealnew(L, Ec0, Eg, phi);
rhon = abs(u).^2;
rho = rhon./norm(rhon);
% solve to obtain discretization error
[lambda,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
% update rho
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho+0.7*abs(u).^2;
scfk = scfk + 1;
errorSCF(scfk) = 1;
errorEc(scfk) = err_post;
% iteration 
[lambdanew,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho+0.7*abs(u).^2;
scfk = scfk + 1;
% to obtain iteration error 
delta = abs(lambdanew - lambda);
lambda = lambdanew;
errorSCF(scfk) = delta;
errorEc(scfk) = err_post;
% to obtain appropriate Ec match iteration error
bieps = 10; % set the parameter of strategy B
[~, ~, ~, Ec0] = algorithm1_2DGPnew(L, Ec0, Neig, l, Eg, delta, bieps, maxk, beta, rho);
%[u] = backrealnew(L, Ec0, Eg, phi);
%rho = 0.7*rho + 0.3*abs(u).^2;
[lambdanew,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
scfk = scfk + 1;
delta = abs(lambdanew - lambda);
lambda = lambdanew;
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
errorSCF(scfk) = delta;
errorEc(scfk) = err_post;
%[u] = backrealnew(L, Ec0, Eg, phi);
while  err_post > delta
% if discretization error dominate, refine
fprintf('At the %d-th step, E_cut = %f, err_post = %e > delta = %e , \n', scfk, Ec0, err_post, delta);
[u] = backrealnew(L, Ec0, Eg, phi);
[~, ~, ~, Ec0] = algorithm1_2DGPnew(L, Ec0, Neig, l, Eg, delta, bieps,  maxk, beta, rho);
rho = 0.3*rho + 0.7*abs(u).^2;
[lambdanew,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
delta = abs(lambdanew - lambda);
scfk = scfk + 1;
if scfk > maxscfk
   fprintf('Exit since the scf iteration steps has been reached!\n');
 break
end
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho + 0.7*abs(u).^2;
errorSCF(scfk) = delta;
errorEc(scfk) = err_post;
fprintf('At the %d-th step ... E_cut = %f, a posteriori error = %e, scf error = %e \n', scfk, Ec0, err_post, delta);
if delta < tol
 fprintf('Has been convergened! \n');
 return
end
end
while delta > tol || delta == tol
% if iteration error > tol, iterate    
fprintf('At the %d-th step, delta > tolin || delta == tolin, \n', scfk);
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho + 0.7*abs(u).^2;
lambda = lambdanew;
[lambdanew,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
[u] = backrealnew(L, Ec0, Eg, phi);
delta = abs(lambdanew - lambda);
lambda = lambdanew;
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
scfk = scfk + 1;
errorSCF(scfk) = delta;
errorEc(scfk) = err_post;
if delta < tol
 fprintf('Has been convergened! \n');
 return
end
while  err_post > delta
% if discretization error dominate, refine
fprintf('At the %d-th step, E_cut = %f, err_post = %e > delta = %e , \n', scfk, Ec0, err_post, delta);
[~, ~, ~, Ec0] = algorithm1_2DGPnew(L, Ec0, Neig, l, Eg, delta, bieps, maxk, beta, rho);
rho = 0.3*rho + 0.7*abs(u).^2;
[lambdanew,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
delta = abs(lambdanew - lambda);
scfk = scfk + 1;
lambda = lambdanew;
if scfk > maxscfk
   fprintf('Exit since the scf iteration steps has been reached!\n');
 break
end
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho + 0.7*abs(u).^2;
errorSCF(scfk) = delta;
errorEc(scfk) = err_post;
fprintf('At the %d-th step ... E_cut = %f, a posteriori error = %e, scf error = %e \n', scfk, Ec0, err_post, delta);
if delta < tol
 fprintf('Has been convergened! \n');
 return
end
end
end
fprintf('Has been convergened! \n');