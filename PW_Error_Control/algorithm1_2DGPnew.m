function [lambda, phi, err_post, Ec0] = algorithm1_2DGPnew(L, Ec0, Neig, l, Eg, eps, bieps, maxk, beta, rho)
% using algorithm 1 to solve linearization GP eigenvalue problem
% (£­Delta+V_{ext}+\beta\rho)\phi=\lambda \phi with periodic boundary condition on [-L:L]
Ecv = zeros(maxk, 1);
indiv = zeros(maxk, 1);
% solve
[lambda,phi] = solve_eigenGPnew(L, Ec0, Neig, Eg, beta, rho);
[u] = backrealnew(L, Ec0, Eg, phi);
rho = abs(u).^2;
Ecv(1) = Ec0;
% estimate
[u] = backrealnew(L, Ec0, Eg, phi);
rho = 0.3*rho+0.7*abs(u).^2;
[err_post, ~, ~, ~] = EcutestiGP(L, Ec0, Neig, Eg, phi, beta, rho);
err_post = err_post(l);
indiv(1) = err_post;
% compare
iter_k = 1;
while err_post > eps
      if iter_k == 1
         Eopt1 = NaN;
      else 
         Ecve = Ecv(Ecv~=0);
         indive = indiv(indiv~=0);
         Eopt1 = strategyA(Ecve,indive,eps);
      end  
         Eopt2 = strategyB(L,Ec0,Eg,Neig,l,phi,eps,bieps,beta,rho);
         Eopt = min(Eopt1, Eopt2);
      % Estimate
         [lambda,phi] = solve_eigenGPnew(L, Eopt, Neig, Eg, beta, rho);
         [err_postnew, ~, ~, ~] = EcutestiGP(L, Eopt, Neig, Eg, phi, beta, rho);
         err_postnew = err_postnew(l);
         iter_k = iter_k + 1;
         Ec0 = Eopt;
         err_post = err_postnew;
         Ecv(iter_k) = Ec0;
         indiv(iter_k) = err_post;
         if iter_k > maxk 
         return
         end
         fprintf('Algorithm1 at the %d-th step ... E_cut = %f, a posteriori error = %e \n', iter_k, Ec0,err_post);
end