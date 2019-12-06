function [lambda, phi, err_post, Ec0] = algorithm1_2DGP(L, Ec0, Neig, l, Eg, eps, bieps, maxk, beta, rho)
Ecv = zeros(maxk, 1);
indiv = zeros(maxk, 1);
% solve
[lambda,phi] = solve_eigenGP(L, Ec0, Neig, Eg, beta, rho);
Ecv(1) = Ec0;
% estimate
[err_post] = residualcrrGP(L, Ec0, Neig, Eg, phi, beta, rho);
err_post = err_post(l);
indiv(1) = err_post;
% compare
iter_k = 1;
while err_post > eps
      % Eopt
      if iter_k == 1;
         Eopt1 = NaN;
      else 
         Ecve = Ecv(Ecv~=0);
         indive = indiv(indiv~=0);
         Eopt1 = strategy1(Ecve,indive,eps);
      end
         Eopt2 = strategy2GP(L,Ec0,Eg,Neig,l,phi,eps,bieps,beta,rho);
         Eopt = min(Eopt1, Eopt2);
         Eopt = min(Eopt, Ec0 + 50);
      % Estimate
         [lambda,phi] = solve_eigenGP(L, Eopt, Neig, Eg, beta, rho);
         [err_postnew] = residualcrrGP(L, Eopt, Neig, Eg, phi, beta, rho);
         err_postnew = err_postnew(l);
         iter_k = iter_k + 1;
         Ec0 = Eopt;
         err_post = err_postnew;
         Ecv(iter_k) = Ec0;
         indiv(iter_k) = err_post;
         if iter_k > maxk 
             break 
         end
         fprintf('Algorithm1 at the %d-th step ... E_cut = %f, a posteriori error = %e \n', iter_k, Ec0,err_post);
end 




