function [lambda, phi, err_post, Ecv, indiv,lambdav] = algorithm1(L, Ec0, Neig, l, Eg, eps, bieps, maxk)
% using algorithm 1 for linear eigenvalue problem
Ecv = zeros(maxk, 1);
indiv = zeros(maxk, 1);
lambdav = zeros(maxk, 1);
% solve
[lambda,phi] = solve_eigen(L, Ec0, Neig, Eg);
lambdav(1) = lambda;
Ecv(1) = Ec0;
% estimate
%[err_post,~] = residualinv(L, Ec0, Neig, Eg, phi);
[err_post,~,~, ~,~,~,~,~] = PMresidual(L, Ec0, Neig, Eg, phi);
err_post = err_post(l);
indiv(1) = err_post;
% compare
iter_k = 1;
while err_post > eps
      % Eopt
      if iter_k == 1
         Eopt1 = NaN;
      else 
         Ecve = Ecv(Ecv~=0);
         indive = indiv(indiv~=0);
         Eopt1 = strategyA(Ecve,indive,eps);
      end
       Eopt2 = strategyBlinear(L,Ec0,Neig,l,phi,eps,bieps);
       Eopt = min(Eopt1, Eopt2);
      % Estimate
         [lambda,phi] = solve_eigen(L, Eopt, Neig, Eg);
         %[err_postnew,~] = residualinv(L, Eopt, Neig, Eg, phi);
         [err_postnew,~,~, ~,~,~,~,~] = PMresidual(L, Eopt, Neig, Eg, phi);
         err_postnew = err_postnew(l);
         iter_k = iter_k + 1;
         Ec0 = Eopt;
         err_post = err_postnew;
         Ecv(iter_k) = Ec0;
         indiv(iter_k) = err_post;
         lambdav(iter_k) = lambda;
         if iter_k > maxk 
             break 
         end
         fprintf('Algorithm1 at the %d-th step ... E_cut = %f, a posteriori error = %e \n', iter_k, Ec0, err_post);
end 
