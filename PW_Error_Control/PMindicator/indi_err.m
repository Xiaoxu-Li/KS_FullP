%clear all
% use the solution (lambda,phi) with Ec=500 as reference 
% Nc = sqrt(Ec), plot graphs of residual and |\lambdaN-\lambda| 
L = pi;
Neig= 1;
Ecr = 500;
Eg = 1000;
%%
[lambdar,~] = solve_eigen(L, Ecr, Neig, Eg);
%%
for Ec = 10:5:120
    Eg = 4*Ec;
    [lambda,phi] = solve_eigen(L, Ec, Neig, Eg);
    [err_post,err_post2,res,invlapres_fftk] = PMresidual(L, Ec, Neig, Eg, phi);
    [err_postinv] = residualinvcrr(L, Ec, Neig, Eg, phi);
    for ll=1:Neig
    Delta =abs(lambda(ll) - lambdar(ll));
    semilogy(sqrt(Ec),Delta,'ro--','LineWidth', 2.5, 'MarkerSize', 15);
    hold on 
    semilogy(sqrt(Ec),err_post(ll),'bs--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_post2(ll),'g*--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_postinv(ll),'kh--','LineWidth', 2.5, 'MarkerSize', 15);
    end
end
s = xlabel('$\sqrt{E_c}$ (log)') ;
 ylabel('Eigenvalue errors (log)') 
a = legend({'$\lambda_{E_{c}} - \lambda$','$\eta_{E_{c}}^{[1]}$','$\eta_{E_{c}}^{[2]}$','$\eta_{E_{c}}$'},'Location','best');
title('1st eigenvalue error for V_1')
set(gca,'Fontsize',25)
set(a,'Interpreter','latex')
set(s,'Interpreter','latex')