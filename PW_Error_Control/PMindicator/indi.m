% clear all
% use the solution (lambda,phi) with Ec=500 as reference 
% Nc = sqrt(Ec), plot graphs of residual and |\lambdaN-\lambda| 
L = 5;
Neig= 1;
Ecr = 100;
Egf = 200; % (to calculate the potential element)
%%
[lambdar,~,dofg] = solve_eigen(L, Ecr, Neig, Egf);
%%
for Ec=1:3:30
    Eg = 4*Ec;
    [lambda,phi,~,vec_k] = solve_eigen(L, Ec, Neig, Egf);
    [err_post,err_post2,err_post12, invlapres_fftk,Hv,Hpiv,res] = PMresidual(L, Ec, Neig, Eg, phi);
    [err_postinv,Hm] = residualinv(L, Ec, Neig, Eg, phi);
    %for ll=1:Neig
    ll = 1;
    Delta =abs(lambda(ll) - lambdar(ll));
    semilogy(sqrt(Ec),Delta,'ro--','LineWidth', 2.5, 'MarkerSize', 15);
    hold on
    semilogy(sqrt(Ec),err_post(ll),'bs--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_post2(ll),'g*--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_post12(ll),'c^--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_postinv(ll),'kh--','LineWidth', 2.5, 'MarkerSize', 15);
    %semilogy(sqrt(Ec),Delta - err_post(ll) ,'r^--','LineWidth', 2.5, 'MarkerSize', 15);
    %semilogy(sqrt(Ec),Delta - err_post2(ll)  ,'b^--','LineWidth', 2.5, 'MarkerSize', 15);
    %semilogy(sqrt(Ec),Delta - err_post12(ll) ,'g^--','LineWidth', 2.5, 'MarkerSize', 15);
    %end
end
s = xlabel('$\sqrt{E_c}$ (log)');
 ylabel('Eigenvalue errors (log)') 
%a = legend({'$\lambda_{E_{c}} - \lambda$','${\|u-u_N\|_a}^2$','$||u-u_N||^2_a - \lambda{\|u-u_N\|_2}^2$','$\eta_{E_{c}}^{[1]}$','$\eta_{E_{c}}^{[2]}$','$\eta_{E_{c}}$'},'Location','best');
a = legend({'$\lambda_{E_{c}} - \lambda$','$\eta_{E_{c}}^{[1]}$','$\eta_{E_{c}}^{[2]}$','$\eta_{E_{c}}^{old,1}$','$\eta_{E_{c}}$'},'Location','best');
title('1st eigenvalue error for V')
set(gca,'Fontsize',25)
set(a,'Interpreter','latex')
set(s,'Interpreter','latex')
axis([2 6 1e-10 1e-2])
