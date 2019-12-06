%clear all
% use the solution (lambda,phi) with Ec=500 as reference 
% Nc = sqrt(Ec), plot graphs of residual and |\lambdaN-\lambda| 
L = 2*pi;
Neig= 1;
Ecr = 1000;
Eg = 4000;
%%
[lambdar,~] = solve_eigen(L, Ecr, Neig, Eg);
%%
for Ec=20:20:400
    [lambda,phi] = solve_eigen(L, Ec, Neig, Eg);
    [err_post,err_post2] = Dresidualcrr(L, Ec, Neig, Eg, phi);
    [err_postinv] = residualinv(L, Ec, Neig, Eg, phi);
    %for ll=1:Neig
    ll = 1;
    Delta =abs(lambda(ll) - lambdar(ll));
    semilogy(sqrt(Ec),Delta,'ro--','LineWidth', 2.5, 'MarkerSize', 15);
    hold on 
    semilogy(sqrt(Ec),err_post(ll),'bs--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_post2(ll),'g*--','LineWidth', 2.5, 'MarkerSize', 15);
    semilogy(sqrt(Ec),err_postinv(ll),'kh--','LineWidth', 2.5, 'MarkerSize', 15);
    %end
end
s = xlabel('$\sqrt{E_c}$ (log)') 
 ylabel('Eigenvalue errors (log)') 
a = legend({'$\lambda_{1,Ec} - \lambda_{1}$','$\eta^{1}_{Ec,1,H^{-1}}$','$\eta^{1}_{Ec,1,GD}$','$\eta^{1}_{Ec,1,inv}$'},'Location','best')
title('1st eigenvalue error for V_3')
set(gca,'Fontsize',25)
set(a,'Interpreter','latex')
set(s,'Interpreter','latex')