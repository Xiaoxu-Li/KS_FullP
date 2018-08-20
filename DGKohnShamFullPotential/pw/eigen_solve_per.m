function [lambda,phi]=eigen_solve_per(L, Nc)
% eigen_solve_per solve the eigenvalue problem with periodic version of
% -1/r       
% veff(r)=-\frac{1}{|\Omega|}*\sum_{k\neq 0}\frac{4\pi*e^{ikr}}{|k|^2}
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain  [-L/2,L/2]^3
% Nc is the number of planewaves used in each direction.
 
%============================= initiate the data =========================%
n=0; % number of basis satisfying |k|<Nc
N=floor(Nc);
p=zeros(8*N^3,3);
for ii=-N:N
    m=sqrt(Nc^2-ii^2);
    m=floor(m);
    for j=-m:m
        l=sqrt(Nc^2-ii^2-j^2);
        l=floor(l);
        for k=-l:l
            n=n+1;
            p(n,:)=[ii,j,k];
        end
    end
end 
fprintf('DOF=%d\n', n)
H = zeros(n,n);
M = zeros(n,n);

%======================== calculate the matrix element ===================%
for ii=1:n
    for j=1:n
        if p(ii,:) == p(j,:)
            H(ii,j) = 0.5*norm(p(ii,:),2)^2 *(2*pi/L)^2;
            M(ii,j) = 1.0;
        else
            H(ii,j)=-4*pi/((L^3)*norm(p(ii,:)-p(j,:),2)^2*(2*pi/L)^2);
        end           
    end
end
%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1, 'sa');
 c=phi;
 nrm=sqrt(c'*M*c);
 save pw_per L  c nrm p n


