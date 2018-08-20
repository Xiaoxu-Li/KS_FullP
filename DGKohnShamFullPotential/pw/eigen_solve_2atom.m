function [lambda,phi]=eigen_solve_2atom(L, Nc)
% eigen_solve_2atom solve the eigenvalue problem with periodic potential
% such as hydrogen H^2
% veff(r)=-\frac{1}{|\Omega|}*\sum_{k\neq 0}\frac{4\pi*e^{ik(r-R1)}}{|k|^2}
%             -\frac{1}{|\Omega|}*\sum_{k\neq 0}\frac{4\pi*e^{ik(r-R2)}}{|k|^2}

% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain [-L/2,L/2]^3
% Nc is the number of planewaves used in each direction.
% Aug 16th, 2018

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
R1=[-2,0,0];  % positions of atoms
R2=[2,0,0];

%======================== calculate the matrix element ===================%
for ii=1:n
    for j=1:n
        dk = p(j,:)-p(ii,:);
        shift1=exp(1i*dk*R1'*2*pi/L);
        shift2=exp(1i*dk*R2'*2*pi/L);
        if p(ii,:) == p(j,:)
            H(ii,j) = 0.5*norm(p(ii,:),2)^2 *(2*pi/L)^2;
            M(ii,j) = 1.0;
        else
            H(ii,j)=-(shift1+shift2)*4*pi/(L^3*norm(dk,2)^2*(2*pi/L)^2);   
        end          
    end
end
%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
H=(H+H')/2;
%[phi,lambda] = eigs(H, M, 1, 'sa');
[phi,lambda] = eigs(H, M, 1, -100);
c=phi;
nrm=sqrt(c'*M*c);
save pw_per L  c nrm p n


