function [energy,lambda,phi]=eigen_solve_nonlinear(L, Nc)
% eigen_solve_nonlinear solve the nonlinear eigenvalue problem with periodic potential
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain  [-L/2,L/2]^3
% Nc is the number of planewaves used in each direction.
% SCF iteration through simple mixing 
% Aug 20th, 2018

%============================= initiate the data =========================%
% t0=cputime;
stablizor=1;
n=0; % number of basis satisfying |k|<Nc
N=floor(Nc);
kk=zeros(8*N^3,3);
for ii=-N:N
    m=sqrt(Nc^2-ii^2);
    m=floor(m);
    for j=-m:m
        l=sqrt(Nc^2-ii^2-j^2);
        l=floor(l);
        for k=-l:l
            n=n+1;
            kk(n,:)=[ii,j,k];
        end
    end
end 
fprintf('DOF=%d\n', n)
kk(n+1:8*N^3,:)=[];    
H = zeros(n,n);
M = zeros(n,n);
A= zeros(n,n);

%======================== calculate the matrix element ===================%
for ii=1:n
    for j=1:n
        if kk(ii,:) == kk(j,:)
            H(ii,j) = 0.5*norm(kk(ii,:),2)^2 *(2*pi/L)^2;
            M(ii,j) = 1.0;
            A(ii,j) =  2*H(ii,j);
        else
            H(ii,j)=-8*pi/((L^3)*norm(kk(ii,:)-kk(j,:),2)^2*(2*pi/L)^2);       % external potential of He -2/r
        end      
    end
end
%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
A=A+stablizor*M;
H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1, -100);
c=phi;
nrm=sqrt(c'*M*c);

% SCF part 
lambda_old=0;
lambda_new=lambda;
f_in=zeros(n,1);    %  \rho_in
iter=0;
tol=10^(-3);

% begin iteration
while abs(lambda_new-lambda_old)>tol && iter<=100   
coef_original=zeros(4*N+1,4*N+1,4*N+1);
 for p=1:n
        k=[kk(p,1),kk(p,2),kk(p,3)]+2*N+1;
        coef_original(k(1),k(2),k(3))= c(p)/nrm;
 end
 
% compute coefficient of planewaves
% DOF=4*N+1  |k|<2*K
 coef_new=zeros(4*N+1,4*N+1,4*N+1);
 for i1=1:2*N+1
     for j1=1:2*N+1
         for k1=1:2*N+1
            for i2=1:2*N+1
                for j2=1:2*N+1
                    for k2=1:2*N+1
                        coef_new(i1-i2+2*N+1, j1-j2+2*N+1, k1-k2+2*N+1)=...
                            coef_new(i1-i2+2*N+1, j1-j2+2*N+1, k1-k2+2*N+1)+...
                            coef_original(i1,j1,k1)*coef_original(i2,j2,k2)';                 
                    end
                end
            end
         end
     end
 end

% -1/2\Delta V_H{\rho}=4*pi*{\rho} 
% assemble stiffness matrix A and load vector f
% solve Ax=f

% element for vector f
% planewaves part
 f_out=zeros(n,1);   %  \rho_out
 
 for ii=1:n
     k_index=kk(ii,:)+2*N+1;
     f_out(ii)=coef_new(k_index(1),k_index(2),k_index(3))*4*pi/sqrt(L^3);
 end
 
 f_out=2*f_out;   % \rho=2*|u|^2  including spin
 
 if iter==0
     f_in=f_out;   % \rho_in=\rho_out in the initial condition
 end
 f_in=9*f_in/10+f_out/10;   % simple mixing for \rho_in=\alpha*\rho_in+(1-alpha)*\rho_out

% A is DG discretization of Laplace operator
 X=A\f_in;
 
 % Hartree potential 
mm=4*N+1;
v_H= zeros(mm,mm,mm);
for ii = 1:mm
    for j = 1:mm
        for k = 1:mm
            x = ii-2*N-1;
            y = j-2*N-1;
            z = k-2*N-1;
            K_index = [x,y,z];
            if norm(K_index,2)>Nc   
                continue;
            else
                target=find(kk(:,1)==-x & kk(:,2)==-y & kk(:,3)==-z);                                    
                v_H(ii, j, k)=X(target)/sqrt(L^3);
            end                                
        end
    end
end
 
 % matrix element for Hartree potential
 H_iter=zeros(n,n);
 
 % planewaves part
 for p=1:n
    for q=1:n
        dk = kk(q,:)-kk(p,:);      
        index=dk+2*N+1;
        H_iter(p,q)=H_iter(p,q)+v_H(index(1),index(2),index(3));
    end
 end
 
 % update  new Hamiltonian
iter=iter+1;
H_new=H+H_iter;
H_new=(H_new+H_new')/2;
[phi,lambda] = eigs(H_new, M, 1,  -10);
lambda_old=lambda_new;
lambda_new=lambda;

c=phi ;
nrm=sqrt(c'*M*c);
error=abs(lambda_new-lambda_old);
fprintf('iteration time=%d, lambda_old=%d , lambda_new=%d ,\n |lambda_new-lambda_old|=%d\n',iter,lambda_old,lambda_new,error)
end

energy=lambda_new-c'*H_iter*c/nrm^2;
fprintf('KS energy is %d\n', energy);
% save pw_nonlinear L R n0 N_p n_r Lm c nrm kk volum_omega delta

% t=(cputime-t0)/3600;
% fprintf('time=%d h\n',t);

return;
