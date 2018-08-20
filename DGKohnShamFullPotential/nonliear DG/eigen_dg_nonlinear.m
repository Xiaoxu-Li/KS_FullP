function [energy,lambda_new,phi]=eigen_dg_nonlinear(L, R, N_p, n_r, Lm, sigma)
% eigen_dg_nonlinear solve the nonlinear eigenvalue problem    
% (-1/2\Delta+V_{ext}+V_H{\rho})u=\lambda u [\lambda,\phi] returns the
% eigenpair and energy
% L is the width of the domain [-L/2,L/2]^3, R is the radious of sphere
% N_p is the number of planewaves used in each direction
% n_r is the number of radious basis, Lm is the number of angular momentum
% sigma is the penalization parameter
% SCF iteration throgh simple mixing for \rho_in=alpha*\rho_in+(1-alpha)*\rho_out
% load jump data 
% Aug 20th, 2018

%============================= initiate the data =========================%
% initial_data=0
% degrees of freedom
% t0=cputime;
load jump_data_R_3
n0=0;
N=floor(N_p);
kk=zeros(8*N^3,3);
for ii=-N:N
    m=sqrt(N_p^2-ii^2);
    m=floor(m);
    for j=-m:m
        l=sqrt(N_p^2-ii^2-j^2);
        l=floor(l);
        for k=-l:l
            n0=n0+1;
            kk(n0,:)=[ii,j,k];
        end
    end
end
kk(n0+1:8*N^3,:)=[];    
n = n0 + n_r*(Lm+1)^2;
fprintf('pw DOF=%d, total DOF=%d\n ',n0 ,n );
H = zeros(n,n);
M = zeros(n,n);
A = zeros(n,n);    % discretization of laplace operator
n_integral = 100;
delta=0.000001;
% lmax = 5;
stablizor=1;

% cut-off for periodic potential
N0=0;
N_cut=2*N_p;    % 2~4 times of N_p
N1=floor(N_cut);
% N_cut=2*N;    %  if N_p is half integer, N1=2*N+1 maybe beyond the index below
% N1=floor(N_cut);
kkk=zeros(8*N1^3,3);
for ii=-N1:N1
      m=sqrt(N_cut^2-ii^2);
      m=floor(m);
      for j=-m:m
            l=sqrt(N_cut^2-ii^2-j^2);
            l=floor(l);
            for k=-l:l
                  N0=N0+1;
                  kkk(N0,:)=[ii,j,k];
            end
      end
end
fprintf('cut-off for periodic potential=%d\n ',N0);

volum_omega = L^3;
volum_in = 4/3*pi*R^3;
volum_out = volum_omega-volum_in;

% periodic potential 
mm=4*N+1;
vext = zeros(mm,mm,mm);
for ii = 1:mm
    for j = 1:mm
        for k = 1:mm
            sum_per=0;
            x = ii-2*N-1;
            y = j-2*N-1;
            z = k-2*N-1;
            K_index = [x,y,z];
            K_norm=norm(K_index,2)*2*pi/L;
                if K_norm==0   
                    for s=1:N0
                        K_per=norm(kkk(s,:),2)*2*pi/L;
                        if K_per==0
                            continue;
                        end
                        sum_per=sum_per+sin(K_per*R)/K_per^5-cos(K_per*R)*R/K_per^4;        
                    end
                    v_per=sum_per*(4*pi)^2/(volum_omega)^2;
                else
                    for s=1:N0
                        KK=norm(kkk(s,:)+K_index,2)*2*pi/L;
                        K_per=norm(kkk(s,:),2)*2*pi/L;
                        if KK==0
                            continue;
                        elseif K_per==0
                            continue;
                        end
                        sum_per=sum_per+1/K_per^2*(sin(KK*R)/KK^3-R*cos(KK*R)/KK^2);
                    end
                    v_per=-4*pi/(K_norm^2*volum_omega)+sum_per*(4*pi)^2/(volum_omega)^2 ...
                            +4*pi*volum_in/(K_norm^2*volum_omega^2);    
                end                    
                vext(ii, j, k)=vext(ii, j, k)+v_per;
        end
    end
end

%==================================================== calculate the matrix element ===============================================%

% [1] ======================== matrix element for planewaves =======================%
matrix_pw=1
for p=1:n0
    for q=1:n0
        dk = kk(q,:)-kk(p,:);      
        K = norm(dk,2) * 2*pi/L;   
        % overlap
        if p == q
            M(p,q) = volum_out/volum_omega; 
        else
            M(p,q) = -4*pi*R^2  * spherical_bessel(1,K*R)/(K*volum_omega);
        end
        % kinetic
        H(p,q) = 0.5 * (kk(p,:)*kk(q,:)') * M(p,q) * (2*pi/L)^2;         
        A(p,q)=2*H(p,q);
        % external potential H_{pq}=vext_fft(p-q)

        % periodic potential   
        index=dk+2*N+1;
        H(p,q)=H(p,q)+2*vext(index(1),index(2),index(3));     % external potential of He
        
         % discontinuous and penalization
         kp = norm(kk(p,:),2) * 2*pi/L;
         kq = norm(kk(q,:),2) * 2*pi/L;
         % case kp=kq=0
         if kp == 0 && kq == 0
             H(p,q) = H(p,q) + 4*pi*R^2 * sigma / volum_omega;
             A(p,q) = A(p,q) + 4*pi*R^2 * sigma / volum_omega;    % penalization not change
         % case kp=0  and case kq=0
         elseif kp > 0 && kq == 0
             djp= (spherical_bessel(0,kp*(R+delta))-spherical_bessel(0,kp*(R-delta)))/(2*delta);
             jp= spherical_bessel(0,kp*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djp+ 4*sigma * jp);
             A(p,q) = A(p,q) + 4*pi * 0.25*R^2/volum_omega * (2*djp+ 4*sigma * jp);  % almost double penalization
         elseif kp == 0 && kq > 0
             djq= (spherical_bessel(0,kq*(R+delta))-spherical_bessel(0,kq*(R-delta)))/(2*delta);
             jq=spherical_bessel(0,kq*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djq + 4*sigma * jq);
             A(p,q) = A(p,q) + 4*pi * 0.25*R^2/volum_omega * (2*djq + 4*sigma * jq);  % almost double penalization
         % case kp>0 and kq>0    
         elseif kp > 0 && kq > 0
             j_index=[kk(p,1),kk(p,2),kk(p,3),kk(q,1),kk(q,2),kk(q,3)]+5+1;
             sum=jump( j_index(1),  j_index(2),  j_index(3),  j_index(4),  j_index(5),  j_index(6));                
             H(p,q) = H(p,q) + 0.25*(4*pi)^2*sum * R^2/volum_omega;
             A(p,q) = A(p,q) + 0.5*(4*pi)^2*sum * R^2/volum_omega;
         end
    end
end

% [2] ====================== matrix elements for radial type basis ==========================%
% calculate the integration of radial basis functions
matrix_element_r=2
int_r = zeros(n_r,n_r);
int_dr = zeros(n_r,n_r);
int = zeros(n_r,n_r);
int_v=zeros(n_r,n_r);
for k1 = 1:n_r
    for k2 = 1:n_r
        int_r(k1,k2) = int_overlap_radius(k1,k2,R,n_integral);
        int(k1,k2) = int_overlap(k1,k2,R,n_integral);
        int_dr(k1,k2) = int_overlap_radius_d(k1,k2,R,n_integral);
        % periodic potential
        sum_per_in=0;
        for j=1:N0
            K_per=norm(kkk(j,:),2)*2*pi/L;
            if K_per==0
                continue;
            end
            sum_per_in=sum_per_in+int_vext_per_sin(K_per,k1,k2,R,n_integral);      
        end            
        int_v(k1,k2) = -4*pi*sum_per_in/volum_omega;
    end
end
% matrix elements
for k1 = 1:n_r
    for l1 = 0:Lm
        for m1 = -l1:l1
            for k2 = 1:n_r
                for l2 = 0:Lm
                    for m2 = -l2:l2
                    %-----------------------------------------------------------------------------------------------------------------%    
                        p = n0 + (k1-1)*(Lm+1)^2 + l1^2 + m1+l1+1;
                        q = n0 + (k2-1)*(Lm+1)^2 + l2^2 + m2+l2+1;
                        % overlap and kinetic
                        if  l1==l2 && m1==m2
                            M(p,q) = int_r(k1,k2);    
                            kinetic_lm = l1*(l1+1) * int(k1,k2);
                            kinetic_r = int_dr(k1,k2);
                            H(p,q) = kinetic_lm + kinetic_r; 
                            A(p,q) = H(p,q);
                            H(p,q) = 0.5 * H(p,q) + 2*int_v(k1,k2);    % external potential of He
                        else
                            M(p,q) = 0.0;
                            H(p,q) = 0.0;
                            A(p,q) = 0.0;
                        end                                                
                        % discontinuous and penalization     
                        if l1==l2 && m1==m2
                            rp = -r_basis(k1,R,R);
                            drp = (r_basis(k1,R+delta,R)-r_basis(k1,R-delta,R))/(2*delta);
                            rq = -r_basis(k2,R,R);
                            drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta); 
                            H(p,q) = H(p,q) + 0.25*(rp*drq +drp*rq + 4*sigma*rp*rq) * R^2; 
                            A(p,q) = A(p,q) + 0.25*(2*rp*drq +2*drp*rq + 4*sigma*rp*rq) * R^2;
                        end
                    end
                end
            end
        end
    end
end

% [3] ====================== matrix elements for planewaves and radial type basis =======================%
matrix_element_pw_r =3
for p = 1:n0
    for k2 = 1:n_r
        for l2 = 0:Lm
            for m2 = -l2:l2
                %-------------------------------------------------------------------------------------------------------------------%
                kp = norm(kk(p,:),2) * 2*pi/L;
                q = n0 + (k2-1)*(Lm+1)^2 + l2^2 + m2+l2+1;
                % overlap is 0
                M(p,q) = 0.0;
                M(q,p) = M(p,q);                        
                % only discontinuous and penalization is considered for H(p,q)
                if kp == 0
                    if l2 == 0
                        rq = -r_basis(k2,R,R);
                        drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);
                        jp = 1;
                        H(p,q) = 0.25*R^2 * sqrt(4*pi)/sqrt(volum_omega) * (jp*drq + 4*sigma*rq*jp);
                    else
                        H(p,q) = 0.0;
                    end
                else
                    ylm_part = spherical_harmonic_xyz(l2,m2,-kk(p,1),-kk(p,2),-kk(p,3));  
                    jp =  spherical_bessel(l2,kp*R);
                    djp = (spherical_bessel(l2,kp*(R+delta))-spherical_bessel(l2,kp*(R-delta)))/(2*delta);
                    rq = -r_basis(k2,R,R);
                    drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);
                    H(p,q) = 0.25*4*pi * 1i^l2 * ylm_part * (jp*drq +djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);  
                    A(p,q)=0.25*4*pi * 1i^l2 * ylm_part * (2*jp*drq +2*djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);    
                end
                H(q,p) = H(p,q)';              
                A(q,p) =A(p,q)';
            end
        end
    end
end
%==================================== end of matrix element computation ==============================%

%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
eigen_solver=4
A=A+stablizor*M;
H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1,  -100);
c=phi ;
nrm=sqrt(c'*M*c);

% SCF part 
iteration=5
lambda_old=0;
lambda_new=lambda;
f_in=zeros(n,1);    %  \rho_in
iter=0;
tol=10^(-3);

% truncation of polynomials
 int_3=zeros(n_r,n_r,n_r);
 int_3r=zeros(n_r,n_r,n_r); 
 for k1 = 1:n_r
    for k2 = 1:n_r
        for k3=1:n_r
            int_3(k1,k2,k3)=int_overlap_3(k1, k2, k3, R, n_integral);
            int_3r(k1,k2,k3)=int_overlap_radius_3(k1, k2, k3, R, n_integral);
        end
    end
 end
 
 % begin iteration
while abs(lambda_new-lambda_old)>tol && iter<=100   
coef_original=zeros(4*N+3,4*N+3,4*N+3);
 for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)]+2*N+2;
        coef_original(k(1),k(2),k(3))= c(p)/nrm;
 end

 % compute coefficient of planewaves
 % DOF=4*N+1  |k|<2*K
 coef_new=zeros(4*N+3,4*N+3,4*N+3);
 for i1=1:2*N+2
     for j1=1:2*N+2
         for k1=1:2*N+2
            for i2=1:2*N+2
                for j2=1:2*N+2
                    for k2=1:2*N+2
                        coef_new(i1-i2+2*N+2, j1-j2+2*N+2, k1-k2+2*N+2)=...
                            coef_new(i1-i2+2*N+2, j1-j2+2*N+2, k1-k2+2*N+2)+...
                            coef_original(i1,j1,k1)*coef_original(i2,j2,k2)';                
                    end
                end
            end
         end
     end
 end
 
%  using DG method to discretize the problem
% -\Delta V_H{\rho}=4*pi*{\rho} 
% assemble stiffness matrix A and load vector f
% solve Ax=f

% element for vector f
% [1] planewaves part
 f_out=zeros(n,1);   %  \rho_out
 
 for ii=1:n0
     sum_per=0;
     for j=1:N0
         k=kkk(j,:)+2*N+2;
         KK=norm(kkk(s,:)-kk(ii,:),2)*2*pi/L;
         if KK==0
                continue;
         end
         sum_per=sum_per+coef_new(k(1),k(2),k(3))*(sin(KK*R)/KK^3-R*cos(KK*R)/KK^2);         
     end
     k_index=kk(ii,:)+2*N+2;
     f_out(ii)=coef_new(k_index(1),k_index(2),k_index(3))*4*pi*volum_out/(volum_omega)^(3/2)...
            - sum_per*(4*pi)^2/(volum_omega)^(3/2);              
 end

 % [2] radial basis functions part
 for k3=1:n_r
     for l3=0:Lm
         for m3=-l3:l3
             ii=n0 + (k3-1)*(Lm+1)^2 + l3^2 + m3+l3+1;
            for k1 = 1:n_r
                for l1 = 0:Lm
                    for m1 = -l1:l1
                         for k2 = 1:n_r
                            for l2 = 0:Lm
                                for m2 = -l2:l2
                                    p = n0 + (k1-1)*(Lm+1)^2 + l1^2 + m1+l1+1;
                                    q = n0 + (k2-1)*(Lm+1)^2 + l2^2 + m2+l2+1;
                                    f_out(ii)=f_out(ii)+c(p)*c(q)'/nrm^2*int_3(k1,k2,k3)*(-1)^(m2)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1))...
                                        /sqrt(4*pi)*Wigner3j(l1,l2,l3,0,0,0)*Wigner3j(l1,l2,l3,m1,-m2,-m3);             
                                end
                            end
                         end                         
                    end
                end
            end  
            f_out(ii)=(-1)^(m3)*f_out(ii);
        end
    end
 end
 
 f_out=2*f_out;   % \rho=2*|u|^2  including spin
 
 if iter==0
     f_in=f_out;   % \rho_in=\rho_out in the initial condition
 end
 f_in=9*f_in/10+f_out/10;   % simple mixing for \rho_in=alpha*\rho_in+(1-alpha)*\rho_out

% A is DG discretization of Laplace operator
 X=A\f_in;

% Hartree potential 
mm=4*N+1;
v_H= zeros(mm,mm,mm);
for ii = 1:mm
    for j = 1:mm
        for k = 1:mm
            sum_per=0;
            x = ii-2*N-1;
            y = j-2*N-1;
            z = k-2*N-1;
            K_index = [x,y,z];
            if norm(K_index,2)>N_p   
                for s=1:n0
                    KK=norm(kk(s,:)+K_index,2)*2*pi/L;
                    sum_per=sum_per+X(s)*(sin(KK*R)/KK^3-R*cos(KK*R)/KK^2);
                end
                v_Hartree=-4*pi*sum_per/(volum_omega)^(3/2);
            else
                target=find(kk(:,1)==-x & kk(:,2)==-y & kk(:,3)==-z);                                    
                for s=1:n0
                    KK=norm(kk(s,:)+K_index,2)*2*pi/L;
                    if KK==0
                        continue;
                    end
                    sum_per=sum_per+X(s)*(sin(KK*R)/KK^3-R*cos(KK*R)/KK^2);
                end
                v_Hartree=-4*pi*sum_per/(volum_omega)^(3/2)+X(target)*volum_out/(volum_omega)^(3/2);    
            end                    
            v_H(ii, j, k)=v_H(ii, j, k)+v_Hartree;
        end
    end
end
 
% matrix element for Hartree potential
 H_iter=zeros(n,n);
 
 % [1] planewaves part
 for p=1:n0
    for q=1:n0
        dk = kk(q,:)-kk(p,:);      
        index=dk+2*N+1;
        H_iter(p,q)=H_iter(p,q)+v_H(index(1),index(2),index(3));
    end
 end
 
 % [2] radial basis functions part
for k1 = 1:n_r
    for l1 = 0:Lm
        for m1 = -l1:l1
            for k2 = 1:n_r
                for l2 = 0:Lm
                    for m2 = -l2:l2
                        p = n0 + (k1-1)*(Lm+1)^2 + l1^2 + m1+l1+1;
                        q = n0 + (k2-1)*(Lm+1)^2 + l2^2 + m2+l2+1;
                        for k3=1:n_r
                            for l3=0:Lm
                                for m3=-l3:l3
                                    ii=n0 + (k3-1)*(Lm+1)^2 + l3^2 + m3+l3+1;
                                    H_iter(p,q)=H_iter(p,q)+X(ii)*int_3r(k1,k2,k3)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1))/sqrt(4*pi)...
                                    *Wigner3j(l1,l2,l3,0,0,0)*Wigner3j(l1,l2,l3,-m1,m2,m3);            
                                end
                            end
                        end
                        H_iter(p,q)=H_iter(p,q)*(-1)^(m1);
                    end
                end
            end           
        end
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
% save MAT L R n0 N_p n_r Lm c nrm kk volum_omega delta

% t=(cputime-t0)/3600;
% fprintf('time=%d h\n',t);

return;
