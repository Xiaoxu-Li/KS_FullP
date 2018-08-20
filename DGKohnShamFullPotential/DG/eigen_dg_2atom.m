function [lambda,phi]=eigen_dg_2atom(L, R, N_p, n_r1,n_r2, Lm1,Lm2, sigma)
% solve the eigenvalue problem of two-atom system with periodic potential
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain [-L/2,L/2]^3, R is the radius of sphere   
% N_p is the number of planewaves used in each direction
% n_r1 and n_r2 is the number of radius basis, Lm1 and Lm2 is the number of angular momentum
% sigma is the penalization parameter
% Aug 17th, 2018

%============================= initiate the data =========================%
% initial_data=0
% degrees of freedom
% t0=cputime;
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
n = n0 + n_r1*(Lm1+1)^2+ n_r2*(Lm2+1)^2;
fprintf('pw DOF=%d, total DOF=%d\n ',n0 ,n );
H = zeros(n,n);
M = zeros(n,n);
n_integral = 100;
delta=0.000001;
lmax = 5;
R1=[-2,0,0];              % locations of atoms
R2=[2,0,0];

% cut-off for periodic potential
N0=0;      % the number of truncated terms for perodic potential
N_cut=2*N_p;    % 2~4 times of N_p
N1=floor(N_cut);
kkk=zeros(8*N1^3,3);
for ii=-N1:N1
      m=sqrt(N_cut^2-ii^2);
      m=floor(m);
      for j=-m:m
            l=sqrt(N_cut^2-ii^2-j^2);
            l=floor(l);
            for n=-l:l
                  N0=N0+1;
                  kkk(N0,:)=[ii,j,n];
            end
      end
end

volum_omega = L^3;
volum_in = 4/3*pi*R^3;
volum_out = volum_omega-2*volum_in;

% evaluate the discrete periodic potential 
mm = 4*N+1;    % need double dofs m=2*N
vext = zeros(mm,mm,mm);
for ii = 1:mm
    for j = 1:mm
        for k = 1:mm
            sum_per=0;
            x = ii-2*N-1;
            y = j-2*N-1;
            z = k-2*N-1;
            k_index = [x,y,z];
            k_norm=norm(k_index,2)*2*pi/L;
            shift1=exp(1i*k_index*R1'*2*pi/L);
            shift2=exp(1i*k_index*R2'*2*pi/L);
                if k_norm==0   
                    for s=1:N0
                        k_per=norm(kkk(s,:),2)*2*pi/L;
                        if k_per==0
                            continue;
                        end
                        sum_per=sum_per+(2+2*cos(kkk(s,:)*(R1-R2)'*2*pi/L)) *(sin(k_per*R)/k_per^5-cos(k_per*R)*R/k_per^4) ;                           
                    end
                    v_per= sum_per*(4*pi)^2/(volum_omega)^2;
                else
                    for s=1:N0
                        KK=norm(kkk(s,:)+k_index,2)*2*pi/L;
                        k_per=norm(kkk(s,:),2)*2*pi/L;
                        if KK==0
                            continue;
                        elseif k_per==0
                            continue;
                        end
                        sum_per=sum_per+1/k_per^2*(sin(KK*R)/KK^3-R*cos(KK*R)/KK^2)*...
                        (shift1+shift2+shift2*exp(1i*kkk(s,:)*(R2-R1)'*2*pi/L)+shift1*exp(1i*kkk(s,:)*(R1-R2)'*2*pi/L));
                    end
                    v_per=4*pi/(volum_omega)^2 * ( 4*pi* sum_per+2*(shift1+shift2)*volum_in/k_norm^2);                            
                    v_per=-4*pi/(k_norm^2*volum_omega) *(shift1+shift2) +v_per;
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
        shift1=exp(1i*dk*R1'*2*pi/L);
        shift2=exp(1i*dk*R2'*2*pi/L);
        % overlap
        if p == q
            M(p,q) = volum_out/volum_omega; 
        else
            M(p,q) = -4*pi*R^2  *(shift1+shift2)* spherical_bessel(1,K*R)/(K*volum_omega);
        end
        % kinetic
        H(p,q) = 0.5 * (kk(p,:)*kk(q,:)') * M(p,q) * (2*pi/L)^2;        
        % periodic potential        
        index=dk+2*N+1;
        H(p,q)=H(p,q)+vext(index(1),index(2),index(3));       
         % discontinuous and penalization
         kp = norm(kk(p,:),2) * 2*pi/L;
         kq = norm(kk(q,:),2) * 2*pi/L;
         
         p_shift1=exp(-1i*kk(p,:)*R1'*2*pi/L);
         p_shift2=exp(-1i*kk(p,:)*R2'*2*pi/L);
         q_shift1=exp(1i*kk(q,:)*R1'*2*pi/L);
         q_shift2=exp(1i*kk(q,:)*R2'*2*pi/L);        
         % case kp=kq=0
         if kp == 0 && kq == 0
             H(p,q) = H(p,q) + 8*pi*R^2 * sigma / volum_omega;
         % case kp=0  and case kq=0
         elseif kp > 0 && kq == 0
             djp= (spherical_bessel(0,kp*(R+delta))-spherical_bessel(0,kp*(R-delta)))/(2*delta);
             jp=  spherical_bessel(0,kp*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djp+ 4*sigma * jp)*(p_shift1+p_shift2);
         elseif kp == 0 && kq > 0
             djq=  (spherical_bessel(0,kq*(R+delta))-spherical_bessel(0,kq*(R-delta)))/(2*delta);
             jq= spherical_bessel(0,kq*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djq + 4*sigma * jq)*(q_shift1+q_shift2);
         % case kp>0 and kq>0    
         elseif kp > 0 && kq > 0
             sum = 0.0;
             for l = 0:lmax
                 jp = spherical_bessel(l,kp*R);
                 jq = spherical_bessel(l,kq*R);
                 djp = (spherical_bessel(l,kp*(R+delta))-spherical_bessel(l,kp*(R-delta)))/(2*delta);
                 djq = (spherical_bessel(l,kq*(R+delta))-spherical_bessel(l,kq*(R-delta)))/(2*delta);                 
                 s_m = 0.0;
                 for m = -l:l
                     ylm1 = spherical_harmonic_xyz(l,m,-kk(p,1),-kk(p,2),-kk(p,3)); 
                     ylm2 = spherical_harmonic_xyz(l,m,kk(q,1),kk(q,2),kk(q,3)); 
                     s_m = s_m + ylm1'*ylm2; 
                 end
                 sum = sum + (-1)^l*s_m*(jp*djq +djp*jq + 4*sigma*jp*jq);
             end
             H(p,q) = H(p,q) + 0.25*(4*pi)^2*sum * R^2/volum_omega* (p_shift1*q_shift1+p_shift2*q_shift2);
         end
    end
end

% [2] ====================== matrix elements for radial type basis1 ==========================%
% calculate the integration of radial basis functions
matrix_element_r=2
g=zeros(n_r1,n_r1,N0);   % if n_r1=n_r2, g can save computational cost; if not, don't use g
int_r1 = zeros(n_r1,n_r1);
int_dr1 = zeros(n_r1,n_r1);
int1= zeros(n_r1,n_r1);
int_v1=zeros(n_r1,n_r1); 
for k1 = 1:n_r1
    for k2 = 1:n_r1
        int_r1(k1,k2) = int_overlap_radius(k1,k2,R,n_integral);
        int1(k1,k2) = int_overlap(k1,k2,R,n_integral);
        int_dr1(k1,k2) = int_overlap_radius_d(k1,k2,R,n_integral);
        % periodic potential term
        sum_per_in=0;
        for j=1:N0
            k_per=norm(kkk(j,:),2)*2*pi/L;
            if k_per==0
                continue;
            end
            g(k1,k2,j)=int_vext_per_sin(k_per,k1,k2,R,n_integral);
            sum_per_in=sum_per_in+(1+exp(1i*kkk(j,:)*(R1-R2)'*2*pi/L))*g(k1,k2,j);       
            % using scattering identity, only one component l=0 and m=0 for e^{ikr} is considered
        end            
        int_v1(k1,k2) = -4*pi*sum_per_in/volum_omega;
    end
end
% matrix elements
for k1 = 1:n_r1
    for l1 = 0:Lm1
        for m1 = -l1:l1
            for k2 = 1:n_r1
                for l2 = 0:Lm1
                    for m2 = -l2:l2
                    %-----------------------------------------------------------------------------------------------------------------%    
                        p = n0 + (k1-1)*(Lm1+1)^2 + l1^2 + m1+l1+1;
                        q = n0 + (k2-1)*(Lm1+1)^2 + l2^2 + m2+l2+1;
                        % overlap and kinetic
                        if  l1==l2 && m1==m2
                            M(p,q) = int_r1(k1,k2);     
                            kinetic_lm = l1*(l1+1) * int1(k1,k2); 
                            kinetic_r = int_dr1(k1,k2);
                            H(p,q) = kinetic_lm + kinetic_r; 
                            H(p,q) = 0.5 * H(p,q) + int_v1(k1,k2);
                        else
                            M(p,q) = 0.0;
                            H(p,q) = 0.0;
                        end                        
                        % discontinuous and penalization     
                        if l1==l2 && m1==m2
                            rp = -r_basis(k1,R,R);
                            drp = (r_basis(k1,R+delta,R)-r_basis(k1,R-delta,R))/(2*delta);
                            rq = -r_basis(k2,R,R);
                            drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta); 
                            H(p,q) = H(p,q) + 0.25*(rp*drq +drp*rq + 4*sigma*rp*rq) * R^2; 
                        end
                    end
                end
            end
        end
    end
end

% [2] ====================== matrix elements for radial type basis2 ==========================%
matrix_element_r=3
int_r2= zeros(n_r2,n_r2);
int_dr2 = zeros(n_r2,n_r2);
int2= zeros(n_r2,n_r2);
int_v2=zeros(n_r2,n_r2); 
for k1 = 1:n_r2
    for k2 = 1:n_r2
        int_r2(k1,k2) = int_overlap_radius(k1,k2,R,n_integral);
        int2(k1,k2) = int_overlap(k1,k2,R,n_integral);
        int_dr2(k1,k2) = int_overlap_radius_d(k1,k2,R,n_integral);
       % periodic potential term 
        sum_per_in=0;
        for j=1:N0
            k_per=norm(kkk(j,:),2)*2*pi/L;
            if k_per==0
                continue;
            end
            sum_per_in=sum_per_in+(1+exp(1i*kkk(j,:)*(R2-R1)'*2*pi/L))*g(k1,k2,j);   
        end            
        int_v2(k1,k2) = -4*pi*sum_per_in/volum_omega;
    end
end
% matrix elements
for k1 = 1:n_r2
    for l1 = 0:Lm2
        for m1 = -l1:l1
            for k2 = 1:n_r2
                for l2 = 0:Lm2
                    for m2 = -l2:l2
                    %-----------------------------------------------------------------------------------------------------------------%    
                        p = n0 + n_r1*(Lm1+1)^2+(k1-1)*(Lm2+1)^2 + l1^2 + m1+l1+1;
                        q = n0 + n_r1*(Lm1+1)^2+(k2-1)*(Lm2+1)^2 + l2^2 + m2+l2+1;
                        % overlap and kinetic
                        if  l1==l2 && m1==m2
                            M(p,q) = int_r2(k1,k2);      
                            kinetic_lm = l1*(l1+1) * int2(k1,k2); 
                            kinetic_r = int_dr2(k1,k2);
                            H(p,q) = kinetic_lm + kinetic_r; 
                            H(p,q) = 0.5 * H(p,q) + int_v2(k1,k2);
                        else
                            M(p,q) = 0.0;
                            H(p,q) = 0.0;
                        end                        
                  
                        % discontinuous and penalization     
                        if l1==l2 && m1==m2
                            rp = -r_basis(k1,R,R);
                            drp = (r_basis(k1,R+delta,R)-r_basis(k1,R-delta,R))/(2*delta);
                            rq = -r_basis(k2,R,R);
                            drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta); 
                            H(p,q) = H(p,q) + 0.25*(rp*drq +drp*rq + 4*sigma*rp*rq) * R^2; 
                        end
                    end
                end
            end
        end
    end
end

% [3] ====================== matrix elements for planewaves and radial type basis1 =======================%
matrix_element_pw_r =4
for p = 1:n0
    for k2 = 1:n_r1
        for l2 = 0:Lm1
            for m2 = -l2:l2
                kp = norm(kk(p,:),2) * 2*pi/L;
                shift1=exp(-1i*kk(p,:)*R1'*2*pi/L);
                q = n0 + (k2-1)*(Lm1+1)^2 + l2^2 + m2+l2+1;
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
                    jp =   spherical_bessel(l2,kp*R);
                    djp = (spherical_bessel(l2,kp*(R+delta))-spherical_bessel(l2,kp*(R-delta)))/(2*delta);
                    rq = -r_basis(k2,R,R);
                    drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);
                    H(p,q) = shift1*0.25*4*pi * 1i^l2 * ylm_part * (jp*drq +djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);
                end
                H(q,p) = H(p,q)';
            end
        end
    end
end

% [3] ====================== matrix elements for planewaves and radial type basis2 =======================%
matrix_element_pw_r =5
for p = 1:n0
    for k2 = 1:n_r2
        for l2 = 0:Lm2
            for m2 = -l2:l2
                kp = norm(kk(p,:),2) * 2*pi/L;
                shift2=exp(-1i*kk(p,:)*R2'*2*pi/L);
                q = n0 + n_r1*(Lm1+1)^2+ (k2-1)*(Lm2+1)^2 + l2^2 + m2+l2+1;
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
                    jp =   spherical_bessel(l2,kp*R);
                    djp = (spherical_bessel(l2,kp*(R+delta))-spherical_bessel(l2,kp*(R-delta)))/(2*delta);
                    rq = -r_basis(k2,R,R);
                    drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);
                    H(p,q) = shift2*0.25*4*pi * 1i^l2 * ylm_part * (jp*drq +djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);
                end
                H(q,p) = H(p,q)';
            end
        end
    end
end

%==================================== end of matrix element computation ==============================%

%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
eigen_solver=6
H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1,-100);
c=phi;
nrm=sqrt(c'*M*c);

% t=(cputime-t0)/60;
% fprintf('time= %d min\n' ,t)
save dg_2atom L R n0 N_p n_r1 n_r2 Lm1 Lm2 R1 R2 c nrm kk volum_omega delta
return;
