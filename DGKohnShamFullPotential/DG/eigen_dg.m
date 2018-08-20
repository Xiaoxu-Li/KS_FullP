  function [lambda,phi]=eigen_dg(L, R, N_p, n_r, Lm, sigma)
% eigen_dg solve the eigenvalue problem with V_smooth
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain [-L/2,L/2]^3, R is the radious of sphere
% N_p is the number of planewaves used in each direction
% n_r is the number of radious basis, Lm is the number of angular momentum
% sigma is the penalization parameter
% 16th Aug, 2018

%============================= initiate the data =========================%
% initial_data=0
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
n = n0 + n_r*(Lm+1)^2;
fprintf('DOF=%d \n',n);
H = zeros(n,n);
M = zeros(n,n);
n_integral = 100;
delta=0.000001;
lmax = 5;

% evaluate the discrete potential
mm = 4*N_p+1;
% numerical integration for FFT needs double dofs m=2*N
vext = zeros(mm,mm,mm);
for ii = 1:mm
    for j = 1:mm
        for k = 1:mm
            x = (ii-1)*L/mm;
            y = (j-1)*L/mm;
            z = (k-1)*L/mm;
            r = ((x-L/2)^2 + (y-L/2)^2 + (z-L/2)^2)^(1/2);
            vext(ii,j,k) = V_smooth(r, 4.5);
        end
    end
end
% FFT for vext
vext_fft = fftn(vext)/(mm^3);

volum_omega = L^3;
volum_in = 4/3*pi*R^3;
volum_out = volum_omega-volum_in;
%==================================================== calculate the matrix element ===============================================%

% [1] ======================== matrix element for planewaves =======================%
matrix_pw=1
for p=1:n0
    for q=1:n0
        dk = kk(q,:)-kk(p,:);      
        K = norm(dk,2) * 2*pi/L;   
        K_sign = exp(1i*dk*[L/2,L/2,L/2]'*2*pi/L);
        % overlap
        if p == q
            M(p,q) = volum_out/volum_omega; 
        else
            M(p,q) = -4*pi*R^2 * K_sign * spherical_bessel(1,K*R)/(K*volum_omega);
        end
        % kinetic
        H(p,q) = 0.5 * (kk(p,:)*kk(q,:)') * M(p,q) * (2*pi/L)^2;         
        % external potential H_{pq}=vext_fft(p-q)
        for k = 1:3
            if dk(k) < 0
                dk(k) = dk(k) + mm;
            end
            dk(k) = dk(k) + 1;
        end
        H(p,q) = H(p,q) + vext_fft(dk(1),dk(2),dk(3)) - 4*pi*K_sign* int_vext_sink(K,R,n_integral)/volum_omega;
         % discontinuous and penalization
         kp = norm(kk(p,:),2) * 2*pi/L;
         kq = norm(kk(q,:),2) * 2*pi/L;
         p_sign = exp(-1i*kk(p,:)*[L/2,L/2,L/2]'*2*pi/L); 
         q_sign = exp(1i*kk(q,:)*[L/2,L/2,L/2]'*2*pi/L); 
         % case kp=kq=0
         if kp == 0 && kq == 0
             H(p,q) = H(p,q) + 4*pi*R^2 * sigma / volum_omega;
         % case kp=0  and case kq=0
         elseif kp > 0 && kq == 0
             djp= p_sign *(spherical_bessel(0,kp*(R+delta))-spherical_bessel(0,kp*(R-delta)))/(2*delta);
             jp= p_sign * spherical_bessel(0,kp*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djp+ 4*sigma * jp);
         elseif kp == 0 && kq > 0
             djq= q_sign * (spherical_bessel(0,kq*(R+delta))-spherical_bessel(0,kq*(R-delta)))/(2*delta);
             jq=q_sign * spherical_bessel(0,kq*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djq + 4*sigma * jq);
             % case kp>0 and kq>0    
         elseif kp > 0 && kq > 0
             sum = 0.0;
             for l = 0:lmax
                 jp =p_sign * spherical_bessel(l,kp*R);
                 jq =q_sign * spherical_bessel(l,kq*R);
                 djp =p_sign * (spherical_bessel(l,kp*(R+delta))-spherical_bessel(l,kp*(R-delta)))/(2*delta);
                 djq =q_sign * (spherical_bessel(l,kq*(R+delta))-spherical_bessel(l,kq*(R-delta)))/(2*delta);                 
                 s_m = 0.0;
                 for m = -l:l
                     ylm1 = spherical_harmonic_xyz(l,m,-kk(p,1),-kk(p,2),-kk(p,3)); 
                     ylm2 = spherical_harmonic_xyz(l,m,kk(q,1),kk(q,2),kk(q,3)); 
                     s_m = s_m + ylm1'*ylm2;
                 end
                 sum = sum + (-1)^l*s_m*(jp*djq +djp*jq + 4*sigma*jp*jq);
             end
             H(p,q) = H(p,q) + 0.25*(4*pi)^2*sum * R^2/volum_omega;
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
        int_v(k1,k2) = int_vext_radius(k1,k2,R,n_integral);
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
                            H(p,q) = 0.5 * H(p,q) + int_v(k1,k2);
                        else
                            M(p,q) = 0.0;
                            H(p,q) = 0.0;
                        end                        
                        % external potential, only one component l=0 and m=0 for v_lm(r) is considered
                        % H(p,q) = H(p,q) + (4*pi)^2*R^4/volum_omega * ClebschGordan(l1,l2,0,m1,m2,0) * int_vext_radius(k1,k2,R,n_integral);
                        
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

% [3] ====================== matrix elements for planewaves and radial type basis =======================%
matrix_element_pw_r =3
for p = 1:n0
    for k2 = 1:n_r
        for l2 = 0:Lm
            for m2 = -l2:l2
                %-------------------------------------------------------------------------------------------------------------------%
                kp = norm(kk(p,:),2) * 2*pi/L;
                p_sign = exp(-1i*kk(p,:)*[L/2,L/2,L/2]'*2*pi/L); 
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
                    jp =  p_sign * spherical_bessel(l2,kp*R);
                    djp = p_sign *(spherical_bessel(l2,kp*(R+delta))-spherical_bessel(l2,kp*(R-delta)))/(2*delta);
                    rq = -r_basis(k2,R,R);
                    drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);
                    H(p,q) = 0.25*4*pi * 1i^l2 * ylm_part * (jp*drq +djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);  
                end
                H(q,p) = H(p,q)';
            end
        end
    end
end
%==================================== end of matrix element computation ==============================%
%********************************************************************************************************************************%

%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
eigen_solver=4
H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1, -100);
c=phi;
nrm=sqrt(c'*M*c);
save dg_Vsmooth L R n0 N_p n_r Lm c nrm kk volum_omega delta
% t=(cputime-t0)/60;
% fprintf('time=%d min\n',t)

return;
