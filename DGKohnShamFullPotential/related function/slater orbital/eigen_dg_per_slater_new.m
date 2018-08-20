function [lambda,phi]=eigen_dg_per_slater_new(L, R, N_p, n_r, Lm, sigma)
% periodic potential 
% slater orbital
% 07/04/2018

% eigen_solve_per solve the eigenvalue problem    [-L/2,L/2]^3
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain, R is the radious of sphere
% N_p is the number of planewaves used in each direction
% n_r is the number of radious basis, Lm is the number of angular momentum
% sigma is the penalization parameter

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
n = n0 + n_r*(Lm+1)^2;
fprintf('pw DOF=%d, total DOF=%d\n ',n0 ,n );
H = zeros(n,n);
M = zeros(n,n);
n_integral = 100;
delta=0.000001;
lmax = 5;

% cut-off for periodic potential
N0=0;
N_cut=2*N_p;    % 2~4 times of N_p
N=floor(N_cut);
a=zeros(8*N^3,3);
for ii=-N:N
      m=sqrt(N_cut^2-ii^2);
      m=floor(m);
      for j=-m:m
            l=sqrt(N_cut^2-ii^2-j^2);
            l=floor(l);
            for n=-l:l
                  N0=N0+1;
                  a(N0,:)=[ii,j,n];
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
                        K_per=norm(a(s,:),2)*2*pi/L;
                        if K_per==0
                            continue;
                        end
                        sum_per=sum_per+sin(K_per*R)/K_per^5-cos(K_per*R)*R/K_per^4;        
                    end
                    v_per=sum_per*(4*pi)^2/(volum_omega)^2;
                else
                    for s=1:N0
                        KK=norm(a(s,:)+K_index,2)*2*pi/L;
                        K_per=norm(a(s,:),2)*2*pi/L;
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

%*********************************************************************************************************************************%
%==================================================== calculate the matrix element ===============================================%

% [1] ======================== matrix element for planewaves =======================%
matrix_pw=1
for p=1:n0
    for q=1:n0
        %-----------------------------------------------------------------%
        dk = kk(q,:)-kk(p,:);      %  conjugate change
        K = norm(dk,2) * 2*pi/L;   
        % K_sigman = exp(1i*dk*[L/2,L/2,L/2]'*2*pi/L);
        % overlap
        if p == q
            M(p,q) = volum_out/volum_omega; 
        else
            M(p,q) = -4*pi*R^2  * spherical_bessel(1,K*R)/(K*volum_omega);
        end
        % kinetic
        H(p,q) = 0.5 * (kk(p,:)*kk(q,:)') * M(p,q) * (2*pi/L)^2;         
        % external potential H_{pq}=vext_fft(p-q)
%         for k = 1:3
%             if dk(k) < 0
%                 dk(k) = dk(k) + mm;
%             end
%             dk(k) = dk(k) + 1;
%         end

         % periodic potential   
        index=dk+2*N+1;
        H(p,q)=H(p,q)+vext(index(1),index(2),index(3));
                
         % discontinuous and penalization
         kp = norm(kk(p,:),2) * 2*pi/L;
         kq = norm(kk(q,:),2) * 2*pi/L;
        % p_sigman = exp(-1i*kk(p,:)*[L/2,L/2,L/2]'*2*pi/L); % change
        % q_sigman = exp(1i*kk(q,:)*[L/2,L/2,L/2]'*2*pi/L); % change
         % case kp=kq=0
         if kp == 0 && kq == 0
             H(p,q) = H(p,q) + 4*pi*R^2 * sigma / volum_omega;
             % case kp=0  and case kq=0
         elseif kp > 0 && kq == 0
             djp= (spherical_bessel(0,kp*(R+delta))-spherical_bessel(0,kp*(R-delta)))/(2*delta);
             jp= spherical_bessel(0,kp*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djp+ 4*sigma * jp);
         elseif kp == 0 && kq > 0
             djq= (spherical_bessel(0,kq*(R+delta))-spherical_bessel(0,kq*(R-delta)))/(2*delta);
             jq=spherical_bessel(0,kq*R);
             H(p,q) = H(p,q) + 4*pi * 0.25*R^2/volum_omega * (djq + 4*sigma * jq);
             % case kp>0 and kq>0    
         elseif kp > 0 && kq > 0
%              sum = 0.0;
%              for l = 0:lmax
%                  jp = spherical_bessel(l,kp*R);
%                  jq =spherical_bessel(l,kq*R);
%                  djp = (spherical_bessel(l,kp*(R+delta))-spherical_bessel(l,kp*(R-delta)))/(2*delta);
%                  djq = (spherical_bessel(l,kq*(R+delta))-spherical_bessel(l,kq*(R-delta)))/(2*delta);                 
%                  s_m = 0.0;
%                  for m = -l:l
%                      ylm1 = spherical_harmonic_xyz(l,m,-kk(p,1),-kk(p,2),-kk(p,3)); % conjugate change
%                      ylm2 = spherical_harmonic_xyz(l,m,kk(q,1),kk(q,2),kk(q,3)); % conjugate change
%                      s_m = s_m + ylm1'*ylm2; % (real1-i*imag1)*(real2+i*imag2);
%                  end
%                  sum = sum + (-1)^l*s_m*(jp*djq +djp*jq + 4*sigma*jp*jq);
%              end
              j_index=[kk(p,1),kk(p,2),kk(p,3),kk(q,1),kk(q,2),kk(q,3)]+5+1;
             sum=jump( j_index(1),  j_index(2),  j_index(3),  j_index(4),  j_index(5),  j_index(6));                
             H(p,q) = H(p,q) + 0.25*(4*pi)^2*sum * R^2/volum_omega;
         end
         %--------------------------------------------------------------------------------------------%
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
        int_r(k1,k2) = int_overlap_radius_slater(k1,k2,R,n_integral);
        int(k1,k2) = int_overlap_slater(k1,k2,R,n_integral);
        int_dr(k1,k2) = int_overlap_radius_d_slater(k1,k2,R,n_integral);
% periodic potential
        sum_per_in=0;
        for j=1:N0
            K_per=norm(a(j,:),2)*2*pi/L;
            if K_per==0
                continue;
            end
            sum_per_in=sum_per_in+int_vext_per_sin_slater(K_per,k1,k2,R,n_integral);      
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
                            M(p,q) = int_r(k1,k2); % 4*pi/(2*l1+1)      
                            kinetic_lm = l1*(l1+1) * int(k1,k2); % *4*pi/(2*l1+1)
                            kinetic_r = int_dr(k1,k2);
                            H(p,q) = kinetic_lm + kinetic_r; % 4*pi/(2*l1+1)
                            H(p,q) = 0.5 * H(p,q) + int_v(k1,k2);
                        else
                            M(p,q) = 0.0;
                            H(p,q) = 0.0;
                        end                        
                        % external potential, only one component l=0 and m=0 for v_lm(r) is considered
                        % H(p,q) = H(p,q) + (4*pi)^2*R^4/volum_omega * ClebschGordan(l1,l2,0,m1,m2,0) * int_vext_radius(k1,k2,R,n_integral);
                        
                        % discontinuous and penalization     
                        if l1==l2 && m1==m2
                            rp = -r_basis_slater(k1,R,R);
                            drp = (r_basis_slater(k1,R+delta,R)-r_basis_slater(k1,R-delta,R))/(2*delta);
                            rq = -r_basis_slater(k2,R,R);
                            drq = (r_basis_slater(k2,R+delta,R)-r_basis_slater(k2,R-delta,R))/(2*delta); 
                            H(p,q) = H(p,q) + 0.25*(rp*drq +drp*rq + 4*sigma*rp*rq) * R^2; % 4*pi/(2*l1+1)
                        end
                    %-------------------------------------------------------------------------------------------------------------------%  
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
                % p_sigman = exp(-1i*kk(p,:)*[L/2,L/2,L/2]'*2*pi/L); % change
                q = n0 + (k2-1)*(Lm+1)^2 + l2^2 + m2+l2+1;
                % overlap is 0
                M(p,q) = 0.0;
                M(q,p) = M(p,q);                        
                % only discontinuous and penalization is considered for H(p,q)
                if kp == 0
                    if l2 == 0
                        rq = -r_basis_slater(k2,R,R);
                        drq = (r_basis_slater(k2,R+delta,R)-r_basis_slater(k2,R-delta,R))/(2*delta);
                        jp = 1;
                        H(p,q) = 0.25*R^2 * sqrt(4*pi)/sqrt(volum_omega) * (jp*drq + 4*sigma*rq*jp);
                    else
                        H(p,q) = 0.0;
                    end
                else
                    ylm_part = spherical_harmonic_xyz(l2,m2,-kk(p,1),-kk(p,2),-kk(p,3));  %  conjugate change
                    jp =  spherical_bessel(l2,kp*R);
                    djp = (spherical_bessel(l2,kp*(R+delta))-spherical_bessel(l2,kp*(R-delta)))/(2*delta);
                    rq = -r_basis_slater(k2,R,R);
                    drq = (r_basis_slater(k2,R+delta,R)-r_basis_slater(k2,R-delta,R))/(2*delta);
                    H(p,q) = 0.25*4*pi * 1i^l2 * ylm_part * (jp*drq +djp*rq + 4*sigma*jp*rq) * R^2/sqrt(volum_omega);  % change
                end
                H(q,p) = H(p,q)';
                %-------------------------------------------------------------------------------------------------------------------%  
            end
        end
    end
end
%==================================== end of matrix element computation ==============================%
%********************************************************************************************************************************%

%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
eigen_solver=4
% lambda = eigs(H, M, 1, -10.0);

H=(H+H')/2;
[phi,lambda] = eigs(H, M, 1,  -100);
c=phi ;
nrm=sqrt(c'*M*c);
% 
% for j=0:R/0.1
%     x(j+1)=0.1*j;
%     y(j+1)=0;
%     for n=1:n_r
%         for l=0:Lm
%             for m=-l:l
%                 if x(j+1) == 0.0
%                    x(j+1) = x(j+1)+delta; 
%                 end
%                 k = n0 + (n-1)*(Lm+1)^2 + l^2 + m+l+1;
%                 y(j+1)=y(j+1) + c(k)/nrm*r_basis(n,x(j+1),R)*spherical_harmonic_xyz(l,m,x(j+1),0,0);
%             end
%         end
%     end
% end
% for j=R/0.1+2:(L/2)/0.1+2
%     x(j)=0.1*(j-2);
%     y(j)=0.0;
%     for p=1:n0
%         k=[kk(p,1),kk(p,2),kk(p,3)];
%         r=[x(j),0,0]; 
%         y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
%     end
% end
% figure(1)
% plot(x,abs(y),'b-')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% for j=0:R/0.1
%    x(j+1)=-0.1*j;
%    y(j+1)=0;
%    for n=1:n_r
%        for l=0:Lm
%            for m=-l:l
%                if x(j+1) == 0.0
%                   x(j+1) = x(j+1)-delta; 
%                end
%                k = n0 + (n-1)*(Lm+1)^2 + l^2 + m+l+1;
%                y(j+1)=y(j+1) + c(k)/nrm*r_basis(n,-x(j+1),R)*spherical_harmonic_xyz(l,m,x(j+1),0,0);
%            end
%        end
%    end
% end
% for j=R/0.1+2:(L/2)/0.1+2
%    x(j)=-0.1*(j-2);
%    y(j)=0.0;
%    for p=1:n0
%        k=[kk(p,1),kk(p,2),kk(p,3)];
%        r=[x(j),0,0]; 
%        y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
%    end
% end
% plot(x,abs(y),'b-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t=(cputime-t0)/60;
% fprintf('time=%d min\n',t);
% save MAT L R n0 N_p n_r Lm c nrm kk volum_omega delta
return;





%**************** functions related to the external potential **************%
% integration y=1/k*\int_0^R r*v(r)*sin(kr) dr
% function y = int_vext_sink(k,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% if k == 0
%     for i=1:N
%         xl = (i-1)*h;
%         xr = i*h;
%         s = s + 0.5 * h * (xl*V_smooth_r(xl,4.5) + xr*V_smooth_r(xr,4.5));
%     end
% else
%     for i=1:N
%         xl = (i-1)*h;
%         xr = i*h;
%         s = s + 0.5 * h * (V_smooth_r(xl,4.5)*sin(k*xl) + V_smooth_r(xr,4.5)*sin(k*xr));
%     end
%     s = s/k;
% end
% y=s;
% return;
% 
% 
% %******************* functions about integrations related to radial basis****************%
% % radial basis functions \chi_i
% function y = r_basis(i,r,R)
% y = r^(i-1) *i/(R^i);
% return;
% 
% % derivative of radial basis functions \chi'_i
% function y = r_basis_d(i,r,R)
% if i == 1
%     y = 0;
% else
%     y = (i-1)*r^(i-2) *i/(R^i);
% end
% return;
% 
% % integration y=\int_0^R r^2*chi_i(r)*chi_j(r) dr
% function y = int_overlap_radius(i,j,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% for k=1:N
%     xl = (k-1)*h;
%     xr = k*h;
%     s = s + 0.5 * h * (xl^2*r_basis(i,xl,R)*r_basis(j,xl,R) +xr^2* r_basis(i,xr,R)*r_basis(j,xr,R));
% end
% y=s;
% return;
% 
% % integration y=\int_0^R chi_i(r)*chi_j(r) dr
% function y = int_overlap(i,j,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% for k=1:N
%     xl = (k-1)*h;
%     xr = k*h;
%     s = s + 0.5 * h * (r_basis(i,xl,R)*r_basis(j,xl,R) + r_basis(i,xr,R)*r_basis(j,xr,R));
% end
% y=s;
% return;
% 
% % integration y=\int_0^R r^2*chi'_i(r)*chi'_j(r) dr
% function y = int_overlap_radius_d(i,j,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% for k=1:N
%     xl = (k-1)*h;
%     xr = k*h;
%     s = s + 0.5 * h * (xl^2*r_basis_d(i,xl,R)*r_basis_d(j,xl,R) + xr^2*r_basis_d(i,xr,R)*r_basis_d(j,xr,R));
% end
% y=s;
% return;
%     
% 
% % integration y=\int_0^R r^2*v_lm(r)*chi_i(r)*chi_j(r) dr
% function y = int_vext_radius(i,j,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% for k=1:N
%     xl = (k-1)*h;
%     xr = k*h;
%     % s = s + 0.5 * h * (- xl*r_basis(i,xl,R)*r_basis(j,xl,R) - xr*r_basis(i,xr,R)*r_basis(j,xr,R));
%     s = s + 0.5 * h * (xl*V_smooth_r(xl,4.5)*r_basis(i,xl,R)*r_basis(j,xl,R) + xr*V_smooth_r(xr,4.5)*r_basis(i,xr,R)*r_basis(j,xr,R));
% end
% y=s;
% return;    
% 
% % integration y=1/k^3*\int_0^R r*sin(kr)*chi_i(r)*chi_j(r) dr
% function y=int_vext_per_sin(K,i,j,R,n_integral)
% N = n_integral; % integration discretization
% h = R/N;
% s = 0;
% for k=1:N
%     xl = (k-1)*h;
%     xr = k*h;
%     s = s + 0.5 * h * (xl*sin(K*xl)*r_basis(i,xl,R)*r_basis(j,xl,R) + xr*sin(K*xr)*r_basis(i,xr,R)*r_basis(j,xr,R));
% end
% y=s/K^3 ;
% return;  
