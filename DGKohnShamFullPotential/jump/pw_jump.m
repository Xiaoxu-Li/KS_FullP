% planewaves jump data on the surface
% Dec 15th, 2017

% fix the following parameters 
N_p=5;                   % N_p determines the size of jump_data
L=10;
R=0.5;  
sigma=10000;
delta=0.000001;
lmax = 5;

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
jump=zeros(2*N+1,2*N+1,2*N+1,2*N+1,2*N+1,2*N+1);

for p=1:n0
    for q=1:n0
         % discontinuous and penalization
         kp = norm(kk(p,:),2) * 2*pi/L;
         kq = norm(kk(q,:),2) * 2*pi/L;
          
         if kp == 0 || kq == 0
             continue;
         else
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
             index=[kk(p,1),kk(p,2),kk(p,3),kk(q,1),kk(q,2),kk(q,3)]+N+1;
             jump(index(1),index(2),index(3),index(4),index(5),index(6))=sum;
         end
    end
end
save jump_data_R_05 jump

