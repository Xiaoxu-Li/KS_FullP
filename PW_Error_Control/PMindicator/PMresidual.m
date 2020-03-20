% Perturbation-based method
function [err_post,err_post2,err_post12, invlapres_fftk,Hv,Hpiv,res] = PMresidual(L, Ec, Neig, Eg, phi)
err_post = zeros(Neig,1);
err_post2 = zeros(Neig,1);
err_post12 = zeros(Neig,1);
Ng = floor((L/pi)*sqrt(2*Eg)); 
Nc = floor((L/pi)*sqrt(2*Ec)); 
Vext = zeros(2*Ng + 1, 2*Ng + 1);
phi_m = zeros(2*Ng + 1, 2*Ng + 1, Neig);
Vphi_m = zeros(2*Ng + 1, 2*Ng + 1, Neig);
u = zeros(2*Ng + 1, 2*Ng + 1, Neig);
Vu = zeros(2*Ng + 1, 2*Ng + 1, Neig);
% initiate the degrees of freedom
dof = 0;
for ii = -Nc:Nc
    for j = -Nc:Nc
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Nc
               dof = dof + 1;
            end
    end
end
dofg = 0;
for ii = -Ng:Ng
    for j = -Ng:Ng
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Ng
               dofg = dofg + 1;
            end
    end
end
HV = zeros(dof,dof);
Hpiv = zeros(dofg,dofg);
vec_k = zeros(dof, 2);
% build up the planewave vectors
n = 0; 
for R = 0:Nc
 for ii = -Nc:Nc
    for j = -Nc:Nc
        r2 = sqrt(ii^2 + j^2 );
        if R-1 < r2 && r2 <= R
           n = n+1;
           vec_k(n,:) = [ii, j];
        end
    end
 end
end
V = @V_gauss_2D;
for ii = 1:2*Ng+1
    for j = 1:2*Ng+1
            x = -L + (ii-1)*2*L/(2*Ng+1);
            y = -L + (j-1)*2*L/(2*Ng+1);
            Vext(ii,j) = V(x,y);
    end
end
% FFT for vext
Vfft = fftn(Vext)/(2*Ng+1)^2;
for n1 = 1:dof
    for n2 = 1:dof
        % external potential H_{pq}=vext_fft(p-q)
        dk = vec_k(n1,:) - vec_k(n2,:);
        for a = 1:2
            if dk(a) < 0
               dk(a) = dk(a) + 2*Ng + 1;
            end
            dk(a) = dk(a) + 1;
        end
        HV(n1,n2) = real(Vfft(dk(1),dk(2)));  
        Hpiv(n1,n2) = real(Vfft(dk(1),dk(2)));  
    end
end
Vn_phi = HV * phi;
res = zeros(2*Ng + 1, 2*Ng + 1, Neig);
for ll = 1:Neig
    % build up the planewave vectors
n = 0; 
for R = 0:Nc
 for ii = -Nc:Nc
    for j = -Nc:Nc
        r2 = sqrt(ii^2 + j^2 );
        if R-1 < r2 && r2 <= R
           n = n+1;
           vec_k(n,:) = [ii, j];
        end
    end
 end
end
n = 0;
      for R = 0:Nc
        for ii = -Ng:Ng
            for j = -Ng:Ng
                    r2 = sqrt(ii^2 + j^2);
                    if R-1 < r2 && r2 <= R
                        n = n + 1;
                        if ii < 0  
                            iik = ii + 2*Ng + 1; 
                        else
                            iik = ii ;  
                        end
                        iik = iik + 1;
                        if  j < 0  
                            jk = j + 2*Ng + 1;   
                        else
                            jk = j ;
                        end
                        jk = jk + 1;
                        phi_m(iik, jk, ll) = phi(n, ll);
                        Vphi_m(iik, jk, ll) = Vn_phi(n, ll);
                    end
            end
        end
      end
        u(:,:,ll)=ifftn(phi_m(:,:,ll))*((2*Ng + 1)^2)/(2*L);
        Vu(:,:,ll)=ifftn(Vphi_m(:,:,ll))*((2*Ng + 1)^2)/(2*L);
        deltaVphi = Vext .* u(:,:,ll) - Vu(:,:,ll);
        res(:,:,ll) = fftn(deltaVphi)*(2*L)/((2*Ng + 1)^2);
        n = 0;
      for R = 0:Nc
        for ii = -Ng:Ng
            for j = -Ng:Ng
                    r2 = sqrt(ii^2 + j^2);
                    if R-1 < r2 && r2 <= R
                        n = n + 1;
                        if ii < 0  
                            iik = ii + 2*Ng + 1; 
                        else
                            iik = ii ;  
                        end
                        iik = iik + 1;
                        if  j < 0  
                            jk = j + 2*Ng + 1;   
                        else
                            jk = j ;
                        end
                        jk = jk + 1;
                        res(iik,jk,ll) = 0;
                    end
            end
        end
      end
        s = size(res(:,:,ll));
        weight = zeros(s(1), s(1));
        weight2 = zeros(s(1), s(1));
        Ns = (s(1)-1)/2;
for R = 0:Ns
 for ii = -Ns:Ns
    for j = -Ns:Ns
            r2 = sqrt(ii^2 + j^2);
            if R-1 < r2 && r2 <= R
                        if ii < 0  
                            iik = ii + s(1); 
                        else
                            iik = ii ;  
                        end
                        iik = iik + 1;
                        if  j < 0  
                            jk = j + s(1);   
                        else
                            jk = j ;
                        end
                        jk = jk + 1;
                weight(iik, jk) = 1.0 / (real(Vfft(1,1)) + r2^2 * (pi/L)^2);
                if r2 ~= 0
                weight2(iik, jk) = 1.0 / (r2^2 * (pi/L)^2);
                else
                weight2(iik, jk) = 0;
                end
            end
    end
 end
end
         invlapres_fft = real(res(:,:,ll)) .* weight;
         invlapres_fft2 = real(res(:,:,ll)) .* weight2;
         res_weight = reshape(real(invlapres_fft.*res(:,:,ll)), s(1)^2, 1);
         res_weight12 = reshape(real(invlapres_fft2.*res(:,:,ll)), s(1)^2, 1);
         invlapres_fftk = zeros(dofg,1);
         invlapres_fftk2 = zeros(dofg,1);
      m=0;
      for R = 0:Ng
        for ii = -Ng:Ng
            for j = -Ng:Ng
                r2 = sqrt(ii^2 + j^2);
                if R-1 < r2 && r2 <= R
                        m = m + 1;
                        if ii < 0  
                            iik = ii + 2*Ng + 1; 
                        else
                            iik = ii ;  
                        end
                        iik = iik + 1;
                        if  j < 0  
                            jk = j + 2*Ng + 1;   
                        else
                            jk = j ;
                        end
                        jk = jk + 1;
                        invlapres_fftk(m) = invlapres_fft(iik, jk);
                        invlapres_fftk2(m) = invlapres_fft2(iik, jk);
                end
            end
        end
      end
vec_km = zeros(dofg, 2);  
m = 0; 
for R = 0:Ng
 for ii = -Ng:Ng
    for j = -Ng:Ng
            r2 = sqrt(ii^2 + j^2 );
            if R-1 < r2 && r2 <= R
                m = m + 1;
            vec_km(m,:) = [ii, j];   
            end
    end
 end
end
Hv = zeros(dofg,dofg);
for n1 = 1:dofg
    for n2 = 1:dofg
        dk = vec_km(n1,:) - vec_km(n2,:);
        for a = 1:2
            if dk(a) < 0
               dk(a) = dk(a) + 2*Ng +1;
            end
            dk(a) = dk(a) + 1;
        end
        Hv(n1,n2) = real(Vfft(dk(1),dk(2)));  
    end
end
         Vinvlapres = (Hv - Hpiv - diag(diag(Hv - Hpiv))) * invlapres_fftk;
         res_weight2 = invlapres_fftk'*Vinvlapres;
         err_post(ll) =err_post(ll) + norm(res_weight, 1); % 1st
         err_post12(ll) =err_post12(ll) + norm(res_weight12, 1); % old 1st
         err_post2(ll) = err_post(ll) - real(res_weight2); % 2nd
end
return