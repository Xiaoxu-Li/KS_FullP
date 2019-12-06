% Use Dusson's formula
function [err_post,err_post2] = Dresidualcrr(L, Ec, Neig, Eg, phi)
err_post = zeros(Neig,1);
err_post2 = zeros(Neig,1);
Nc = floor(sqrt(Ec));
Ng = floor(sqrt(Eg));
Vfftdelta = zeros(2*Ng + 1, 2*Ng + 1);
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
HV = zeros(dof,dof);
vec_k = zeros(dof, 2);
% build up the planewave vectors (insert in initiation?)
n = 0; 
for ii = -Nc:Nc
    for j = -Nc:Nc
            r2 = sqrt(ii^2 + j^2 );
            if r2 <= Nc
                n=n+1;
            vec_k(n,:) = [ii, j];   
            end
    end
end
%V = @V_gauss_2D;
%V = @V_osc_2D;
V = @V_hydr_2D;
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
    end
end
Vn_phi = HV * phi;
res = zeros(2*Ng + 1, 2*Ng + 1, Neig);
for ll = 1:Neig
        n = 0;
        for ii = -Ng:Ng
            for j = -Ng:Ng
                    r2 = sqrt(ii^2 + j^2);
                    if r2 <= Nc
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
                        Vfft(iik,jk) = 0;
                    end

                end
        end 
        u(:,:,ll)=ifftn(phi_m(:,:,ll))*((2*Ng + 1)^2)/(2*L);
        Vu(:,:,ll)=ifftn(Vphi_m(:,:,ll))*((2*Ng + 1)^2)/(2*L);
        deltaVphi = Vext .* u(:,:,ll) - Vu(:,:,ll);
        res(:,:,ll) = fftn(deltaVphi)*(2*L)/((2*Ng + 1)^2);
        s = size(res(:,:,ll));
        weight = zeros(s(1), s(1));
        Ns = (s(1)-1)/2;
for ii = -Ns:Ns
    for j = -Ns:Ns
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Ns
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
                weight(iik, jk) = 1.0 / (1.0 + r2^2 * (pi/L)^2);
            end
     end
end
         invlapres_fft = real(res(:,:,ll)) .* weight;
         res_weight = reshape(real(invlapres_fft.*res(:,:,ll)), s(1)^2, 1);
         trunvext = ifftn(Vfft)*(s(1)^2)/(2*L);
         invlapres_fft = real(res(:,:,ll)) .* weight;
         invlapres_ifft = ifftn(invlapres_fft)*(s(1)^2)/(2*L);
         Vinvlapres_ifft = trunvext .* invlapres_ifft ;
         Vinvlapres_fft = fftn(Vinvlapres_ifft)*(2*L)/(s(1)^2);
         invlapVinvlapres_fft = real(Vinvlapres_fft) .* weight;
         res_weight2 = reshape(real(invlapVinvlapres_fft.*res(:,:,ll)), s(1)^2, 1);
         err_post(ll) =err_post(ll) + norm(res_weight, 1);
         err_post2(ll) =err_post2(ll) + norm(res_weight, 1) - norm(res_weight2, 1);
end
return
