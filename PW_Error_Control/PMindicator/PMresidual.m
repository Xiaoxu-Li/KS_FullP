% Using Perturbation-based-method Calculate Error Indicator
function [err_post,err_post2,res,invlapres_fftk] = PMresidual(L, Ec, Neig, Eg, phi)
err_post = zeros(Neig,1);
err_post2 = zeros(Neig,1);
Ng = floor(L/pi*sqrt(2*Eg)); % Gg radius
Nc = floor(L/pi*sqrt(2*Ec)); % Gc radius
phi_m = zeros(2*Ng + 1, 2*Ng + 1, Neig);
u = zeros(2*Ng + 1, 2*Ng + 1, Neig);
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
% build up the planewave vectors
vec_k = zeros(dof, 2);
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
vec_km = zeros(dofg, 2);  
m = 0; 
for ii = -Ng:Ng
    for j = -Ng:Ng
            r2 = sqrt(ii^2 + j^2 );
            if r2 <= Ng
                m = m + 1;
            vec_km(m,:) = [ii, j];   
            end
    end
end
% initiate the discrete potential
%V = @V_gauss_2D;
V = @V_osc_2D;
%V = @V_hydr_2D;
Vext = zeros(2*Ng+1,2*Ng+1);
for ii = 1:2*Ng+1
    for j = 1:2*Ng+1
            x = -L + (ii-1)*2*L/(2*Ng+1);
            y = -L + (j-1)*2*L/(2*Ng+1);
            Vext(ii,j) = V(x,y);
    end
end
% FFT for vext
Vfft = fftn(Vext)/(2*Ng+1)^2;
HV = zeros(dof,dof);
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

res = zeros(2*Ng + 1, 2*Ng + 1, Neig);
for ll = 1:Neig
        n = 0;
        for ii = -Nc:Nc
            for j = -Nc:Nc
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
                    end
            end
        end 
        u(:,:,ll)=ifftn(phi_m(:,:,ll))*((2*Ng + 1)^2)/(2*L);
        deltaVphi = Vext .* u(:,:,ll);
        
        res(:,:,ll) = (fftn(deltaVphi)*(2*L))/((2*Ng + 1)^2);
        n = 0;
        for ii = -Nc:Nc
            for j = -Nc:Nc
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
                        res(iik,jk,ll) = 0;
                    end
            end
        end        
weight = zeros(2*Ng + 1, 2*Ng + 1);
for ii = -Ng:Ng
    for j = -Ng:Ng
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Ng
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
                if r2~=0
                weight(iik, jk) = 1.0 / (r2^2 * (pi/L)^2);
                else     
                weight(iik, jk) = 0;    
                end
            end
     end
end

Hm = zeros(dofg,dofg);
for n1 = 1:dofg
    for n2 = 1:dofg
        dk = vec_km(n1,:) - vec_km(n2,:);
        for a = 1:2
            if dk(a) < 0
                dk(a) = dk(a) + 2*Ng +1;
            end
            dk(a) = dk(a) + 1;
        end
        Hm(n1,n2) = Hm(n1,n2) + real(Vfft(dk(1),dk(2)));  
    end
end
         invlapres_fft = real(res(:,:,ll)) .* weight;
         res_weight = reshape(real(invlapres_fft.*res(:,:,ll)), (2*Ng + 1)^2, 1);
% back to reciprocal Space
        invlapres_fftk = zeros(dofg,1);
        m=0;
        for ii=-Ng:Ng
            for j = -Ng:Ng
                r2 = sqrt(ii^2 + j^2);
                if r2 <= Ng
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
                end
            end
        end
        Vinvlapres = Hm * invlapres_fftk;
        res_weight2 = invlapres_fftk'*Vinvlapres;
        err_post(ll) = err_post(ll) + norm(res_weight, 1);
        err_post2(ll) = err_post2(ll) + norm(res_weight, 1) - norm(res_weight2,1);
end
return
