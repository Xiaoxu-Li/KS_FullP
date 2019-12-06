% calculate the a posteriori indicator by dusson's formula
% we need res, laplace, V-\PiV\Pi
function [err_post] = residual(L, Ec, Neig, Eg, phi)
err_post=zeros(Neig,1);
Nc = floor(sqrt(Ec));
Ng = floor(sqrt(Eg));
Vfftdelta = zeros(2*Ng + 1, 2*Ng + 1);
phi_m = zeros(2*Ng + 1, 2*Ng + 1, Neig);
%res = zeros(2*Ng + 1, 2*Ng + 1, Neig);
for ll = 1:Neig
        n = 0;
        for ii = -Ng:Ng
            for j = -Ng:Ng
                    r2 = sqrt(ii^2 + j^2);
                    if r2 <= Nc
                        n = n + 1;
                        if ii <= 0  
                            iik = ii + 2*Ng + 1; 
                        else
                            iik = ii ;  
                        end
                        if  j <= 0  
                            jk = j + 2*Ng + 1;   
                        else
                            jk = j ;
                        end
                        Vfftdelta(iik , jk) = 0;
                        phi_m(iik, jk, ll) = phi(n, ll);
                    else
                        if ii <= 0  
                            iik = ii + 2*Ng + 1; 
                        else
                            iik = ii ;  
                        end
                        if  j <= 0  
                            jk = j + 2*Ng + 1;   
                        else
                            jk = j ;
                        end
                        Vfftdelta(iik , jk) = 20/norm([ii,j],2)^2;
                        phi_m(iik, jk, ll) = 0;
                    end
                end
        end 
        phi_m_ll = phi_m(:,:,ll);
        res(:,:,ll) = conv2(Vfftdelta, phi_m_ll);
        s = size(res(:,:,ll));
        weight = zeros(s(1), s(1));
        Ns = (s(1)-1)/2;
for ii = -Ns:Ns
    for j = -Ns:Ns
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Ns
                        if ii <= 0  
                            iik = ii + s(1); 
                        else
                            iik = ii ;  
                        end
                        if  j <= 0  
                            jk = j + s(1);   
                        else
                            jk = j ;
                        end
                weight(iik, jk) = 1.0 / (1.0 + r2^2 * (pi/L)^2);
            end
     end
end
         res_weight = reshape( real(res(:,:,ll)).^2 .* weight , s(1)^2, 1);
         err_post(ll) =err_post(ll) + norm(res_weight, 1);
end
return
        