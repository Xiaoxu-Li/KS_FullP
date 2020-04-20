function [u] = backrealnew(L, Ec, Eg, phi)
ll = 1;
Ng = floor((L/pi)*sqrt(2*Eg)); 
Nc = floor((L/pi)*sqrt(2*Ec)); 
phi_m = zeros(2*Ng + 1, 2*Ng + 1);
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
                    end
            end
        end
      end
      u = ifftn(phi_m)*((2*Ng + 1)^2)/(2*L);
      return