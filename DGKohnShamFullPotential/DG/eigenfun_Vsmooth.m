% plot for eigenfunction of one atom with V_smooth
% Aug 17th, 2018

load dg_Vsmooth
for j=0:R/0.1
    x(j+1)=0.1*j;
    y(j+1)=0;
    for n=1:n_r
        for l=0:Lm
            for m=-l:l
                if x(j+1) == 0.0
                   x(j+1) = x(j+1)+delta; 
                end
                k = n0 + (n-1)*(Lm+1)^2 + l^2 + m+l+1;
                y(j+1)=y(j+1) + c(k)/nrm*r_basis(n,x(j+1),R)*spherical_harmonic_xyz(l,m,x(j+1),0,0);
            end
        end
    end
end
for j=R/0.1+2:(L/2)/0.1+2
    x(j)=0.1*(j-2);
    y(j)=0.0;
    for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)];
        r=[x(j),0,0]+L/2;    % L/2 here because of  using FFT in the eigen_dg
        y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
    end
end
plot(x,abs(y),'r-')

hold on
for j=0:R/0.1
   x(j+1)=-0.1*j;
   y(j+1)=0;
   for n=1:n_r
       for l=0:Lm
           for m=-l:l
               if x(j+1) == 0.0
                  x(j+1) = x(j+1)-delta; 
               end
               k = n0 + (n-1)*(Lm+1)^2 + l^2 + m+l+1;
               y(j+1)=y(j+1) + c(k)/nrm*r_basis(n,-x(j+1),R)*spherical_harmonic_xyz(l,m,x(j+1),0,0);
           end
       end
   end
end
for j=R/0.1+2:(L/2)/0.1+2
   x(j)=-0.1*(j-2);
   y(j)=0.0;
   for p=1:n0
       k=[kk(p,1),kk(p,2),kk(p,3)];
       r=[x(j),0,0]+L/2; 
       y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
   end
end
plot(x,abs(y),'r-')