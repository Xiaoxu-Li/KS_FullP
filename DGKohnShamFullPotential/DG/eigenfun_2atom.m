% plot for eigenfunction of two-atom system
% Aug 17th, 2018

load dg_2atom
r1=norm(R2,2)-R;   % guarantee r1, r3>0
r2=norm(R2,2)+R;
r3=norm(R1,2)-R;
r4=norm(R1,2)+R;

for j=0:r1/0.1
    x(j+1)=0.1*j; 
    y(j+1)=0;
    for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)];
        r=[x(j+1),0,0]; 
        y(j+1)=y(j+1) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
    end
end

for j=r1/0.1+2:r2/0.1+2
    x(j)=0.1*(j-2);
    y(j)=0.0;
    r=x(j)-R2(1);
    for n=1:n_r2
        for l=0:Lm2
            for m=-l:l
                if  r== 0.0
                   x(j) = x(j)+delta; 
                    r=x(j)-R2(1);
                end
                k = n0 +n_r1*(Lm1+1)^2+ (n-1)*(Lm2+1)^2 + l^2 + m+l+1;
                y(j)=y(j) + c(k)/nrm*r_basis(n,abs(r),R)*spherical_harmonic_xyz(l,m,r, 0, 0);
            end
        end
    end
end

for j=r2/0.1+3:(L/2)/0.1+3
    x(j)=0.1*(j-3);
    y(j)=0.0;
    for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)];
        r=[x(j),0,0]; 
        y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
    end
end
plot(x,abs(y),'bx-')

hold on
clear x y
for j=0:r3/0.1
   x(j+1)=-0.1*j;
   y(j+1)=0;
   for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)];
        r=[x(j+1),0,0]; 
        y(j+1)=y(j+1) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
    end
end

for j=r3/0.1+2:r4/0.1+2
   x(j)=-0.1*(j-2);
   y(j)=0.0;
   r=x(j)-R1(1);
   for n=1:n_r1
       for l=0:Lm1
           for m=-l:l
               if  r == 0.0
                   x(j) = x(j)+delta; 
                   r=x(j)-R1(1);
                end
               k = n0 + (n-1)*(Lm1+1)^2 + l^2 + m+l+1;
               y(j)=y(j) + c(k)/nrm*r_basis(n,abs(r),R)*spherical_harmonic_xyz(l,m,r, 0, 0);
           end
       end
   end
end

for j=r4/0.1+3:(L/2)/0.1+3
   x(j)=-0.1*(j-3);
   y(j)=0;
   for p=1:n0
        k=[kk(p,1),kk(p,2),kk(p,3)];
        r=[x(j),0,0]; 
        y(j)=y(j) + c(p)/nrm*exp(1i*k*r'*2*pi/L)/sqrt(volum_omega);
    end
end
plot(x,abs(y),'bx-')