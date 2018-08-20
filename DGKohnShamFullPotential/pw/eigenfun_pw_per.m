% plot for eigen_solve with periodic potential
% one atom or two atom 
% Aug 16th, 2018

load pw_per
for j=0:50
    xx(j+1)=0.1*j;
    yy(j+1)=0.0;
    for l=1:n
        kk=[p(l,1),p(l,2),p(l,3)];
        r=[xx(j+1),0,0]; 
        yy(j+1)=yy(j+1) + c(l)/nrm *exp(1i*kk*r'*2*pi/L)/sqrt(L^3);
    end
end
plot(xx,abs(yy),'r--')

hold on
for j=0:50
    xx(j+1)=-0.1*j;
    yy(j+1)=0.0;
    for l=1:n
        kk=[p(l,1),p(l,2),p(l,3)];
        r=[xx(j+1),0,0]; 
        yy(j+1)=yy(j+1) + c(l)/nrm *exp(1i*kk*r'*2*pi/L)/sqrt(L^3);
    end
end
plot(xx,abs(yy),'r--')