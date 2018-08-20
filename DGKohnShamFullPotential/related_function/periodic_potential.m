% plot for periodic potential and -1/r
% May 31st, 2018

clear;
n0=0;
N_p=5;        % cutoff
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
fprintf('cutoff=%d\n ',n0  );

L=10;
x=-L/2:0.1:L/2;
y=zeros(size(x));
for k=1:101
    r=[x(k),0,0];
    sum=0;
    for j=1:n0
        K=norm(kk(j,:),2)*2*pi/L;
         if K==0
            continue;
         end
        sum=sum+exp(1i*kk(j,:)*r'*2*pi/L)/K^2;
    end
    y(k)=-4*pi/L^3*sum;
end
yy=-1./abs(x);
plot(x,real(y),'b-',x,yy,'r:')