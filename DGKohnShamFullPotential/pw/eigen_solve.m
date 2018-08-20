function [lambda,phi]=eigen_solve(L, Nc)
% eigen_solve solve the eigenvalue problem
% (-1/2\Delta+V_{ext})u=\lambda u [\lambda,\phi] returns the eigenpair
% L is the width of the domain     [-L/2,L/2]^3
% Nc is the number of planewaves used in each direction.
 

%============================= initiate the data =========================%
n=0; % number of basis satisfying |k|<Nc
N=floor(Nc);
p=zeros(8*N^3,3);
for ii=-N:N
    m=sqrt(Nc^2-ii^2);
    m=floor(m);
    for j=-m:m
        l=sqrt(Nc^2-ii^2-j^2);
        l=floor(l);
        for k=-l:l
            n=n+1;
            p(n,:)=[ii,j,k];
        end
    end
end 
fprintf('DOF=%d\n', n)
H = zeros(n,n);
M = zeros(n,n);

% evaluate the discrete potential
m = 4*N+1;  % numerical integration for FFT needs double dofs m=2*N
vext = zeros(m,m,m);
for ii = 1:m
    for j = 1:m
        for k = 1:m
           x = (ii-1)*L/m;
            y = (j-1)*L/m;
            z = (k-1)*L/m;
            r = ((x-L/2)^2 + (y-L/2)^2 + (z-L/2)^2)^(1/2);
            %vext(ii,j,k) = V_hydrogen(r);            
            vext(ii,j,k) = V_smooth(r,4.5);
           % vext(ii,j,k) = V_osc(r);
        end
    end
end
% FFT for vext
vext_fft = fftn(vext)/(m^3);

%======================== calculate the matrix element ===================%
for ii=1:n
    for j=1:n
        if p(ii,:) == p(j,:)
            H(ii,j) = 0.5*norm(p(ii,:),2)^2 *(2*pi/L)^2;
            M(ii,j) = 1.0;
        end
        % external potential H_{pq}=vext_fft(p-q)
       %  dk = p(ii,:)-p(j,:);
       dk = p(j,:)-p(ii,:);
        for k = 1:3
            if dk(k) < 0
                dk(k) = dk(k) + m;
            end
            dk(k) = dk(k) + 1;
        end
        H(ii,j) = H(ii,j) + vext_fft(dk(1),dk(2),dk(3));        
    end
end
%=============== solve the eigenvalue problem Hu=\lamba Mu ===============%
[phi,lambda] = eigs(H, M, 1, -100);
c=phi;
nrm=sqrt(c'*M*c);
save pw_Vsmooth L  c nrm p n

