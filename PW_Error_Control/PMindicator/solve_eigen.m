% solve the eigenvalue problem :
% (£­Delta+V_{ext})\phi=\lambda \phi with periodic boundary condition on [-L:L]
% as discrete version H \phi=\lambda \phi with Ec plane-wave discretization.
% Input:
% V: V_{ext}
% L: width of the domain
% Ec: energy cutoff of the planewaves
% Neig: number of eigenpairs to compute
% Eg: reference grid to calculate the residual
% Output:
% lambda: eigenvalues
% phi: eigenfunctions
function [lambda,phi] = solve_eigen(L, Ec, Neig, Eg)
% initiate the degrees of freedom
Ng = floor(L/pi*sqrt(2*Eg)); % Gg radius (to calculate the potential element)
Nc = floor(L/pi*sqrt(2*Ec)); % Gc radius
dof = 0;
dofg = 0;
for ii = -Nc:Nc
    for j = -Nc:Nc
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Nc
               dof = dof + 1;
            end
    end
end
for ii = -Ng:Ng
    for j = -Ng:Ng
            r2 = sqrt(ii^2 + j^2);
            if r2 <= Ng
               dofg = dofg + 1;
            end
    end
end
fprintf('With Ecut = %f, number of planewaves DOF = %d \n', Ec, dof);
H = zeros(dof, dof);
%HV = zeros(dof, dof);
vec_k = zeros(dof, 2);
% build up the planewave vectors
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
% initiate the discrete potential
Vext = zeros(2*Ng+1,2*Ng+1);
%V = @V_gauss_2D;
V = @V_osc_2D;
%V = @V_hydr_2D;
for ii = 1:2*Ng+1
    for j = 1:2*Ng+1
            x = -L + (ii-1)*2*L/(2*Ng+1);
            y = -L + (j-1)*2*L/(2*Ng+1);
            Vext(ii,j) = V(x,y);
     end
end

% FFT for vext
Vfft = fftn(Vext)/(2*Ng+1)^2;
% calculate the matrix element
for n1 = 1:dof
    for n2 = 1:dof
        if n1 == n2
            H(n1,n2) = norm(vec_k(n1,:),2)^2 * (pi/L)^2;
        end
        % external potential H_{pq}=vext_fft(p-q)
        dk = vec_k(n1,:) - vec_k(n2,:);
        for a = 1:2
            if dk(a) < 0
               dk(a) = dk(a) + 2*Ng + 1;
            end
            dk(a) = dk(a) + 1;
        end
        H(n1,n2) = H(n1,n2) + real(Vfft(dk(1),dk(2)));
        %HV(ii,j) = real(Vfft(dk(1),dk(2)));    
    end
end
% solve the eigenvalue problem Hu=\lamba Mu 
[phi, lambda] = eigs(H, Neig, 'SA');
lambda = diag(lambda);
return;
