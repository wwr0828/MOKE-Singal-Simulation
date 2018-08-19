% calculates polar complex Kerr angle (incident angle is 0 degree)
% uses formalism of 1998_You_JAP

function [Reflectance, Phi_S, Phi_P] = MOKE_MLmodel_SUB3(lambda,h,n,Q,mz)


Af = [1                               0                    1                                  0;...
      0                               1                    0                                 -1;...
      -1i*n(end)*Q(end)*mz(end)/2   -n(end)              1i*n(end)*Q(end)*mz(end)/2     -n(end);...
      n(end)               1i*n(end)*Q(end)*mz(end)/2       -n(end)                1i*n(end)*Q(end)*mz(end)/2]; %The substrate's boundary matrix

Told = Af;%Put the substrate's the matrix as the last term

for j = length(mz)-1:-1:2 %calculate the boundary matices for magnetic layers (total is n-2 layers)
    A = [1                               0                    1                                  0;...
         0                               1                    0                                 -1;...
      -1i*n(j)*Q(j)*mz(j)/2            -n(j)              1i*n(j)*Q(j)*mz(j)/2                -n(j);...
      n(j)                   1i*n(j)*Q(j)*mz(j)/2       -n(j)                1i*n(j)*Q(j)*mz(j)/2];
     
    U = exp(-1i*2*pi/lambda*n(j)*h(j));
    deltai = -pi*n(j)*Q(j)*h(j)*mz(j)/lambda; 
    deltar = deltai;  %identical
    
    D = [U*cos(deltai)   U*sin(deltai)  0              0;...
         -U*sin(deltai)  U*cos(deltai)  0              0;...
         0               0              cos(deltar)/U  sin(deltar)/U;...
         0               0              -sin(deltar)/U  cos(deltar)/U];
     
    Tnew = A*D*A^(-1)*Told; %multiply one layer per iteration
    Told = Tnew;
end

j = 1;
Ai = [1                               0                    1                                  0;...
      0                               1                    0                                 -1;...
     -1i*n(j)*Q(j)*mz(j)/2            -n(j)              1i*n(j)*Q(j)*mz(j)/2                -n(j);...
      n(j)                   1i*n(j)*Q(j)*mz(j)/2       -n(j)                1i*n(j)*Q(j)*mz(j)/2];%Calculate the boundary matrix of air
        
T = Ai^(-1)*Told;

G = [T(1,1) T(1,2); T(2,1) T(2,2)];
I = [T(3,1) T(3,2); T(4,1) T(4,2)];

R = I/G;

Reflectance = abs(R(1,1))^2 + abs(R(2,1))^2;
Phi_S = R(2,1)/R(1,1);
Phi_P = R(1,2)/R(2,2);