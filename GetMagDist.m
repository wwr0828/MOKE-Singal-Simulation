%Input the spin chemical potential distribution from step 1, and output the
%SOT induced magnetization distribution in a singal layer of FM
function [m1, m2] = GetMagDist(H_ex, M_s, H_a, n, a, J_ex, Js_x, Js_z, mu)

M_eff = M_s - H_a; %Calculate the effective field from demagnetizaztion and surface anisotropy

% construct the matrix
A_ii = (H_ex + M_s)*mu*a*M_s + 2*J_ex; %create the diagonal element (Inside layers are not affected by the surface anisotropy)
Diag = linspace(A_ii, A_ii, n); %create the vector for matrix diagonal 
subDiag_1 = linspace(-J_ex, -J_ex, n-1);%create the subdiagonal vectors
subDiag_2 = linspace(-J_ex, -J_ex, n-1); 

%Calculate the effective field caused by spin current
for i = 1:1:n 
    
    dJs_x(i) = Js_x(i)-Js_x(i+1); %Calculate the spin current loss in each layers (total n layers, compared to n+1 interfaces)
    dJs_z(i) = Js_z(i)-Js_z(i+1);

end 

M = diag(Diag) + diag(subDiag_1, 1) + diag(subDiag_2, -1); %construct the Matrix in euqation 5
M(1,1) = (H_ex + M_eff)*mu*a*M_s + J_ex;
M(n,n) = (H_ex + M_eff)*mu*a*M_s + J_ex;
% calculate magnetization distribution 
m1 = M\dJs_x.';
m2 = M\dJs_z.';










