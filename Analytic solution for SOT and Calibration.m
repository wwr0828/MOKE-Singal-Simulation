clear all

a = 4e-10;
tau = -2.4e-6;
c = 3e+8; %speed of light
A_ex = 0.64*2*19.1e-12; %exchange stiffness 
lambda = 780e-9; %wavelength
omega = 2*pi*c/lambda; %angular frequency
mu_Ms = 1.09; %(tesla)
mu = 4*pi*1e-7; %permeability (N/A^2)
%H_a = 79.77*27074/2/(a/1e-9); %surface anisotropy effective field (A/m)
M_s = mu_Ms/mu;
%M_eff = M_s - H_a;
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of (0.0036-1i*0.011)
n1 = 2.38+1i*4.36; % refractive index of the FM layer(for Py,2.2+1i*4.2)
n0 = 1.00 + 1i*0.00;
n2 = 1.43+1i*0.00; %refractive index of silica
delta_L = -c/(2*1i*n1*omega); %complex penetration depth
m_z = 0.85e-3;
lambda2 = sqrt(A_ex/(mu_Ms*M_s));
m_0 = tau/(mu_Ms*lambda2*M_s);



theta_k_cali = 1i*n0*n1*Q_FM*m_z/(n1^2-n0^2); %Calibration MOKE from front
theta_k_cali2 = -1i*n2*n1*Q_FM*m_z/(n1^2-n2^2); %Calibration MOKE from back
theta_k_SOT = -1i*n0*n1*Q_FM*m_0/(n1^2-n0^2) / (1+delta_L/lambda2);


