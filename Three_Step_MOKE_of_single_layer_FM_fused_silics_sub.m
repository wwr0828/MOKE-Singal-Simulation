%% Three steps of getting SOT induced MOKE response of a singal layer FM  (SI units)

%% Step 1 and 2: Calculate the spin current and magnetization distribution for a certain thickness
clear all
% Define all the parameters and constant 
%Parameters for step one
tau_T = -0*1.6e-6; %Spin Hall angle for x-spin spin current (conventional)
tau_B = 1.6e-6; %Spin Hall angle for z-spin spin current (spin rotation)
%J_e = 1.25*0.8e+11; %charge current density  A/m^2)
%sigma = 9.43e+6; %electrical conductivity (S/m) Pt
%N_s = 1.6e+47; % electron density of states (/J/m^3)
T = 32e-9; %thickness of FM (m) 80nm
a = 0.4e-10; %lattice constant (m) 4A
n = round(T/a); %iteration number (number of layers, not the number of interfaces)
%l_sf = 1000e-9; %l_sf; spin-flip length, large! 1000nm (m)
%l_dp = 1e-9; %l_dp; dephasing length, approx. 1 nm; (m)
%l_ex = 0.4e-9; %l_ex; precession length, due to exchange coupling, approx 0.5 nm; (m)

%Parameters for step two
mu = 4*pi*1e-7; %permeability (N/A^2)
%h_bar = 1.0546e-34; %planck constant (J*s)
%e = 1.6e-19; % electron charge (C)
H_ex = 0; %external magnetic field (in-plane) (A/m)
M_s = 1.09/mu; %saturation magnetization (A/m)
H_a = 0;%79.77*22074/2/(a/1e-9); %surface anisotropy effective field (A/m)
J_ex = 0.64*2*19.1e-12/a; %interlayer exchange strength (J/m^2)

% Calculate SOT induced magnetization distribution 

m1 = GetMagDist_2(H_ex, M_s, H_a, n, a, J_ex, tau_T, tau_B, mu); %m1 represent the magnetization tilt caused by SOT in first 1 layer
m2 = GetMagDist_3(H_ex, M_s, H_a, n, a, J_ex, tau_T, tau_B, mu); %m2 represent the magnetization tilt caused by SOT in first 5 layers
%%
tau = -1.5e-6;
c = 3e+8; %speed of light
A_ex = 0.64*2*19.1e-12; %exchange stiffness 
lambda = 780e-9; %wavelength
omega = 2*pi*c/lambda; %angular frequency
mu_Ms = 1.06; %(tesla)
mu = 4*pi*1e-7; %permeability (N/A^2)
%H_a = 79.77*27074/2/(a/1e-9); %surface anisotropy effective field (A/m)
M_s = mu_Ms/mu;
%M_eff = M_s - H_a;
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of (0.0036-1i*0.011)
n1 = 2.38+1i*4.36; % refractive index of the FM layer(for Py,2.2+1i*4.2)
n0 = 1.00 + 1i*0.00;
delta_L = c/(2*1i*n1*omega); %complex penetration depth
m_z = 0.85e-3;
lambda2 = sqrt(A_ex/(mu_Ms*M_s));
m_0 = tau/(mu_Ms*lambda2*M_s);
i = 1;

for z = a:a:32e-9
    m(i) = m_0*exp(-z/lambda2);
    
    i = i + 1;
end 

%% Plot the spin current distribution
t = 0:T/n:T;
tt = a:a:T;

%data = [t/1e-9; Js_x];
%fileID1 = fopen('Spin current distribution.txt','w');
%fprintf(fileID1, '%f %e\n',data);
%fclose(fileID1);

%data2 = [tt/1e-9; m1'/max(m1)];
%fileID2 = fopen('Magnetization distribution.txt','w');
%fprintf(fileID2, '%f %e\n',data2);
%fclose(fileID2);

figure(1)
plot(tt/1e-9,m1/max(m1),'r-','LineWidth', 4); %spin current distribution (two curves)
%set(gca,'FontSize',20)
xlabel('Thickness (nm)', 'FontSize', 20);
ylabel('m_z/m_(max) (rad)','FontSize', 20);
title('Magnetization distribution in 32 nm Py');
%axis([0 32 -1.2 1.2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1, 2 & 3-- MOKE response as a fucntion of FM thickness

clear all
% Define all the parameters and constant 
%Parameters for step one
tau_T = -1e-6;%-1.9e-6; %torque at top surface (conventional) (J/(T*m^2))
tau_B = 0;%1e-6;%-tau_T; %torque at bottom surface 
J_e = 1.25*0.2e+11; %charge current density  A/m^2)
%sigma = 9.43e+6; %electrical conductivity (S/m) Pt
%N_s = 1.6e+47; % electron density of states (/J/m^3)
%l_sf = 1000e-9; %l_sf; spin-flip length, large! 1000nm (m)
%l_dp = 1e-9; %l_dp; dephasing length, approx. 1 nm; (m)
l_ex = 4e-9; %l_ex; exchange length, due to exchange coupling, approx 4 nm; (m)
a = 4e-10; %lattice constant--thickness of each sublayer (m) 0.4nm

%Parameters for step two
mu = 4*pi*1e-7; %permeability (N/A^2)
h_bar = 1.0546e-34; %planck constant (J*s)
e = 1.6e-19; % electron charge (C)
H_ex = 0; %external magnetic field (in-plane) (A/m)
M_s = 1.09/mu; %saturation magnetization (A/m) (from mu*Ms = 1.06T)
H_a = 79.77*27074/2/(a/1e-9); %surface anisotropy effective field (A/m) 
J_ex = 0.64*2*19.1e-12/a; %interlayer exchange strength (J/m^2) %0.64 is from the scaling of thickness factor 0.8

%Parameters for step three
lambda = 780e-9;  % wavelength (m)
refrac_layer = 2.38+1i*4.36; %refractive index of the FM layer(for Py,2.2+1i*4.2) (14.78nm penetration depth)
refrac_Si = 3.7+1i*0.008; %refractive index of the substrate (Si) 3.7+1i*0.007
refrac_SiO2 = 1.47+1i*0.0; %refractive index of the oxidation (SiO2) 1.47+1i*0
refrac_Pt = 2.73+1i*6.04; %refractive index of the (Pt)
refrac_Cu = 0.153 + 1i*4.84; %refractive index of copper 0.153 + 1i*4.84
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of Co (0.043+1i*0.007)
T_max = 40e-9; %The maximum thickness of the simulation (m)


i = 1;
Reflectance = zeros(round(T_max/a),1);
Phi_S = zeros(round(T_max/a),1);
Phi_P = zeros(round(T_max/a),1);
for T = a: a: T_max %FM layer total thickness in (m)

n = round(T/a); %iteration number (number of sublayers, not the number of interfaces)

% Calculate SOT induced magnetization distribution (Second Step)

m1 = GetMagDist_2(H_ex, M_s, H_a, n, a, J_ex, tau_T, tau_B, mu); %m1, m2 represent the magnetization tilt caused by the two spin current, respectively

% sample layer (first is air, last is substrate, thickness of both do not matter)
refrac_temp = linspace(refrac_layer, refrac_layer, n); % Make a refractive index vector for all the sublayers
refrac = [1  refrac_Pt  refrac_Cu  refrac_temp  refrac_SiO2];  % Construct the total refractive index vector 
Q_temp = linspace(Q_FM, Q_FM, n); % Make a Q vector for all the sublayers
Q = [0  0  0  Q_temp  0];  % Construct the total magneto-optical constants vector 
h_temp = linspace(a, a, n); % Make a thickness vector for all the sublayers
h = [inf  3e-9  1.5e-9  h_temp  inf];  % Construct the overall thickness vector
mz = [0  0  0  m1'  0]; % Construct the overall magnetization vector in z-direction


[Reflectance(i), Phi_S(i), Phi_P(i)] = MOKE_MLmodel_SUB3(lambda, h, refrac, Q, mz); %Calculate the MOKE response for one thickness


i = i + 1; 

end 

SHA = 2*e*tau_T/(h_bar*J_e);

%%
TT = a:a:T_max;
%Thickness = 0.8*[5, 10, 15, 20, 30, 40, 60, 80]; %experiemental thickness points
%SOT = [17*1.62, 323*1.24, 856*1.11, 934*1.05, 1102*1.00, 1184*0.99, 1298*0.96, 1310*0.94]/5;
%err = [7*1.62, 10*1.24, 6*1.11, 23*1.05, 2*9*1.00, 2*8*0.99, 2*26*0.96, 2*12*0.94]/5; %error bar
%SOT = [17, 323, 856, 934, 1102, 1184, 1298, 1310]/10;%Kerr rotation
%SOT2 = [-113, -55, -112, 0, -249 ,-299 ,-308, 0]; 
%err = [7*1.62, 10*1.24, 6*1.11, 23*1.05, 2*9*1.00, 2*8*0.99, 2*26*0.96, 2*12*0.94]/5; %error bar

%data = [TT/1e-9; real(Phi_S')];
%fileID2 = fopen('tau_B_MOKE.txt','w');
%fprintf(fileID2, '%f %e\n',data);
%fclose(fileID2);

%data2 = [TT/1e-9; imag(Phi_S')*1e9];
%fileID2 = fopen('thickness dependence ellipticity.txt','w');
%fprintf(fileID2, '%f %e\n',data2);
%fclose(fileID2);

figure(2)
clf
plot(TT/1e-9,real(Phi_S)*1e9,'r','LineWidth', 4);
hold on
plot(TT/1e-9,imag(Phi_S)*1e9,'b','LineWidth', 4);
hold on
%scatter(Thickness, SOT, 'k', 'linewidth', 4);
%hold on
%scatter(Thickness, SOT2, 'b', 'linewidth', 4);
%hold on
%errorbar(Thickness, SOT, err, 'LineStyle','none');
hold off
%figure(4)
%plot(TT, Reflectance, 'r')
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
ylabel('Kerr angle (nrad)')
xlabel('FM thickness (nm)')
title('Singal FM layer MOKE signal at 780 nm wavelength')
%axis([10 80 -500 500]);
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])