%% Step 1, 2 & 3-- MOKE response as a fucntion of FM thickness

clear all
% Define all the parameters and constant 
%Parameters for step one
SHA_x = 0.16; %Spin Hall angle for x-spin spin current (conventional)
SHA_z = 0.01; %Spin Hall angle for z-spin spin current (spin rotation)
J_e = 1.25*0.4e+11; %charge current density  A/m^2)
%sigma = 9.43e+6; %electrical conductivity (S/m) Pt
%N_s = 1.6e+47; % electron density of states (/J/m^3)
l_sf = 1000e-9; %l_sf; spin-flip length, large! 1000nm (m)
l_dp = 1e-9; %l_dp; dephasing length, approx. 1 nm; (m)
l_ex = 0.3e-9; %l_ex; precession length, due to exchange coupling, approx 0.5 nm; (m)
a = 4e-10; %lattice constant--thickness of each sublayer (m) 0.4nm

%Parameters for step two
mu = 4*pi*1e-7; %permeability (N/A^2)
h_bar = 1.0546e-34; %planck constant (J*s)
e = 1.6e-19; % electron charge (C)
H_ex = 0; %external magnetic field (in-plane) (A/m)
M_s = 1.09/mu; %saturation magnetization (A/m) (from mu*Ms = 1.06T)
H_a = 79.77*22074/2/(a/1e-9); %surface anisotropy effective field (A/m) 
J_ex = 0.64*2*19.1e-12/a; %interlayer exchange strength (J/m^2)

%Parameters for step three
lambda = 405e-9;  % wavelength (m)
refrac_layer = 1.4+1i*2.6; % refractive index of the FM layer(for Py,2.2+1i*4.2)
refrac_Si = 5.49+1i*0.336; %refractive index of the substrate (Si) 3.7+1i*0.008
refrac_SiO2 = 1.47+1i*0.000; %refractive index of the oxidation (SiO2) 1.43+1i*0
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of Co (0.043+1i*0.007)
T_max = 80e-9; %The maximum thickness of the simulation (m)


i = 1;
Reflectance = zeros(round(T_max/a),1);
Phi_S = zeros(round(T_max/a),1);
Phi_P = zeros(round(T_max/a),1);
for T = a: a: T_max %FM layer total thickness in (m)

n = round(T/a); %iteration number (number of sublayers, not the number of interfaces)
% Calculate spin current and spin chemical potential distribution (First Step)

[Js_x, Js_z] = SpinCurrDist(SHA_x, SHA_z, J_e, T, l_sf, l_dp, l_ex, n, h_bar, e);

% Calculate SOT induced magnetization distribution (Second Step)

[m1, m2] = GetMagDist(H_ex, M_s, H_a, n, a, J_ex, Js_x, Js_z, mu); %m1, m2 represent the magnetization tilt caused by the two spin current, respectively

% sample layer (first is air, last is substrate, thickness of both do not matter)
refrac_temp = linspace(refrac_layer, refrac_layer, n); % Make a refractive index vector for all the sublayers
refrac = [1  refrac_temp  refrac_SiO2  refrac_Si];  % Construct the total refractive index vector 
Q_temp = linspace(Q_FM, Q_FM, n); % Make a Q vector for all the sublayers
Q = [0   Q_temp  0  0];  % Construct the total magneto-optical constants vector 
h_temp = linspace(a, a, n); % Make a thickness vector for all the sublayers
h = [inf   h_temp  1.0e-6  inf];  % Construct the overall thickness vector
mz = [0  m1'  0  0]; % Construct the overall magnetization vector in z-direction


[Reflectance(i), Phi_S(i), Phi_P(i)] = MOKE_MLmodel_SUB3(lambda, h, refrac, Q, mz); %Calculate the MOKE response for one thickness


i = i + 1; 

end 

%%
TT = a:a:T_max;
Thickness = 0.8*[5, 10, 15, 20, 30, 40, 60, 80]; %experiemental thickness points
%SOT = [17*1.62, 323*1.24, 856*1.11, 934*1.05, 1102*1.00, 1184*0.99, 1298*0.96, 1310*0.94]/5;
%err = [7*1.62, 10*1.24, 6*1.11, 23*1.05, 2*9*1.00, 2*8*0.99, 2*26*0.96, 2*12*0.94]/5; %error bar
SOT = [24, -59, -163, -316, -459, -549, -601, 0]/2;
%err = [7*1.62, 10*1.24, 6*1.11, 23*1.05, 2*9*1.00, 2*8*0.99, 2*26*0.96, 2*12*0.94]/5; %error bar

%data = [TT/1e-9; -real(Phi_S')*1e9];
%fileID = fopen('thickness dependence.txt','w');
%fprintf(fileID, '%f %e\n',data);
%fclose(fileID);

figure(2)
clf
plot(TT/1e-9,real(Phi_S)*1e9,'r','LineWidth', 4);
hold on
plot(TT/1e-9,imag(Phi_S)*1e9,'b','LineWidth', 4);
hold on
scatter(Thickness, SOT, 'k', 'linewidth', 4);
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
title('Singal FM layer MOKE signal at 405 nm wavelength')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])