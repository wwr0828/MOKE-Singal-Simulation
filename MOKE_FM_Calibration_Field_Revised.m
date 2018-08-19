%% Step 1, 2 & 3-- MOKE response as a fucntion of FM thickness

clear all


%Parameters for step two
%mu = 4*pi*1e-7; %permeability (N/A^2)
%H_ex = 1e+5; %external magnetic field (in-plane) (A/m)
%M_s = 1/mu; %saturation magnetization (A/m) (from mu*Ms = 1T)
%H_a = 0; %surface anisotropy effective field (A/m) assume to be 0


%Parameters for step three
a = 0.4e-10; %lattice constant (m)
lambda = 780e-9;  % wavelength (m)
refrac_layer = 2.2+1i*4.2; % refractive index of the FM layer(for Py,2.2+1i*4.2)
refrac_Si = 3.7+1i*0.007; %refractive index of the substrate (Si) 3.7+1i*0.007
refrac_SiO2 = 1.45+1i*0; %refractive index of the oxidation (SiO2) 1.47+1i*0
%refrac_Al2O3 = 1.6716+1i*0; %refractive index of the oxidation (Al2O3) 1.6716+1i*0
Q_FM = 0.036-1i*0.020; % magneto-optic coefficient of Co (0.043+1i*0.007)
T_max = 100e-9; %The maximum thickness of the simulation (m)
m_z = 0.95e-3; %normalized magnetization tilt in z-direction, is actually m_z/M_s.(rad)

i = 1;
Reflectance1 = zeros(round(T_max/a),1);
Phi_S1 = zeros(round(T_max/a),1);
Phi_P1 = zeros(round(T_max/a),1);
Reflectance2 = zeros(round(T_max/a),1);
Phi_S2 = zeros(round(T_max/a),1);
Phi_P2 = zeros(round(T_max/a),1);
Reflectance3 = zeros(round(T_max/a),1);
Phi_S3 = zeros(round(T_max/a),1);
Phi_P3 = zeros(round(T_max/a),1);
for T = a: a: T_max %FM layer total thickness in (m)

% sample layer (first is air, last is substrate, thickness of both do not matter)
refrac1 = [1  refrac_layer  refrac_SiO2  refrac_Si];  % Construct the total refractive index vector for Front(Si) 
refrac2 = [1  refrac_SiO2  refrac_layer  1];  % Construct the total refractive index vector for BACK refrac_SiO2
refrac3 = [1  refrac_layer  refrac_SiO2  1];  % Construct the total refractive index vector for FRONT 
Q1 = [0   Q_FM  0  0];  % Construct the total magneto-optical constants vector for FRONT(Si) 
Q2 = [0   0  Q_FM  0];  % Construct the total magneto-optical constants vector for BACK
Q3 = [0   Q_FM  0  0];  % Construct the total magneto-optical constants vector for FRONT
h1 = [inf  T  1e-6  1e-3];  % Construct the overall thickness vector for FRONT(Si)
h2 = [inf  0.5e-3  T  inf];  % Construct the overall thickness vector for BACK
h3 = [inf  T  0.5e-3  inf];  % Construct the overall thickness vector for FRONT
mz1 = [0  m_z 0 0]; % Construct the overall magnetization vector in z-direction for FRONT(Si)
mz2 = [0  0  -m_z  0]; % Construct the overall magnetization vector in z-direction for BACK
mz3 = [0  m_z  0  0]; % Construct the overall magnetization vector in z-direction for FRONT



[Reflectance1(i), Phi_S1(i), Phi_P1(i)] = MOKE_MLmodel_SUB3(lambda, h1, refrac1, Q1, mz1); %Calculate the MOKE response for one thickness for FRONT(Si)
[Reflectance2(i), Phi_S2(i), Phi_P2(i)] = MOKE_MLmodel_SUB3(lambda, h2, refrac2, Q2, mz2); %Calculate the MOKE response for one thickness for BACK
[Reflectance3(i), Phi_S3(i), Phi_P3(i)] = MOKE_MLmodel_SUB3(lambda, h3, refrac3, Q3, mz3); %Calculate the MOKE response for one thickness for FRONT


i = i + 1; 

end 

%%
TT = a:a:T_max;
Thickness1 = [5, 10, 15, 20, 30, 40, 60, 80]; %experiemental thickness points
Thickness2 = [40, 40, 40]; %experiemental thickness points
Calib1 = [2.85, -28.72, -22.88, -18.12, -14.07, -10.88, -7.01, -5.33]; %calibration Kerr rotation (urad)
Calib2 = [-10.88, 11.26, -9.46]; %calibration Kerr rotation (urad)

figure(3)
clf
plot(TT/1e-9,-real(Phi_S1)*1e9,'r','LineWidth', 4);
hold on
plot(TT/1e-9,-real(Phi_S2)*1e9,'b','LineWidth', 4);
hold on
plot(TT/1e-9,-real(Phi_S3)*1e9,'k','LineWidth', 4);
hold on
scatter(40, -10.88*1000, 'r', 'linewidth', 4);
hold on
scatter(40, 11.26*1000, 'b', 'linewidth', 4);
hold on
scatter(40, -9.46*1000, 'k', 'linewidth', 4);
hold off
%figure(4)
%plot(TT/1e-9, Reflectance, 'b','LineWidth', 4);
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
%axis([40,50, -2e4, 2e4]);
ylabel('Kerr rotation (nrad)')
xlabel('FM thickness (nm)')
title('Singal FM layer calibration MOKE signal at 780 nm wavelength')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])