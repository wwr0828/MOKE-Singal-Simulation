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
refrac_layer = 2.85+1i*4.06; % refractive index of the FM layer(for Py,2.2+1i*4.2)
%refrac_Cu = 0.2+1i*4.9; % refractive index of the Cu layer(for Py,2.2+1i*4.2)
refrac_Si = 3.7+1i*0.008; %refractive index of the substrate (Si) 
refrac_SiO2 = 1.47+1i*0.0; %refractive index of the oxidation (SiO2) 1.43+1i*0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%VERY SENSITIVE%%%%%%%%%%%%%%%%
%refrac_Al2O3 = 1.6716+1i*0; %refractive index of the oxidation (Al2O3) 1.6716+1i*0
%refrac_Pt = 2.73+1i*6.04; %refractive index of the (Pt)
%refrac_Cu = 0.153 + 1i*4.84; %refractive index of copper 0.153 + 1i*4.84
Q_FM = 0.012-1i*0.012; % magneto-optic coefficient of (0.0036-1i*0.011)
T_max = 60e-9; %The maximum thickness of the simulation (m)
m_z = 0.55e-3; %normalized magnetization tilt in z-direction, is actually m_z/M_s.(rad) (assume M_eff = 1.5T)

i = 1;
Reflectance = zeros(round(T_max/a),1);
Phi_S = zeros(round(T_max/a),1);
Phi_P = zeros(round(T_max/a),1);
for T = a: a: T_max %FM layer total thickness in (m)

% sample layer (first is air, last is substrate, thickness of both do not matter)
refrac = [1  refrac_layer  refrac_SiO2  refrac_Si];  % Construct the total refractive index vector 
Q = [0  Q_FM  0  0];  % Construct the total magneto-optical constants vector  (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)
h = [inf  T  1.0e-6  inf];  % Construct the overall thickness vector (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)
mz = [0  m_z  0  0]; % Construct the overall magnetization vector in z-direction (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)


[Reflectance(i), Phi_S(i), Phi_P(i)] = MOKE_MLmodel_SUB3(lambda, h, refrac, Q, mz); %Calculate the MOKE response for one thickness


i = i + 1; 

end 

%%
TT = a:a:T_max;
Thickness = 40; %experiemental thickness points

Calib = -1.4; % -5.33*0.94 (64Py) -18.12*1.05 (16Py)calibration Kerr rotation (urad)
Calib2 = -0.06; %calibration Kerr ellipticity (urad)
%err = [0.01*1.62, 0.05*1.24, 0.03*1.11, 0.03*1.05, 0.02*1.00, 0.01*0.99, 0.01*0.96, 0.01*0.94]/5; %error bar

%data = [TT/1e-9; real(Phi_S')*1e6];
%fileID = fopen('thickness calibration rotation.txt','w');
%fprintf(fileID, '%f %e\n',data);
%fclose(fileID);

%data2 = [TT/1e-9; imag(Phi_S')*1e6];
%fileID2 = fopen('thickness calibration ellipticity.txt','w');
%fprintf(fileID2, '%f %e\n',data2);
%fclose(fileID2);

figure(3)
clf
plot(TT/1e-9,real(Phi_S)*1e9,'r','LineWidth', 4);
hold on
plot(TT/1e-9,imag(Phi_S)*1e9,'b','LineWidth', 4);
hold on
scatter(Thickness, Calib*1000, 'k', 'linewidth', 4);
hold on
scatter(Thickness, Calib2*1000, 'b', 'linewidth', 4);
%hold on
%errorbar(Thickness, Calib*1000, err*1000, 'LineStyle','none');
hold off
%figure(4)
%plot(TT/1e-9, Reflectance, 'b','LineWidth', 4);
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
ylabel('Kerr rotation (nrad)')
xlabel('FM thickness (nm)')
title('Singal FM layer calibration MOKE signal at 780 nm wavelength (Ha = 10k)')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])

%% BACK----Calibration MOKE response as a fucntion of FM thickness

clear all


%Parameters for step two
%mu = 4*pi*1e-7; %permeability (N/A^2)
%H_ex = 1e+5; %external magnetic field (in-plane) (A/m)
%M_s = 1/mu; %saturation magnetization (A/m) (from mu*Ms = 1T)
%H_a = 0; %surface anisotropy effective field (A/m) assume to be 0


%Parameters for step three
a = 0.4e-10; %lattice constant (m)
lambda = 780e-9;  % wavelength (m)
refrac_layer = 2.38+1i*4.36; % refractive index of the FM layer(for Py,2.2+1i*4.2)
%refrac_Cu = 0.2+1i*4.9; % refractive index of the Cu layer(for Py,2.2+1i*4.2)
refrac_Si = 3.7+1i*0.008; %refractive index of the substrate (Si) 
refrac_SiO2 = 1.43+1i*0.0; %refractive index of the oxidation (SiO2) 1.43+1i*0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%VERY SENSITIVE%%%%%%%%%%%%%%%%
%refrac_Al2O3 = 1.6716+1i*0; %refractive index of the oxidation (Al2O3) 1.6716+1i*0
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of (0.0036-1i*0.011)
T_max = 80e-9; %The maximum thickness of the simulation (m)
m_z = 0.85e-3; %normalized magnetization tilt in z-direction, is actually m_z/M_s.(rad) (assume M_eff = 1T)

i = 1;
Reflectance = zeros(round(T_max/a),1);
Phi_S = zeros(round(T_max/a),1);
Phi_P = zeros(round(T_max/a),1);
for T = a: a: T_max %FM layer total thickness in (m)

% sample layer (first is air, last is substrate, thickness of both do not matter)
refrac = [refrac_SiO2  refrac_layer  1];  % Construct the total refractive index vector 
Q = [0  Q_FM  0];  % Construct the total magneto-optical constants vector  (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)
h = [inf  T  inf];  % Construct the overall thickness vector (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)
mz = [0  m_z  0]; % Construct the overall magnetization vector in z-direction (air/Sio2/Al2O3/Py/Al2O3/SiO2/Si)


[Reflectance(i), Phi_S(i), Phi_P(i)] = MOKE_MLmodel_SUB3(lambda, h, refrac, Q, mz); %Calculate the MOKE response for one thickness


i = i + 1; 

end 

%%
TT = a:a:T_max;
Thickness = 0.8*[5, 10, 15, 20, 30, 40, 60, 80, 100]; %experiemental thickness points

Calib = [2.85*1.62, -28.72*1.24, -22.88*1.11, -18.12*1.05, -14.07*1.00, -10.88*0.99, -7.01*0.96, -5.33*0.94, -6.85*0.94]/5; %calibration Kerr rotation (urad)
Calib2 = [-8.7*1.62, -5.72*1.24, -2.59*1.11, 0*1.05, -0.68*1.00, -0.24*0.99, 0.40*0.96, 0*0.94, 1.37*0.94]; %calibration Kerr ellipticity (urad)
%err = [0.01*1.62, 0.05*1.24, 0.03*1.11, 0.03*1.05, 0.02*1.00, 0.01*0.99, 0.01*0.96, 0.01*0.94]/5; %error bar

data = [TT/1e-9; real(Phi_S')*1e6];
fileID = fopen('thickness calibration rotation.txt','w');
fprintf(fileID, '%f %e\n',data);
fclose(fileID);

data2 = [TT/1e-9; imag(Phi_S')*1e6];
fileID2 = fopen('thickness calibration ellipticity.txt','w');
fprintf(fileID2, '%f %e\n',data2);
fclose(fileID2);

figure(3)
clf
plot(TT/1e-9,real(Phi_S)*1e9,'r','LineWidth', 4);
hold on
plot(TT/1e-9,imag(Phi_S)*1e9,'b','LineWidth', 4);
hold on
scatter(Thickness, Calib*1000, 'k', 'linewidth', 4);
hold on
scatter(Thickness, Calib2*1000, 'b', 'linewidth', 4);
%hold on
%errorbar(Thickness, Calib*1000, err*1000, 'LineStyle','none');
hold off
%figure(4)
%plot(TT/1e-9, Reflectance, 'b','LineWidth', 4);
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
ylabel('Kerr rotation (nrad)')
xlabel('FM thickness (nm)')
title('Singal FM layer calibration MOKE signal at 780 nm wavelength (Ha = 10k)')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])