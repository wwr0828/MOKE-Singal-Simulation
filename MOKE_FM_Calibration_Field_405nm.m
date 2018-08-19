%% Step 1, 2 & 3-- MOKE response as a fucntion of FM thickness

clear all


%Parameters for step two
%mu = 4*pi*1e-7; %permeability (N/A^2)
%H_ex = 1e+5; %external magnetic field (in-plane) (A/m)
%M_s = 1/mu; %saturation magnetization (A/m) (from mu*Ms = 1T)
%H_a = 0; %surface anisotropy effective field (A/m) assume to be 0


%Parameters for step three
a = 4e-10; %lattice constant (m)
lambda = 405e-9;  % wavelength (m)
refrac_layer = 1.4+1i*2.6; % refractive index of the FM layer(for Py,2.2+1i*4.2)
refrac_Si = 5.44+1i*0.34; %refractive index of the substrate (Si) 5.46+1i*0.97
refrac_SiO2 = 1.47+1i*0.000; %refractive index of the oxidation (SiO2) 1.43+1i*0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%VERY SENSITIVE%%%%%%%%%%%%%%%%
%refrac_Al2O3 = 1.6716+1i*0; %refractive index of the oxidation (Al2O3) 1.6716+1i*0
Q_FM = 0.0036-1i*0.011; % magneto-optic coefficient of Py(0.0036-1i*0.011)
T_max = 80e-9; %The maximum thickness of the simulation (m)
m_z = 0.85e-3; %normalized magnetization tilt in z-direction, is actually m_z/M_s.(rad) (assume M_eff = 1T)

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
Thickness = 0.8*[5, 10, 15, 20, 30, 40, 60, 80]; %experiemental thickness points

Calib = [2.40*1.62, 1.66*1.24, 2.61*1.11, 3.08*1.05, 3.03*1.00, 2.71*0.99, 2.54*0.96, 0*0.94]; %calibration Kerr rotation (urad)
%err = [0.01*1.62, 0.05*1.24, 0.03*1.11, 0.03*1.05, 0.02*1.00, 0.01*0.99, 0.01*0.96, 0.01*0.94]/5; %error bar

%data = [TT/1e-9; real(Phi_S')*1e9];
%fileID = fopen('thickness dependence_calibration.txt','w');
%fprintf(fileID, '%f %e\n',data);
%fclose(fileID);

figure(3)
clf
plot(TT/1e-9,-real(Phi_S)*1e9,'r','LineWidth', 4);
hold on
%plot(TT/1e-9,imag(Phi_S)*1e9,'b','LineWidth', 4);
%hold on
scatter(Thickness, Calib*1000, 'k', 'linewidth', 4);
hold on
%errorbar(Thickness, Calib*1000, err*1000, 'LineStyle','none');
%hold off
%figure(4)
%plot(TT/1e-9, Reflectance, 'b','LineWidth', 4);
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
ylabel('Kerr rotation (nrad)')
xlabel('FM thickness (nm)')
title('Singal FM layer calibration MOKE signal at 405 nm wavelength')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])