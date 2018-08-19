% This program allows you to calculate polar complex Kerr angle of singal FM layer with inuniform magnetiation distribution at normal incidence
% The theory can be found in 1998_You_JAP 

%% USER INPUT
lambda = 780e-9;  % wavelength (m)

% sample layer (first is air, last is substrate, thickness of both do not matter)
n = [1   2.1605+1i*0    1.6543+1i*3.1335  1.3810+1i*3.2074  1.4892+1i*0 5.4376+1i*0.3421];  % refractive indices (for Co,3.542+1i*4.6)
Q = [0   0   0   0.0063-1i*0.062  0   0];  % magneto-optical constants (again, not accurate value) 1e-6-1i*1e-7 *0.062 0.0063

%% SOLVE MULTILAYER MODEL
% example below calculates polar complex Kerr angle as function of angle of incidence

Vec_h = (0 : 0.1 : T);  %change FM thickness

Reflectance = zeros(length(Vec_h),1);
Phi_S = zeros(length(Vec_h),1);
Phi_P = zeros(length(Vec_h),1);
for i = 1:n
    h = [0   Vec_h(i)  0]*1e-9;  % layer thicknesses (nm)
    [Reflectance(i), Phi_S(i), Phi_P(i)] = MOKE_MLmodel_SUB3(lambda,h,n,Q,m1);
end

%%
figure(325)
clf
plot(Vec_h,real(Phi_S)*1e9,'k')
hold on
plot(Vec_h,real(Phi_P)*1e9,'k--')
figure(326)
plot(Vec_h, Reflectance, 'r')
%plot(Vec_h,Reflectance,'r')
%plot(Vec_h,imag(Phi_S)*1e9,'r')
%plot(Vec_h,imag(Phi_P)*1e9,'r--')  
ylabel('Kerr angle (nrad)')
xlabel('HfO2 thickness (nm)')
title('5Pt/3Py with Co/Pt Q')
%text(40,8,'Kerr rotation')
%text(5,-0.028,'Kerr ellipticity','Color',[1 0 0])
%text(45,22,'S')
%text(50,15,'P')
%text(50,-0.02,'S','Color',[1 0 0])
%text(50,-0.045,'P','Color',[1 0 0])