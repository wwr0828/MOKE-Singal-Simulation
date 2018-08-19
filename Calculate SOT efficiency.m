clear all

mu = 4*pi*1e-7; %permeability (N/A^2)
h_bar = 1.0546e-34; %planck constant (J*s)
e = 1.6e-19; % electron charge (C)
h_sot = 3.70*80; %A/m
j_e = 0.5*1e+11; %A/m^2
mu0_Ms = 1.09; %T
d = 5e-9; %m
beta = h_sot/j_e % m

sigma = 2*e*beta*mu0_Ms*d/h_bar;
tau = h_sot*mu0_Ms*d