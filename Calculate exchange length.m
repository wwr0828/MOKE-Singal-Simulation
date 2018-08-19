clear all

A_ex = 19.7e-12; %J/m
mu_0M_s = 2.24; %Tesla
mu = 4*pi*1e-7; %permeability (N/A^2)
Ms = mu_0M_s/mu;

lambda_ex = sqrt(2*A_ex/(mu_0M_s*Ms));% exchange length