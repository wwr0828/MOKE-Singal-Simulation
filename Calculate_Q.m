clear all

n0 = 1;
n1 = 2.7+1i*3.7;
mz = 0.19; %mrad
theta_k = -0.00063-1i*0.00081; %mrad

Q = theta_k * (n1^2 - n0^2)/(1i*n0*n1*mz)