clear all

n2 = 2.38+1i*4.36; % refractive index of the FM layer(for Py,2.2+1i*4.2)
n1 = 1.43+1i*0.0; % n of FS
n0 = 1.00 + 1i*0.00;

R00 = (n0-n1)/(n0+n1);
P00 = abs(R00)^2;

R_out = (2*n0/(n0+n1))*((n1-n2)/(n1+n2))*(2*n1/(n0+n1));
P_out = abs(R_out)^2;

R_front = (n0-n2)/(n0+n2);
P_front = abs(R_front)^2;