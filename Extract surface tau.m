 clear all

B = 1e-9*[5-32*1i; (68-122*1i)/0.944]; %2*Measured Kerr angle from [Front;Back]
%B_E = 1e-9*[-134; -188]/2;

A = 1e-7*[0.69-0.81*1i -0.22-0.17*1i; -0.49-0.16*1i 1.26-2.20*1i]; %Simulated Front-back transformation matrix
%A_E = 1e-7*[-1.256 -0.2598; -0.2853 -1.9811];

tau = A\B; %Extracted surface torques
%tau_E = A_E\B_E;