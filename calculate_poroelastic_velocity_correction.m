function[Qp1,Qp2,Qs,Cp1,Cp2,Cs] = calculate_poroelastic_velocity_correction(freq,poro,den_s,den_f,Ks,Kf,Kb,miu_b,c,k,visco)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%by Haorui Peng 31-10-2018
%correcting Qp1 Qp2 calculation in equation (201) in Christina Morency's 2008 GJI
%paper, the code uses equation(198) and (199), parameters in Table 2 but
%the velocity Cp2 does not match case a, possibly because of the frequency-dependent
%term ufr(w) and b(w) are considered as constants, or the velocity in Table
%2 for case a is not right.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % clear
% freq = 30;
% 
% poro = 0.3; %porosity
% den_f = 880;
% den_s = 2650;
% 
% Ks = 12.2e9;
% Kf = 1.985e9;
% Kb = 9.6e9; %Bulk modulus of solid skeleton
% miu_b = 0;
% c = 1.2; %structure factor
% k = 1e-10;  %permeability;
% visco = 1e-3;
% % visco = 0;
b = visco/k;
omega = 2*pi*freq;
den_p = (1-poro)*den_s + poro*den_f;


m = c*den_f/poro; %mass coupling effect

D = Ks*(1+poro*(Ks/Kf - 1));
H = (Ks - Kb)^2/(D - Kb) + Kb + 4/3*miu_b;

C = Ks*(Ks - Kb)/(D - Kb);
M = Ks^2/(D - Kb);
B = H - 4/3*miu_b;
ii = sqrt(-1);
A1 = 1;
fc = visco*poro/2/pi/c/den_f/k

%
y0 = poro/c;
y1=(den_f*c*den_p-poro*den_f^2)/poro/den_p;
y2 = den_f/den_p;
y3=(den_p-y0*den_f)*omega^2;
y4=y0*C-H;
y5=C-y0*M;
y6=ii*y0*b*omega;
y7=ii*b*omega;
y8=y1*omega^2;
y9=C-y2*H;

x2 = y9*y5+ M*y4 - y2*C*y4;
x1 = -y8*y4-y9*y6+M*y3-y2*C*y3+y7*y4;
x0=-y8*y3+y7*y3;

L1 = (-x1+sqrt(x1^2-4*x2*x0))/(2*x2);
L2 = (-x1-sqrt(x1^2-4*x2*x0))/(2*x2);
Zp1 = omega^2/L1;
Zp2 = omega^2/L2;


Qp1= real(Zp1)/imag(Zp1);
Qp2= real(Zp2)/imag(Zp2);


Zp1 = sqrt(Zp1);
Zp2 = sqrt(Zp2);
Cp1 = 1/real(1/Zp1);
Cp2= 1/real(1/Zp2);


B2 =1;
y0=1i*poro/c*b*B2*omega;
y1=(den_p-poro*den_f/c)*omega^2;
y2=(den_f*c*den_p-poro*den_f^2)/(poro*den_p)*B2*omega^2;
y3 = -den_f/den_p*miu_b;
y4 = 1i*b*B2*omega;

Ls_square = (y2*y1-y4*y1)/(y2*miu_b-y3*y0-y4*miu_b);
Zs  = omega/sqrt(Ls_square);
Qs= real(Zs^2)/imag(Zs^2);
Cs = 1/real(1/Zs);

% L = omega^2/Zp1;

% A2 = ((den_p-poro*den_f/c)*omega^2-(B+4/3*miu_b)*L1+poro/c*C*L1)/(C*L1-poro/c*M*L1-ii*poro/c*b*omega)*A1
% 
% 
% left = (den_f*c*den_p-poro*den_f^2)/poro/den_p*A2*omega^2
% right = C*A1*L1+M*A2*L1-den_f/den_p*((B+4/3*miu_b)*A1*L1+C*A2*L1)+ii*b*A2*omega


% B1 = y0/(miu_b*Ls_square - y1)
% left  = y2
% right=y3*Ls_square*B1+y4
