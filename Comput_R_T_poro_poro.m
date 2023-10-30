

% Time 2019-7-8 by Haorui Peng in Utrecht Seismological group
% algorithm by Jun Yang,1992, paper:"IMPORTANCE OF FLOW CONDITION ON
%SEISMIC WAVES AT A SATURATED POROUS SOLID BOUNDARY"
% calculate relection and transmission coefficients of displacement at interface
%of elastic solid and porous solid saturated by viscous liquid
% Also "(1998_Jun Yang)Influence of Viscous Coupling
% on Seismic Reflection and Transmission in Saturated Porous Media"
%comparing with (Herbert Deresiewicz_1963)On Uniqueness in Dynamic Poroelasticity
% of equation (19), the normal stress equation should consider the fluid
% stress as well, so I changed the equation (26) in jun_Yang's paper
%Tzz(poro) +Pf = Tzz(elastic)
%(Wu 1990)Reflection and Transmission of Elastic Waves from a Fluid-Saturated Porous Solid Boundary


clear;
close all;


format long

ii = sqrt(-1);

freq = 15;
omega = 2*pi*freq;

%input parameters
%material property can be found in Table 5 of Spectral element simulation
%of wave propagation in porous media, Christina Morency et al., GJI, 2008.
%% top porous medium
poro1 = 0.2; %porosity
Kb1 = 2.2e9; %Bulk modulus of solid skeleton
Ks1 = 37e9;
Kf1 = 1.7e9;
miu_b1 = 4.4e9;
den_f1 = 750;
den_s1 = 2650;
den_p1 = (1-poro1)*den_s1 + poro1*den_f1;
k1 = 1e-10;  %permeability;
% visco1 = 1e-3;
visco1 = 0;
b1 = visco1/k1;
Kd1 = Ks1*(1+poro1*(Ks1/Kf1 - 1));
alpha1 = 1 - Kb1/Ks1;
M1 = Ks1^2/(Kd1 - Kb1);
lambda_b1 = Kb1 - 2/3*miu_b1;
lambda_c1 = lambda_b1+alpha1^2*M1;
c1 = 2; %structure factor
m1 = c1*den_f1/poro1; %mass coupling effect

N1 = miu_b1;
P1 = ((1 - poro1)*(1 - poro1 - Kb1/Ks1)*Ks1 + poro1*Ks1*Kb1/Kf1)/(1-poro1 - Kb1/Ks1 +poro1*Ks1/Kf1) + 4/3*N1;
Q1 = (1 - poro1 - Kb1/Ks1)*poro1*Ks1/(1 - poro1 - Kb1/Ks1 + poro1*Ks1/Kf1);
R1  = poro1^2*Ks1/(1 - poro1 - Kb1/Ks1 + poro1*Ks1/Kf1);

[Qp11,Qp12,Qs1,Vp11,Vp12,Vs1] = calculate_poroelastic_velocity_correction(freq,poro1,den_s1,den_f1,Ks1,Kf1,Kb1,miu_b1,c1,k1,visco1);


%% bottom porous medium
poro2 = 0.4; %porosity
Kb2 = 6.7e9; %Bulk modulus of solid skeleton
Ks2 = 6.9e9;
Kf2 = 2.0e9;
miu_b2 = 3.0e9;
den_f2 = 950;
den_s2 = 2200;
den_p2 = (1-poro2)*den_s2 + poro2*den_f2;
k2 = 1e-10;  %permeability;
% visco2 = 1e-3;
visco2 = 0;
b2 = visco2/k2;
Kd2 = Ks2*(1+poro2*(Ks2/Kf2 - 1));
alpha2 = 1 - Kb2/Ks2;
M2 = Ks2^2/(Kd2 - Kb2);
lambda_b2 = Kb2 - 2/3*miu_b2;
lambda_c2 = lambda_b2+alpha2^2*M2;
c2 = 2; %structure factor
m2 = c2*den_f2/poro2; %mass coupling effect

N2 = miu_b2;
P2 = ((1 - poro2)*(1 - poro2 - Kb2/Ks2)*Ks2 + poro2*Ks2*Kb2/Kf2)/(1-poro2 - Kb2/Ks2 +poro2*Ks2/Kf2) + 4/3*N2;
Q2 = (1 - poro2 - Kb2/Ks2)*poro2*Ks2/(1 - poro2 - Kb2/Ks2 + poro2*Ks2/Kf2);
R2  = poro2^2*Ks2/(1- poro2 - Kb2/Ks2 + poro2*Ks2/Kf2);

[Qp21,Qp22,Qs2,Vp21,Vp22,Vs2] = calculate_poroelastic_velocity_correction(freq,poro2,den_s2,den_f2,Ks2,Kf2,Kb2,miu_b2,c2,k2,visco2);


%critical angles
theta_critical_plus = asin(Vp11/Vp21);


%% get L1 L2 L3; equation 14,15, poro1
A1  = M1*(lambda_b1 + 2*miu_b1 + alpha1^2*M1) - (alpha1*M1)^2;
B1 = 2*den_f1*omega^2*alpha1*M1 - den_p1*omega^2*M1 - ...
    (lambda_b1 + 2*miu_b1 + alpha1^2*M1)*(m1*omega^2 - ii*b1*omega);
C1 = den_p1*omega^2*(m1*omega^2 - ii*b1*omega) - (den_f1*omega^2)^2;

L1_square1 = (-B1 - sqrt(B1^2 - 4*A1*C1))/(2*A1);
L2_square1 = (-B1 + sqrt(B1^2 - 4*A1*C1))/(2*A1);


delta1(1) = ((lambda_c1+2*miu_b1)*L1_square1 - den_p1*omega^2)/...
    (den_f1*omega^2 - alpha1*M1*L1_square1);
delta1(2) = ((lambda_c1+2*miu_b1)*L2_square1 - den_p1*omega^2)/...
    (den_f1*omega^2 - alpha1*M1*L2_square1);

G_plus1 = -(delta1(1)/poro1 + 1);
G_minus1 = -(delta1(2)/poro1 + 1);



den12_1 = - (c1-1)*poro1*den_f1;
den11_1 = (1-poro1)*den_s1 + (c1-1)*poro1*den_f1;
den22_1 = c1*den_f1/poro1;

Ls_square1 = (den_p1*omega^2 - (den_f1*omega^2)^2/(m1*omega^2 - ii*b1*omega))/miu_b1;
delta1(3) = (miu_b1*Ls_square1 - den_p1*omega^2)/(den_f1*omega^2);

H1 = delta1(3)/poro1 + 1;


%
%% get L1 L2 L3; equation 14,15, poro2
A2  = M2*(lambda_b2 + 2*miu_b2 + alpha2^2*M2) - (alpha2*M2)^2;
B2 = 2*den_f2*omega^2*alpha2*M2 - den_p2*omega^2*M2 - ...
    (lambda_b2 + 2*miu_b2 + alpha2^2*M2)*(m2*omega^2 - ii*b2*omega);
C2 = den_p2*omega^2*(m2*omega^2 - ii*b2*omega) - (den_f2*omega^2)^2;

L1_square2 = (-B2 - sqrt(B2^2 - 4*A2*C2))/(2*A2);
L2_square2 = (-B2 + sqrt(B2^2 - 4*A2*C2))/(2*A2);


delta2(1) = ((lambda_c2+2*miu_b2)*L1_square2 - den_p2*omega^2)/...
    (den_f2*omega^2 - alpha2*M2*L1_square2);
delta2(2) = ((lambda_c2+2*miu_b2)*L2_square2 - den_p2*omega^2)/...
    (den_f2*omega^2 - alpha2*M2*L2_square2);

G_plus2 = -(delta2(1)/poro2 + 1);
G_minus2 = -(delta2(2)/poro2 + 1);

den12_2 = - (c2-1)*poro2*den_f2;
den11_2 = (1-poro2)*den_s2 + (c2-1)*poro2*den_f2;
den22_2 = c2*den_f2/poro2;

Ls_square2 = (den_p2*omega^2 - (den_f2*omega^2)^2/(m2*omega^2 - ii*b2*omega))/miu_b2;
delta2(3) = (miu_b2*Ls_square2 - den_p2*omega^2)/(den_f2*omega^2);

H2 = delta2(3)/poro2 + 1;




for n = 1:90*1+1
    % for n = 1
    theta = (n-1)*pi/2/90;
    %wavenumber
    sigma = omega/Vp11*sin(theta); % Kx
    garma_plus1 = omega/Vp11*cos(theta); %Kz for incident fast P
    garma_minus1 = omega/Vp12*cos(theta); %Kz for reflected slow P
    garma_sh1 = sqrt((omega/Vs1)^2 - sigma^2); %Kz for reflected S
    
    theta_plus1 = asin(sigma/(omega/Vp11));
%     theta_minus1 = asin(sigma/(omega/Vp12));
    theta_sh1 = asin(sigma/(omega/Vs1));
    
    theta_plus2 = asin(sigma/(omega/Vp21));
    theta_minus2 = asin(sigma/(omega/Vp22));
    theta_sh2 = asin(sigma/(omega/Vs2));
    
    Lsz = sqrt(Ls_square2 - sigma^2);
    L2z = sqrt(L2_square2 - sigma^2);
    L1z = sqrt(L1_square2 - sigma^2);
    
    garma_plus2 = L1z;
    garma_minus2 = L2z;
    garma_sh2 = Lsz;
    
    X = zeros(6,6); %inpermeable interface
    X(1,1) = (P1 - 2*N1 + Q1)*(-sigma^2) + (P1+Q1)*(-garma_plus1^2) + (Q1 + R1)*G_plus1*(sigma^2 + garma_plus1^2); %B_plus1 reflection for fast P wave
    X(1,2) = (P1 - 2*N1 + Q1)*(-sigma^2) + (P1+Q1)*(-garma_minus1^2) + (Q1 + R1)*G_minus1*(sigma^2 + garma_minus1^2);   %B_mius1 reflection for slow P wave
    X(1,3) =  -2*N1*sigma*garma_sh1;   %C1 reflection for S wave
    X(1,4) = (P2 - 2*N2 + Q2)*(sigma^2) + (P2+Q2)*(garma_plus2^2) - (Q2 + R2)*G_plus2*(sigma^2 + garma_plus2^2); %B_plus2 transmission for fast P wave
    X(1,5) = (P2 - 2*N2 + Q2)*(sigma^2) + (P2+Q2)*(garma_minus2^2) - (Q2 + R2)*G_minus2*(sigma^2 + garma_minus2^2); %B_minus2 transmission for slow P wave
    X(1,6) =  -2*N2*sigma*garma_sh2;   %C2 transmission for S wave
    
    X(2,1) = 2*N1*garma_plus1*sigma;  %B_plus1 reflection for fast P wave
    X(2,2) = 2*N1*garma_minus1*sigma;  %B_mius1 reflection for slow P wave
    X(2,3) = N1*(sigma^2-garma_sh1^2); %C1 reflection for S wave
    X(2,4) = 2*N2*garma_plus2*sigma;   %B_plus2 transmission for fast P wave
    X(2,5) = 2*N2*garma_minus2*sigma;  %B_mius2 transmission for slow P wave
    X(2,6) = N2*(-sigma^2+garma_sh2^2); %C2 transmission for S wave
    
    X(3,1) = -1/poro1*Q1*(sigma^2+ garma_plus1^2) +1/poro1*R1*G_plus1*(sigma^2+ garma_plus1^2); %B_plus1 reflection for fast P wave
    X(3,2) = -1/poro1*Q1*(sigma^2+ garma_minus1^2) +1/poro1*R1*G_minus1*(sigma^2+ garma_minus1^2);  %B_mius1 reflection for slow P wave
    X(3,3) = 0;  %C1 reflection for S wave
    X(3,4) = 1/poro2*Q2*(sigma^2+ garma_plus2^2) -1/poro2*R2*G_plus2*(sigma^2+ garma_plus2^2);  %B_plus2 transmission for fast P wave
    X(3,5) = 1/poro2*Q2*(sigma^2+ garma_minus2^2) -1/poro2*R2*G_minus2*(sigma^2+ garma_minus2^2); %B_minus2 transmission for slow P wave
    X(3,6) = 0; %C2 transmission for S wave
    
    X(4,1) = garma_plus1;
    X(4,2) = garma_minus1;
    X(4,3) = sigma;
    X(4,4) = garma_plus2;
    X(4,5) = garma_minus2;
    X(4,6) = -sigma;
    
    X(5,1) = sigma;
    X(5,2) = sigma;
    X(5,3) = -garma_sh1;
    X(5,4) = -sigma;
    X(5,5) = -sigma;
    X(5,6) = -garma_sh2;
    
    X(6,1) = poro1*(G_plus1*garma_plus1 + garma_plus1);
    X(6,2) = poro1*(G_minus1*garma_minus1 + garma_minus1);
    X(6,3) = poro1*(-H1*sigma + sigma);
    X(6,4) = poro2*(G_plus2*garma_plus2 + garma_plus2);
    X(6,5) = poro2*(G_minus2*garma_minus2 + garma_minus2);
    X(6,6) = poro2*(H2*sigma - sigma);
    
    R_e =  zeros(6,1);
    R_e(1,1) = (P1 - 2*N1 + Q1)*sigma^2 + (P1 +Q1)*garma_plus1^2 - (Q1 + R1)*G_plus1*(sigma^2 + garma_plus1^2);
    R_e(2,1) = 2*N1*garma_plus1*sigma;
    R_e(3,1) = 1/poro1*Q1*(sigma^2 + garma_plus1^2) - 1/poro1*R1*G_plus1*(sigma^2 + garma_plus1^2);
    R_e(4,1) = garma_plus1;
    R_e(5,1) = -sigma;
    R_e(6,1) = poro1*(G_plus1*garma_plus1 + garma_plus1);
    
    if n-1 ==0
        Q_e = pinv(X)*R_e;
    else
        Q_e = X\R_e;
    end
    
    Ai = 1;
    B_plus1 = abs(Q_e(1,1));
    B_minus1 = abs(Q_e(2,1));
    C1 = abs(Q_e(3,1));
    

    B_plus2 = abs(Q_e(4,1));
    B_minus2 = abs(Q_e(5,1));
    
    C2 = abs(Q_e(6,1));
    
    
    %calculate Poynting vector Wu's paper equation (28)
    ux_in = ii*sigma;
    uz_in = -ii*garma_plus1;
    
    ux_plus1 = ii*B_plus1*sigma;
    ux_minus1 = ii*B_minus1*sigma;
    ux_sh1 = -ii*C1*garma_sh1;
    
    uz_plus1 = ii*B_plus1*garma_plus1;
    uz_minus1 = ii*B_minus1*garma_minus1;
    uz_sh1 = ii*sigma*C1;
    
    ux_plus2 = ii*B_plus2*sigma;
    ux_minus2 = ii*B_minus2*sigma;
    ux_sh2 = ii*garma_sh2*C2;
    
    uz_plus2= -ii*B_plus2*garma_plus2;
    uz_minus2 = -ii*B_minus2*garma_minus2;
    uz_sh2 = ii*sigma*C2;
 
    
    R_P1_disp_ratio_u(n,1) = abs(sqrt((uz_plus1^2+ux_plus1^2)/(uz_in^2+ux_in^2)));
    R_P2_disp_ratio_u(n,1) = abs(sqrt((uz_minus1^2+ux_minus1^2)/(uz_in^2+ux_in^2)));
    R_S_disp_ratio_u(n,1) = abs(sqrt((uz_sh1^2+ux_sh1^2)/(uz_in^2+ux_in^2)));
    T_P1_disp_ratio_u(n,1) = abs(sqrt((uz_plus2^2+ux_plus2^2)/(uz_in^2+ux_in^2)));
    T_P2_disp_ratio_u(n,1) = abs(sqrt((uz_minus2^2+ux_minus2^2)/(uz_in^2+ux_in^2)));
    T_S_disp_ratio_u(n,1) = abs(sqrt((uz_sh2^2+ux_sh2^2)/(uz_in^2+ux_in^2)));
end


R_T_disp_poro_poro = [R_P1_disp_ratio_u R_P2_disp_ratio_u R_S_disp_ratio_u T_P1_disp_ratio_u T_P2_disp_ratio_u T_S_disp_ratio_u];


line_width = 8
circle_size =30
font_size = 40;
lw =4
close all

figure(3);
plot(0:90,R_P1_disp_ratio_u(:,1),'r-','LineWidth',lw); hold on;
plot(0:90,R_P2_disp_ratio_u(:,1),'k-','LineWidth',lw); hold on;
plot(0:90,R_S_disp_ratio_u(:,1),'y-','LineWidth',lw); hold on;
plot(0:90,T_P1_disp_ratio_u(:,1),'g-','LineWidth',lw); hold on;
plot(0:90,T_P2_disp_ratio_u(:,1),'m-','LineWidth',lw); hold on;
plot(0:90,T_S_disp_ratio_u(:,1),'c-','LineWidth',lw); hold on;
xlabel('P wave Incidence angle');
ylabel('displacement ratio');
legend('Reflected fast P wave','Reflected slow P wave','Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave','FontSize',15);


% legend( 'Transmitted P1 wave',  'Transmitted S wave');
%         title('displacement ratio for open-pore acoustic-poro interface');

