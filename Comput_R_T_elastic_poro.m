

% Time 2017-12-5 by Haorui Peng in Utrecht Seismological group
% algorithm by Jun Yang,1992, paper:"IMPORTANCE OF FLOW CONDITION ON
%SEISMIC WAVES AT A SATURATED POROUS SOLID BOUNDARY"
% calculate amplitude relection and transmission coefficient between interface
%of elastic solid and porous solid saturated by viscous liquid
% Also "(1998_Jun Yang)Influence of Viscous Coupling
% on Seismic Reflection and Transmission in Saturated Porous Media"
%comparing with (Herbert Deresiewicz_1963)On Uniqueness in Dynamic Poroelasticity
% of equation (19), the normal stress equation should consider the fluid
% stress as well, so I changed the equation (26) in jun_Yang's paper
%Tzz(poro) +Pf = Tzz(elastic)
%(Wu 1990)Reflection and Transmission of Elastic Waves from a Fluid-Saturated Porous Solid Boundary


% plot parameter_type: 1 = displacement ratio; 2 = energy ratio; 3 =
% potential ratio
parameter_type = 1;
% plot parameter_type: 1 = z component; 2 = x component; 3 = norm of
% displacement x and z
% displacement vecter
parameter_wave = 3;

%input parameters
%material property can be found in Table 5 of Spectral element simulation
%of wave propagation in porous media, Christina Morency et al., GJI, 2008.
ii = sqrt(-1);
%input parameters
freq = 15;
omega = 2*pi*freq;

%top elastic medium
Vp = 2219;
Vs = 1325;
den_e = 2600;
miu_e = den_e*Vs^2;
lambda_e = den_e*Vp^2 - 2*miu_e;
poisson_ratio = (Vp^2-2*Vs^2)/(2*(Vp^2-Vs^2));

%% bottom porous medium
poro = 0.4; %porosity
Kb = 6.7e9; %Bulk modulus of solid skeleton
Ks = 6.9e9;
Kf = 2.0e9;
miu_b = 3.0e9;
den_f = 950;
den_s = 2200;
den_p = (1-poro)*den_s + poro*den_f;
k = 1e-10;  %permeability;
visco = 0;
b = visco/k;

Kd = Ks*(1+poro*(Ks/Kf - 1));
alpha = 1 - Kb/Ks;
M = Ks^2/(Kd - Kb);
lambda_b = Kb - 2/3*miu_b;
lambda_c = lambda_b+alpha^2*M;
c = 2; %structure factor
m = c*den_f/poro; %mass coupling effect

[Qp1,Qp2,Qs,Vp1,Vp2,Vs_p] = calculate_poroelastic_velocity_correction(freq,poro,den_s,den_f,Ks,Kf,Kb,miu_b,c,k,visco);

%critical angles
theta_critical_plus = asin(Vp/Vp1);
theta_critical_plus_d = theta_critical_plus/pi*180;

%%
N = miu_b;
P = ((1 - poro)*(1 - poro - Kb/Ks)*Ks + poro*Ks*Kb/Kf)/(1-poro - Kb/Ks +poro*Ks/Kf) + 4/3*N;
Q = (1 - poro - Kb/Ks)*poro*Ks/(1 - poro - Kb/Ks + poro*Ks/Kf);
R  = poro^2*Ks/(1 - poro - Kb/Ks + poro*Ks/Kf);


P_e = lambda_e + 2*miu_e;
N_e = miu_e;

% get L1 L2 L3; equation 14,15

A  = M*(lambda_b + 2*miu_b + alpha^2*M) - (alpha*M)^2;
B = 2*den_f*omega^2*alpha*M - den_p*omega^2*M - ...
    (lambda_b + 2*miu_b + alpha^2*M)*(m*omega^2 - ii*b*omega);
C = den_p*omega^2*(m*omega^2 - ii*b*omega) - (den_f*omega^2)^2;

L1_square = (-B - sqrt(B^2 - 4*A*C))/(2*A);
L2_square = (-B + sqrt(B^2 - 4*A*C))/(2*A);


delta(1) = ((lambda_c+2*miu_b)*L1_square - den_p*omega^2)/...
    (den_f*omega^2 - alpha*M*L1_square);
delta(2) = ((lambda_c+2*miu_b)*L2_square - den_p*omega^2)/...
    (den_f*omega^2 - alpha*M*L2_square);

G_plus = -(delta(1)/poro + 1);
G_minus = -(delta(2)/poro + 1);

den12 = - (c-1)*poro*den_f;
den11 = (1-poro)*den_s + (c-1)*poro*den_f;
den22 = c*den_f/poro;

G_plus_2 = (Vp1^2*den11 - P)/(Vp1^2*den12 - Q);
G_minus_2 = (Vp2^2*den11 - P)/(Vp2^2*den12 - Q);

G_plus - G_plus_2

Ls_square = (den_p*omega^2 - (den_f*omega^2)^2/(m*omega^2 - ii*b*omega))/miu_b;
delta(3) = (miu_b*Ls_square - den_p*omega^2)/(den_f*omega^2);

for n = 1:90*1+1
% for n = 58
    theta = (n-1)*pi/2/90;
    %wavenumber 
    sigma = omega/Vp*sin(theta);
    garma_f = omega/Vp*cos(theta);
    garma_S = sqrt((omega/Vs)^2 - sigma^2);
    sigma = omega/Vp*sin(theta);
    theta_S = asin(sigma/(omega/Vs));
    theta_plus = asin(sigma/(omega/Vp1));
    theta_minus = asin(sigma/(omega/Vp2));
    theta_sh = asin(sigma/(omega/Vs_p));
    
    Lsz = sqrt(Ls_square - sigma^2);
    L2z = sqrt(L2_square - sigma^2);
    L1z = sqrt(L1_square - sigma^2);

    garma_plus = L1z;
    garma_minus = L2z;
    garma_sh = Lsz;
    
    X = zeros(5,5); %inpermeable interface
    X(1,1) = (P_e - 2*N_e)*sigma^2 + P_e*garma_f^2; %Ar reflection for P wave
    X(1,2) = (2*N_e)*garma_S*sigma; %D reflection for S wave
    X(1,3) = -(P - 2*N + Q)*sigma^2 - (P + Q)*garma_plus^2 + (Q+R)*G_plus*omega^2/Vp1^2; %B_plus
    X(1,4) = -(P - 2*N + Q)*sigma^2 - (P + Q)*garma_minus^2 + (Q+R)*G_minus*omega^2/Vp2^2; %B_minus
    X(1,5) = 2*N*garma_sh*sigma; %D transmission for S
    X(2,1) = 2*N_e*garma_f*sigma; 
    X(2,2) = N_e*(sigma^2-garma_S^2);
    X(2,3) = 2*N*sigma*garma_plus;
    X(2,4) = 2*N*sigma*garma_minus;
    X(2,5) = N*(garma_sh^2 - sigma^2);
    X(3,1) = garma_f;
    X(3,2) = sigma;
    X(3,3) = garma_plus;
    X(3,4) = garma_minus;
    X(3,5) = -sigma;
    X(4,1) = sigma;
    X(4,2) = -garma_S;
    X(4,3) = -sigma;
    X(4,4) = -sigma;
    X(4,5) = -garma_sh;
    X(5,1) = 0;
    X(5,2) = 0;
    X(5,3) = (1+G_plus)*garma_plus;
    X(5,4) = (1+G_minus)*garma_minus;
    X(5,5) = sigma*((delta(3)/poro + 1) - 1);
    
    R_e =  zeros(5,1);
    R_e(1,1) = -(P_e - 2*N_e)*sigma^2 - P_e*garma_f^2;
    R_e(2,1) = 2*N_e*garma_f*sigma;
    R_e(3,1) = garma_f;
    R_e(4,1) = -sigma;
    R_e(5,1) = 0;
    if n-1 ==0
        Q_e = pinv(X)*R_e;
    else
        Q_e = X\R_e;
    end
    
     R_e =  zeros(5,1);
    R_e(1,1) = -(P_e - 2*N_e)*sigma^2 - P_e*garma_f^2;
    R_e(2,1) = 2*N_e*garma_f*sigma;
    R_e(3,1) = garma_f;
    R_e(4,1) = -sigma;
    R_e(5,1) = 0;
    if n-1 ==0
        Q_e = pinv(X)*R_e;
    else
        Q_e = X\R_e;
    end
   
    Ai = 1;
    P_ratio(n,:) = abs(Q_e(:));
    Ar = P_ratio(n,1);
    A_R(n,1) = P_ratio(n,1);
    As = P_ratio(n,2);
    A_S(n,1) = P_ratio(n,2);

    B_plus = P_ratio(n,3);
    B_PLUS(n,1) = P_ratio(n,3);

    
    B_minus = P_ratio(n,4);
    B_MINUS(n,1) = P_ratio(n,4);
    Bs = P_ratio(n,5);
    B_S(n,1) = P_ratio(n,5);
    
    %calculate Poynting vector Wu's paper equation (28)
    uz_r = -ii*Ar*garma_f;
    uz_in = ii*garma_f;
    ux_r = ii*Ar*sigma;
    ux_in = ii*sigma;
    I_in = omega^4*den_e/Vp*Ai^2*cos(theta);
    
    ux_S = 1i*garma_S*As;
    uz_S = 1i*sigma*As;
    
    T1_S = -2*N_e*sigma*garma_S*As;
    T3_S = 2*N_e*sigma*garma_S*As;
    T5_S = N_e*(garma_S^2 - sigma^2)*As;
    
    I_S_x = -1i*omega*(ux_S*T1_S + uz_S*T5_S);
    I_S_z = -1i*omega*(ux_S*T5_S + uz_S*T3_S);
    I_S = sqrt(I_S_x^2 + I_S_z^2)*cos(theta_S);
    
    ux_plus = ii*sigma*B_plus;
    uz_plus = -ii*garma_plus*B_plus;
    Ux_plus = -ii*G_plus*sigma*B_plus;
    Uz_plus = ii*G_plus*garma_plus*B_plus;
    T1_plus = ((Q*G_plus - P + 2*N)*omega^2/Vp1^2 - 2*N*sigma^2)*B_plus;
    T3_plus = ((Q*G_plus - P + 2*N)*omega^2/Vp1^2 - 2*N*garma_plus^2)*B_plus;
    T5_plus = 2*N*sigma*garma_plus*B_plus;
    S_plus = (R*G_plus - Q)*omega^2/Vp1^2*B_plus;
    
    I_plus_c(n,1) = -ii*omega*(ux_plus*T1_plus + uz_plus*T5_plus + Ux_plus*S_plus);
    I_plus_c(n,2) = -ii*omega*(ux_plus*T5_plus + uz_plus*T3_plus + Uz_plus*S_plus);
    I_plus = norm(I_plus_c(n,:))*cos(theta_plus);
    %     I_plus =  ((2*Q - R*G_plus)*G_plus - P)*omega^4/(Vp1^3)*B_plus^2*cos(theta_plus);
    
    
    ux_minus = ii*sigma*B_minus;
    uz_minus = -ii*garma_minus*B_minus;
    Ux_minus = -ii*G_minus*sigma*B_minus;
    Uz_minus = ii*G_minus*garma_minus*B_minus;
    T1_minus = ((Q*G_minus - P + 2*N)*omega^2/Vp2^2 - 2*N*sigma^2)*B_minus;
    T3_minus = ((Q*G_minus - P + 2*N)*omega^2/Vp2^2 - 2*N*garma_minus^2)*B_minus;
    T5_minus = 2*N*sigma*garma_minus*B_minus;
    S_minus = (R*G_minus - Q)*omega^2/Vp2^2*B_minus;
    
    I_minus_c(n,1) = -ii*omega*(ux_minus*T1_minus + uz_minus*T5_minus +Ux_minus*S_minus);
    I_minus_c(n,2) = -ii*omega*(ux_minus*T5_minus + uz_minus*T3_minus +Uz_minus*S_minus);
    I_minus = norm(I_minus_c(n,:))*cos(theta_minus);
    
    
    
    ux_sh = 1i*garma_sh*Bs;
    uz_sh = 1i*sigma*Bs;
    Ux_sh = 1i*garma_sh*(delta(3)/poro + 1)*Bs;
    Uz_sh = 1i*sigma*(delta(3)/poro + 1)*Bs;
    T1_sh = -2*N*sigma*garma_sh*Bs;
    T3_sh = 2*N*sigma*garma_sh*Bs;
    T5_sh = N*(garma_sh^2 - sigma^2)*Bs;
    S_sh = 0;
    
    I_sh_x = -1i*omega*(ux_sh*T1_sh + uz_sh*T5_sh + Ux_sh*S_sh);
    I_sh_z = -1i*omega*(ux_sh*T5_sh + uz_sh*T3_sh + Uz_sh*S_sh);
    I_sh = sqrt(I_sh_x^2 + I_sh_z^2)*cos(theta_sh);
    %     I_sh = N*omega^4/Vs_p^3*Bs^2*cos(theta_sh);
    
    %compute energy ratio
    %Wu's paper B+ is B_plus for fast P wave; B- is B_minus for slow P wave
    R_P_energy_ratio(n,1) = abs(Ar)^2; %reflected P
    R_S_energy_ratio(n,1) =  abs(I_S/I_in); % S wave
    T_P1_energy_ratio(n,1) =  abs(I_plus/I_in); % P1 wave
    T_P2_energy_ratio(n,1) =  abs(I_minus/I_in); % P2 wave
    T_S_energy_ratio(n,1) =  abs(I_sh/I_in); % S wave
    Energy_total_ratio(n,1) = R_P_energy_ratio(n,1) + R_S_energy_ratio(n,1) + T_P1_energy_ratio(n,1) + T_P2_energy_ratio(n,1) + T_S_energy_ratio(n,1);
    
    R_P_disp_ratio_uz(n,1) = abs(uz_r/uz_in);
    R_S_disp_ratio_uz(n,1) = abs(uz_S/uz_in);
    T_P1_disp_ratio_uz(n,1) = abs(uz_plus/uz_in);
    T_P2_disp_ratio_uz(n,1) = abs(uz_minus/uz_in);
    T_S_disp_ratio_uz(n,1) = abs(uz_sh/uz_in);
    
    R_P_disp_ratio_ux(n,1) = abs(ux_r/ux_in);
    R_S_disp_ratio_ux(n,1) = abs(ux_S/ux_in);
    T_P1_disp_ratio_ux(n,1) = abs(ux_plus/ux_in);
    T_P2_disp_ratio_ux(n,1) = abs(ux_minus/ux_in);
    T_S_disp_ratio_ux(n,1) = abs(ux_sh/ux_in);
    
    R_P_disp_ratio_u(n,1) = abs(sqrt((uz_r^2+ux_r^2)/(uz_in^2+ux_in^2)));
    R_S_disp_ratio_u(n,1) = abs(sqrt((uz_S^2+ux_S^2)/(uz_in^2+ux_in^2)));
    T_P1_disp_ratio_u(n,1) = abs(sqrt((uz_plus^2+ux_plus^2)/(uz_in^2+ux_in^2)));
    T_P2_disp_ratio_u(n,1) = abs(sqrt((uz_minus^2+ux_minus^2)/(uz_in^2+ux_in^2)));
    T_S_disp_ratio_u(n,1) = abs(sqrt((uz_sh^2+ux_sh^2)/(uz_in^2+ux_in^2)));
end


lw =4;

% image
if  parameter_type == 1
    if parameter_wave == 1
        figure(1);
        plot(0:90,R_P_disp_ratio_uz(:,1),'r-','LineWidth',lw); hold on;
        plot(0:90,R_S_disp_ratio_uz(:,1),'y-','LineWidth',lw); hold on;
        plot(0:90,T_P1_disp_ratio_uz(:,1),'g-','LineWidth',lw); hold on;
        plot(0:90,T_P2_disp_ratio_uz(:,1),'m-','LineWidth',lw); hold on;
        plot(0:90,T_S_disp_ratio_uz(:,1),'c-','LineWidth',lw); hold on;
        xlabel('P wave Incidence angle');
        ylabel('z displacement ratio');
        legend('Reflected P wave', 'Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave','FontSize',15);
        title('z displacement ratio for open-pore acoustic-poro interface');
    elseif parameter_wave == 2
        figure(2);
        plot(0:90,R_P_disp_ratio_ux(:,1),'r-','LineWidth',lw); hold on;
        plot(0:90,R_S_disp_ratio_ux(:,1),'y-','LineWidth',lw); hold on;
        plot(0:90,T_P1_disp_ratio_ux(:,1),'g-','LineWidth',lw); hold on;
        plot(0:90,T_P2_disp_ratio_ux(:,1),'m-','LineWidth',lw); hold on;
        plot(0:90,T_S_disp_ratio_ux(:,1),'c-','LineWidth',lw); hold on;
        xlabel('P wave Incidence angle');
        ylabel('x displacement ratio');
        legend('Reflected P wave','Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave','FontSize',15);
        title('x displacement ratio for open-pore acoustic-poro interface');
    elseif parameter_wave == 3
        figure(3);
        plot(0:90,R_P_disp_ratio_u(:,1),'r-','LineWidth',lw); hold on;
        plot(0:90,R_S_disp_ratio_u(:,1),'y-','LineWidth',lw); hold on;
        plot(0:90,T_P1_disp_ratio_u(:,1),'g-','LineWidth',lw); hold on;
        plot(0:90,T_P2_disp_ratio_u(:,1),'m-','LineWidth',lw); hold on;
        plot(0:90,T_S_disp_ratio_u(:,1),'c-','LineWidth',lw); hold on;
        xlabel('P wave Incidence angle');
        ylabel('displacement ratio');
        legend('Reflected P wave','Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave','FontSize',15);
% legend( 'Transmitted P1 wave',  'Transmitted S wave');
%         title('displacement ratio for open-pore acoustic-poro interface');
    end
end
if parameter_type == 2
    figure(4);
    plot(0:90,R_P_energy_ratio(:,1),'r-','LineWidth',lw); hold on;
    plot(0:90,R_S_energy_ratio(:,1),'b-','LineWidth',lw); hold on;
    plot(0:90,T_P1_energy_ratio(:,1),'g-','LineWidth',lw); hold on;
    plot(0:90,T_P2_energy_ratio(:,1),'m-','LineWidth',lw); hold on;
    plot(0:90,T_S_energy_ratio(:,1),'c-','LineWidth',lw); hold on;
    plot(0:90,Energy_total_ratio(:,1),'y-','LineWidth',lw); hold on;
    xlabel('P wave Incidence angle');
    ylabel('Energy ratio');
    legend('Reflected P wave','Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave', 'Total energy ratio','FontSize',15);
%     title('open-pore acoustic-poro interface');
end


if parameter_type == 3
    figure(4);
    plot(0:90,A_R(:,1),'r-','LineWidth',lw); hold on;
    plot(0:90,A_S(:,1),'b-','LineWidth',lw); hold on;
    plot(0:90,B_PLUS(:,1),'g-','LineWidth',lw); hold on;
    plot(0:90,B_MINUS(:,1),'m-','LineWidth',lw); hold on;
    plot(0:90,B_S(:,1),'c-','LineWidth',lw); hold on;
    xlabel('P wave Incidence angle');
    ylabel('Energy ratio');
    legend('Reflected P wave','Reflected S wave', 'Transmitted P1 wave', 'Transmitted P2 wave', 'Transmitted S wave','FontSize',15);
%     title('Jun Yang: open-pore acoustic-poro interface');
end

R_T_disp = [R_P_disp_ratio_u R_S_disp_ratio_u T_P1_disp_ratio_u T_P2_disp_ratio_u T_S_disp_ratio_u];

