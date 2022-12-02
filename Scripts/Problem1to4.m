%% Author: Jack Lambert
% ASEN 3128
% Homework 6
% Problems 1-4
% Last Edited: 3/11/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%% Problem 1 - B747 Dimensional Derivatives at 40,000 ft
% Nondimensional Derivatives
% Table 6.1 -
Cx = [-.108, .2193, 0, 0];
Cz = [-.106, -4.92, -5.921, 5.896];
Cm = [.1043, -1.023, -23.92, -6.314];

% Table E.1 B747 Case 3
Alt = 40000*(0.3048); % Altitude [ft] -> [m]
[T, a, P, rho] = atmosisa(Alt); % Standard Atmosphere Properties at Alt.
Vel = 871*(0.3048);% Velocity [ft/s] -> [m/s]
u0 = Vel; % Initial Velocity in x-coord - Stability Axis Frame
W = 6.366*10^5*4.44822; % Weight [lb]->[N]
Ix_PA = 1.82*10^7*1.35581795; % Moment of Interia x-PA [slug ft^2]-> [kg m^2]
Iy_PA = 3.31*10^7*1.35581795; % Moment of Interia y-PA [slug ft^2]-> [kg m^2]
Iz_PA = 4.97*10^7*1.35581795; % Moment of Interia z-PA [slug ft^2]-> [kg m^2]
Izx_PA = 9.70*10^5*1.35581795; % Moment of Interia zx-PA [slug ft^2]-> [kg m^2]
zeta = -2.4; % Angle between Stability Axis and PA [degrees] 
I = [Ix_PA, 0,-Izx_PA;...
    0, Iy_PA,0;...
    -Izx_PA, 0, Iz_PA]; % Inertia Matrix in PA
Q_PA_SA = [cosd(zeta), 0, -sind(zeta);...
    0, 1, 0;...
    sind(zeta), 0, cosd(zeta)]; % Transformation Matrix [PA-SA]
I_SA = Q_PA_SA * I * Q_PA_SA'; % MOI in Stability axis Frame
Ix = I_SA(1,1); % Moment of Interia x-SA [kg m^2]
Iy = I_SA(2,2); % Moment of Interia y-SA [kg m^2]
Iz = I_SA(3,3); % Moment of Interia z-SA [kg m^2]
Izx = (1/2)*(Ix-Iz)*sind(2*zeta)+Izx_PA*...
    (sind(zeta)^2-cosd(zeta)^2); % Moment of Interia zx-SA [kg m^2]

CD = .043; % Coefficient of Drag
theta0 = 0; % Initial Pitch Angle [deg]
cbar = 27.31*(0.3048); % Mean Chord Length [ft]->[m]
S = 5500*(0.3048)^2; % Surface Area [ft^2]->[m^2]
Cw0 = W/(.5*rho*S*u0^2); 
g = 9.81; % Gravity Constant [m/s^2]
m = W/g; % Mass of Plane [kg]

% Function that Computes Dimensional Derivatives from Non-Dimenional derivatives
[X, Z, M ] = NonDimLong(rho,u0,S,W,theta0,Cx,Cz,Cm,cbar); 
T = table(X',Z',M');
T.Properties.VariableNames = {'X' 'Z' 'M'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2 - A Matrix for Linearized Longitudinal Dynamics  
[A,theta0,u0] = A_Matrix(); % Function that Computes the A Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3 - Short Period and Phugoid Modes
[eigVec,eigVal] = eig(A);

modes = diag(eigVal);
max_real = max(abs(real(modes)));
% Short Mode has Larger Real Part
 j = 1;
 k = 1;
for i = 1:length(modes)
    if abs(real(modes(i))) == max_real
        SP_Mode(j) = modes(i); % Short Period Mode
        SP_vec(:,j) = eigVec(:,i); % Short Period Eigen Vec 
        j = j+1;
    else
        Phu_Mode(k) = modes(i); % Phugoid Mode
        Phu_vec(k,:) = eigVec(:,i); % Short Period Eigen Vec
        k = k+1;
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Natural Frequency and Dampeing Ratio

% Phugoid Mode

Wn_PM = ( real(Phu_Mode(1))^2+imag(Phu_Mode(1))^(2) )^(1/2); % Natural Frequency
zeta_PM = -real(Phu_Mode(1))/Wn_PM; % Dampening Coefficient
Period_PM = (2*pi) / imag(Phu_Mode(1)); % Period

% Short Period Mode

Wn_SP = ( real(SP_Mode(1))^2+imag(SP_Mode(1))^(2) )^(1/2); % Natural Frequency
zeta_SP = -real(SP_Mode(1))/Wn_SP; % Dampening Coefficient
Period_SP = 1 / imag(SP_Mode(1)); % Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 4 

% Linearized Short Period Approximation
A_Lin = [M(3)/Iy u0*M(2)/Iy; 1 0]; % State Variable Matrix

[eigVec_Lin,eigVal_Lin] = eig(A_Lin);
Mode_SP_Lin = diag(eigVal_Lin); % EigenValues

% Lanchester Approximation for Phugoid Period

Period_PM_Lan = pi*(2)^(1/2)*(u0/g);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




