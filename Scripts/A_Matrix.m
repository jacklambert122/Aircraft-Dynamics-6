%% Author: Jack Lambert
% ASEN 3128
% Homework 6
% Purpose: To keep all constants in one function so they are not defined
% more than once and then Compute the A Matrix and provide Trim States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,theta0,u0] = A_Matrix()
%% Airplane Parameters
% Nondimensional Derivatives
% Table 6.1 -
Cx = [-.108, .2193, 0, 0];
Cz = [-.106, -4.92, -5.921, 5.896];
Cm = [.1043, -1.023, -23.92, -6.314];

% Table E.1 B747 Case 3
Alt = 40000*(0.3048); % Altitude [ft] -> [m]
[T, a, P, rho] = atmosisa(Alt); % Standard Atmosphere Properties at Alt.
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
cbar = 27.31*(0.3048); % Mean Chord Length [ft]->[m]
S = 5500*(0.3048)^2; % Surface Area [ft^2]->[m^2]
g = 9.81; % Gravity Constant [m/s^2]
m = W/g; % Mass of Plane [kg]

%% Trim States
Vel = 871*(0.3048);% Velocity [ft/s] -> [m/s]
u0 = Vel; % Initial Velocity in x-coord - Stability Axis Frame (Trim State)
theta0 = 0; % Initial Pitch Angle [deg]
Cw0 = W/(.5*rho*S*u0^2); 
%% Function that Computes Dimensional Derivatives from Non-Dimenional derivatives
[X, Z, M ] = NonDimLong(rho,u0,S,W,theta0,Cx,Cz,Cm,cbar); 

%% State Variable Matrix A
deltaU_dot = [X(1)/m, X(2)/m, 0, -g*cosd(theta0)];
w_dot = [Z(1)/(m-Z(4)), Z(2)/(m-Z(4)), (Z(3)+m*u0)/(m-Z(4)), (-W*sind(theta0))/(m-Z(4))];
q_dot = [(1/Iy)*(M(1) + ((M(4)*Z(1))/(m-Z(4))) ),...
        (1/Iy)*(M(2) + ((M(4)*Z(2))/(m-Z(4))) ),...
        (1/Iy)*(M(3) + ((M(4)*(Z(3)+m*u0))/(m-Z(4))) ),...
        -((M(4)*W*sind(theta0))/(Iy*(m-Z(4))))];
deltaTheta_dot = [0, 0, 1, 0];

A = [deltaU_dot;w_dot;q_dot;deltaTheta_dot];

end