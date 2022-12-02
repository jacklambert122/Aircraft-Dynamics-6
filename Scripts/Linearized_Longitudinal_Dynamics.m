%% Author: Jack Lambert
% ASEN 3128
% Problem 5
% Purpose: Function for ODE45 to call to calculate the State variables xE,
% zE, u_dot, w_dot, q_dot, and theta_dot. This function uses the simplified
% assumptions for the Linearized Longitudinal Dynamics Set
% Last Edited: 3/11/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dzdt] = Linearized_Longitudinal_Dynamics(t,z)

X_E = z(1); % x-position, Inerital Frame 
Z_E = z(2); % z-position, Inerital Frame 
u_dot = z(3); % x-component of Velocity, Body Frame
w_dot = z(4); % z-component of Velocity, Body Frame
q_dot = z(5); % y-component of Angular Velocity, Body Frame
theta_dot = z(6); % Pitch Angle 

%% State Variable Matrix for Linearized Longitudinal Set
[A,theta0,u0] = A_Matrix(); % A matrix function based on plane and parameters
State = [u_dot, w_dot, q_dot, theta_dot]'; % Couple State Variables in Long. Set
var = A*State; % % Couple State Variables in Long. Set
%% Solving for Inertial Position
dzdt(1) = u_dot*cosd(theta0) + w_dot*sind(theta0) - u0*theta_dot*sind(theta0); % xE
dzdt(2) = -u_dot*sind(theta0) + w_dot*cosd(theta0)-u0*theta_dot*cosd(theta0); % zE
%% Solving for State Variables in the Linearized Longitudinal Set
dzdt(3) = var(1); % uE
dzdt(4) = var(2); % wE
dzdt(5) = var(3); % q
dzdt(6) = var(4); % theta

dzdt = dzdt'; % Inverts for ODE45   
end