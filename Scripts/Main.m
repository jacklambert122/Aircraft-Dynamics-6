%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Jack Lambert
% Dale Lawrence
% Aircraft Dynmaics Homework 6
% Problem 5
% Purpose: Sets Initial COnditions for each Pertubation Case and Calls ODE45
% to plot the State Variables vs time
% Date Modefied: 2/12/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Variable Allocation
%                     X_E = z(1); % z-position, Inerital Frame 
%                     Z_E = z(2); % z-position, Inerital Frame 
%                     u_dot = z(3); % x-component of Velocity, Body Frame
%                     z_dot = z(4); % x-component of Velocity, Body Frame
%                     q_dot = z(6); % Angular Velocity about the y-axis [rad/s]
%                     theta_dot = z(5); % Pitch Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inital Conditions:
%                     i = 1 ----> 10 [m/s] - U
%                     i = 2 ----> 10 [m/s] - W
%                     i = 3 ----> 0.1 [rad/s] - q
%                     i = 4 ----> 0.1 [rad] - Pitch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Conditions
c1 = [0 0 0 0]; % xE: Location in Inertial Coordinates [m]
c2 = [0 0 0 0]; % zE: Location in Inertial Coordinates [m]
c3 = [10 0 0 0]; % Delta U: x-comp, BF Interial Velocity [m/s]
c4 = [0 10 0 0]; % Delta W: z-comp, BF Interial Velocity [m/s]
c5 = [0 0 0.1 0]; % Delta q: y-comp, BF Angular Velocity [rad/s]
c6 = [0 0 0 0.1]; % Delta Theta: Pitch Angle

for i = 1:4
    condition{i}= [c1(i) c2(i) c3(i) c4(i) c5(i) c6(i)]; 
end

%% State Variables vs. Time
t_Phugoid = [0 1200]; % Longer time to see Phugoid Mode
t_ShortP = [0 7]; % Shorter time to see Short Phase Mode

string = ["U = 10 [m/s]","W = 10 [m/s]","q = 0.1 [rad/s]",...
    "theta = 0.1 [rad] "]; % Title for Varying IC's
% Phugoid Response (Longer Time)
for i = 1:4
    % Calling ODE45 
    [t,z] = ode45('Linearized_Longitudinal_Dynamics',t_Phugoid,condition{i});
    
    % Plotting Conditions
    figure
    % U_E vs time
    subplot(4,1,1)
    plot(t ,z(:,3),'Linewidth',1)
    tit = sprintf('%s %s %s','State Variable of a B 747,',string(i),'(Phugoid Mode)');
    title(tit)
    ylabel('u_E [m/s]')
    
    
    % W_E vs time
    subplot(4,1,2)
    plot(t ,z(:,4),'Linewidth',1)
    ylabel('w_E [m/s]')
    
    % q vs time
    subplot(4,1,3)
    plot(t ,z(:,5),'Linewidth',1)
    ylabel('q [rad/s]')
    
    % Theta vs time
    subplot(4,1,4)
    plot(t ,z(:,6),'Linewidth',1)
    ylabel('Theta [rad]')
    xlabel('Time [s]')

end

% Short Phase Response (shorter Time)
for i = 1:4
    % Calling ODE45 
    [t,z] = ode45('Linearized_Longitudinal_Dynamics',t_ShortP,condition{i});
    
    % Plotting Conditions
    figure
    % U_E vs time
    subplot(4,1,1)
    plot(t ,z(:,3),'Linewidth',1)
    tit = sprintf('%s %s %s','State Variable of a B 747,',string(i),'(Short Phase Mode)');
    title(tit)
    ylabel('u_E [m/s]')
    
    
    % W_E vs time
    subplot(4,1,2)
    plot(t ,z(:,4),'Linewidth',1)
    ylabel('w_E [m/s]')
    
    % q vs time
    subplot(4,1,3)
    plot(t ,z(:,5),'Linewidth',1)
    ylabel('q [rad/s]')
    
    % Theta vs time
    subplot(4,1,4)
    plot(t ,z(:,6),'Linewidth',1)
    ylabel('Theta [rad]')
    xlabel('Time [s]')

end

%% Plotting Position
for i = 1:4
    [t,z] = ode45('Linearized_Longitudinal_Dynamics',t_Phugoid,condition{i});

    % xE vs zE
    figure
    plot(z(:,1) ,z(:,2),'Linewidth',1)
    tit = sprintf('%s %s %s','Position of a B747,',string(i));
    title(tit)
    xlabel('xE [m]')
    ylabel('zE [m]')
    axis equal
end