% Undergraduate Student: Arturo Burgos
% Professor: Aristeu da Silveira Neto
% Federal University of Uberlândia - UFU, Fluid Mechanics Laboratory - MFLab, Block 5P, Uberlândia, MG, Brazil


close all
clear all
clc

% Disk and flat surface

% Variables

Patm = 101250; % Atmospheric pressure [Pa]
mp = 9; % Disk mass [kg]
g = 9.81; % Gravity [m/s^2]
dp = 0.02; % Disk diameter [m]
Ap = pi * ( dp/2 )^2; % Disk area [m^2]
mu = 17.72e-6; % Dynamic Viscosity [N*s/m^2]

Cp = 1004.8; % Specific heat [J/kg*K]
R = 287; % Gas constant [J/kg*K]
mu_oil = 2.9e-1; % Oil viscosity [N*s/m^2]
e = 0.8e-6; % Oil gap [m]
hp = 2e-3; % Disk thickness [m]
smax = 0.00092; % Maximum amplitude of the surface [m]
f = 9; % Vibration frequency of the surface [Hz]




% Initial Conditions

L = 0.01; % Initial lenght of the disk [m]
sb = 0; % Initial lenght of the surface [m] 
rho = 1.22; % Initial specific mass [kg/m^3]

t = 0; % Initial time [s]
tf = 3; % Final time [s]
dt = 0.0000004; % Time step [s], in the beginning it was 0.00000005
np = (tf-t)/dt; % Number of steps

T = 15 + 273; % Initial Temperature [K]
V = 0; % Initial disk speed [m/s]
Vb = 0; % Initial surface speed [m/s]  

mg = rho*Ap*L; % Mass of gas contained in the container [kg] (SHOULD BE CLOSE TO CONSTANT)




% Creation of the initial and final dy/dt arrays

yn_1 = [ V rho L T mg sb Vb]; % Initial array
yn = [ 0 0 0 0 0 0 0 ]; % Final array




% Criation of the matrix m 

m = zeros(8, np); 


for i = 1:np
    
    m(1,i) = t; % the matrix m gets in the row 1 the time
    m(2,i) = yn_1(1); % a matrix m gets in the row 2 the speed of the disk 
    m(3,i) = yn_1(2); % a matrix m gets in the row 3 the rho  
    m(4,i) = yn_1(3); % a matrix m gets in the row 4 the L 
    m(5,i) = yn_1(4); % a matrix m gets in the row 5 the T 
    m(6,i) = yn_1(5); % a matrix m gets in the row 6 the mg 
    m(7,i) = yn_1(6); % a matrix m gets in the row 7 the sb 
    m(8,i) = yn_1(7); % a matrix m gets in the row 8 the Vb 
    
    % yn_1(1) = V
    % yn_1(2) = rho
    % yn_1(3) = L
    % yn_1(4) = T
    % yn_1(5) = mg
    % yn_1(6) = sb
    % yn_1(7) = Vb
    
        
        yn(1) = yn_1(1) + dt * ((((yn_1(2) * R * yn_1(4) - Patm) * Ap) / mp) - g - (mu_oil * yn_1(1) * pi * dp * hp / (e*mp)));
        yn(2) = yn_1(2) - (dt * (yn_1(2) * (yn(1) - yn_1(7)))) / yn_1(3);
        yn(3) = yn_1(3) - dt * (yn_1(7) - yn(1));
        yn(4) = yn_1(4) +  dt * ((4/3 * mu * (((yn(1) - yn_1(7)) / yn(3))^2) / (yn(2) * Cp)) + (((yn(2) * R *yn_1(4)-Patm) * (yn(2) - yn_1(2))) / (( yn(2) * Cp) * (yn(2) *dt))));
        %yn(4) = yn_1(4);
        yn(5) = yn(2)*yn(3)*Ap;
        yn(6) = smax * sin(2 * pi * f * (t + pi/2));
        yn(7) = smax * 2 * pi * f * cos(2 * pi * f * t);
        
        
    t = t + dt;
    yn_1 = yn;

end

figure
plot(m(1,:),m(2,:))
title('Disk Speed');
xlabel('Time [s]');
ylabel('Speed [m/s]');
grid minor


figure
plot(m(1,:),m(3,:))
title('Disk Specific Mass Variation');
xlabel('Time [s]');
ylabel('\rho [kg/m^3]');
grid minor


figure
plot(m(1,:),m(4,:),m(1,:),m(7,:))
title('Height Variation');
xlabel('Time [s]');
ylabel('H [m]');
grid minor


figure
plot(m(1,:),m(5,:))
title('Temperature Variation');
xlabel('Time [s]');
ylabel('Temperature [K]');
grid minor
