close all
clear all
clc

% Disco receptor e oscilador

% Dados Fixos

Patm = 101250; % pressão atmosférica [Pa]

mp = 9; % massa pistão [kg]

g = 9.81; % aceleração da gravidade [m/s^2]

dp = 0.02; % diametro do pistão [m]

Ap = pi * ( dp/2 )^2; % area do pistão [m^2]

mu = 17.72e-6; % viscosidade dinâmica [N*s/m^2]

Cp = 1004.8; % calor especifico [J/kg*K]

R = 287; % Constante do gás [J/kg*K]

mu_oleo = 2.9e-1; % Viscosidade do Oleo [N*s/m^2]

e = 0.8e-6; % Gap de Oleo [m]

hp = 2e-3; % Espessura do pistão [m]

smax = 0.00092;

% Dados Iniciais

L = 0.01; % Comprimento inicial receptor [m]

sb = 0; % Comprimento inicial oscilador [m] 

rho = 1.22; % massa específica inicial [kg/m^3]

t = 0; % tempo inicial [s]
tf = 3; % tempo final [s]
dt = 0.0000004; % passo de tempo[s] inicialmente 0.00000005

np = (tf-t)/dt; % numero de passos

T = 15 + 273; % Temperatura inicial [K]
V = 0; %Velociade inicial receptor [m/s]
Vb = 0; %Velocidade inicial do excitador [m/s]  

mg = rho*Ap*L; % massa de gás contida no recipiente [kg] (DEVERÁ APRESENTAR UM COMPORTAMENTO MAIS OU MENOS FIXO)

% Criação dos vetores dy/dt final e inicial

yn_1 = [ V rho L T mg sb Vb]; % vetor inicial
yn = [ 0 0 0 0 0 0 0 ]; % vetor final 

% Criação da matriz que em cada divisão armazena os valores das iterações

m = zeros(8, np); % matriz m de 4 linhas e np colunas

f= 9;
for i = 1:np
    
    m(1,i) = t; % a matriz recebe na linha 1 os tempos nos intervalos np
    m(2,i) = yn_1(1); % a matriz recebe na linha 2 as velocidades do receptor nos intervalos np
    m(3,i) = yn_1(2); % a matriz recebe na linha 3 os rho nos intervalos np 
    m(4,i) = yn_1(3); % a matriz recebe na linha 4 os L nos intervalos np
    m(5,i) = yn_1(4); % a matriz recebe na linha 5 os T nos intervalos np
    m(6,i) = yn_1(5); % a matriz recebe na linha 6 os mg nos intervalos np
    m(7,i) = yn_1(6); % a matriz recebe na linha 7 os sb nos intervalos np
    m(8,i) = yn_1(7); % a matriz recebe na linha 8 os Vb nos intervalos np
    
    % yn_1(1) = V
    % yn_1(2) = rho
    % yn_1(3) = L
    % yn_1(4) = T
        
        yn(1) = yn_1(1) + dt * ( ( (  ( yn_1(2) * R * yn_1(4) - Patm ) * Ap ) / mp ) - g - ( mu_oleo * yn_1(1) * pi * dp * hp / (e*mp) ) );
        yn(2) = yn_1(2) - ( dt * ( yn_1(2) * ( yn(1) - yn_1(7) ) ) ) / yn_1(3);
        yn(3) = yn_1(3) - dt * ( yn_1(7) - yn(1) );
        yn(4) = yn_1(4) +  dt * ( (4/3 * mu * ( ( ( yn(1) - yn_1(7) ) / yn(3) )^2 ) / ( yn(2) * Cp ) ) + ( ( (yn(2) * R *yn_1(4)-Patm) * ( yn(2) - yn_1(2) ) ) / (( yn(2) * Cp ) * ( yn(2) *dt )) ) ) ;
        %yn(4) = yn_1(4);
        yn(5) = yn(2)*yn(3)*Ap;
        yn(6) = smax * sin( 2 * pi * f * (t + pi/2));
        yn(7) = smax * 2 * pi * f * cos( 2 * pi * f * t);
        
        
    t = t + dt;
    yn_1 = yn;

end

figure
plot(m(1,:),m(2,:))
title(' Velocidade do Disco');
xlabel('Tempo [s]');
ylabel('Velocidade [m/s]');
grid minor

figure
plot(m(1,:),m(3,:))
title('Variação da Massa específica do Disco');
xlabel('Tempo [s]');
ylabel('\rho [kg/m^3]');
grid minor


figure
plot(m(1,:),m(4,:),m(1,:),m(7,:))
title('Variação da Altura');
xlabel('Tempo [s]');
ylabel('H [m]');
grid minor


figure
plot(m(1,:),m(5,:))
title('Variação da Temperatura');
xlabel('Tempo [s]');
ylabel('Temperatura [K]');
grid minor
