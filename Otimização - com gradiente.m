close all; clear all; clc;
% Script principal - untitled3.m
global m_garrafa V_garrafa A_bocal A_frontal P_inicial h_rampa rho_agua rho_ar P_atm C_a g

% Inicialização das constantes
m_garrafa = 0.1;         % kg
V_garrafa = 0.003;       % m^3
A_bocal = 7.85398e-5;    % m^2
A_frontal = 7.853982e-3; % m^2
P_inicial = 300000;      % Pa
h_rampa = 1;             % m
rho_agua = 1000;         % kg/m^3
rho_ar = 1.4;            % kg/m^3
P_atm = 101325;          % Pa
C_a = 0.5;               % coeficiente de arrasto
g = 9.8;                 % m/s^2

% Configuração de otimização
options = optimset('Display','iter','Diagnostics','on','Algorithm','sqp');

% Ângulo da rampa (pode ser ajustado conforme necessário)
theta = 45; % graus

% Definindo os limites para a massa da água (por exemplo, entre 0.1 e 1 kg)
lb = 0.1;  % Limite inferior da massa de água
ub = 1;    % Limite superior da massa de água

% Chamada de otimização para encontrar a massa de água que maximiza a distância
[m_agua_opt, d_max] = fmincon(@(m_agua) -calcular_distancia_maxima(m_agua, theta), 0.5, [], [], [], [], lb, ub, [], options);

% Exibe o resultado
disp(['Massa de água ótima: ', num2str(m_agua_opt), ' kg']);
disp(['Distância máxima horizontal: ', num2str(-d_max), ' m']);

% Função local calcular_distancia_maxima
function d_max = calcular_distancia_maxima(m_agua, theta)
    global m_garrafa V_garrafa A_bocal A_frontal P_inicial h_rampa rho_agua rho_ar P_atm C_a g
    
    % Ângulo em radianos
    theta_rad = deg2rad(theta);
    
    % Condições iniciais
    v_x = 0;    % Velocidade horizontal inicial
    v_y = 0;    % Velocidade vertical inicial
    x = 0;      % Posição horizontal inicial
    y = h_rampa; % Posição vertical inicial (altura da rampa)
    
    % Intervalo de tempo para integração numérica
    dt = 0.05;
    tempo_max = 20; % Tempo máximo de simulação em segundos
    t = 0;          % Inicializa o contador de tempo
    
    % Loop de simulação de movimento
    while y >= 0 && t < tempo_max
        % Cálculo do empuxo
        E = (P_inicial - P_atm) * A_bocal;
        
        % Cálculo da força de arrasto
        v = sqrt(v_x^2 + v_y^2);  % Magnitude da velocidade
        F_a = 0.5 * C_a * rho_ar * A_frontal * v^2;
        
        % Evitar divisão por zero na direção das forças
        if v ~= 0
            F_x = E * cos(theta_rad) - F_a * (v_x / v);
            F_y = E * sin(theta_rad) - (m_garrafa + m_agua) * g - F_a * (v_y / v);
        else
            F_x = E * cos(theta_rad);
            F_y = E * sin(theta_rad) - (m_garrafa + m_agua) * g;
        end
        
        % Atualização das velocidades
        v_x = v_x + (F_x / (m_garrafa + m_agua)) * dt;
        v_y = v_y + (F_y / (m_garrafa + m_agua)) * dt;
        
        % Atualização das posições
        x = x + v_x * dt;
        y = y + v_y * dt;
        
        % Verificação se a água foi ejetada completamente
        if m_agua > 0
            m_agua = m_agua - rho_agua * A_bocal * v * dt;  % Perda de massa devido à ejeção da água
            if m_agua < 0
                m_agua = 0;
            end
        else
            E = 0;  % Sem empuxo após ejetar toda a água
        end
        
        % Incrementa o tempo
        t = t + dt;
    end
    
    % Distância máxima horizontal percorrida
    d_max = x;
end
