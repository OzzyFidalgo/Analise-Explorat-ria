close all; clear; clc;

% Script principal
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

% Configuração de limites para massa de água e ângulo
lb = [0.1, 10];  % Limite inferior: [massa de água, ângulo]
ub = [1, 80];    % Limite superior: [massa de água, ângulo]

% Configuração do algoritmo genético
options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 50, ...
    'MaxGenerations', 100, 'FunctionTolerance', 1e-6);

% Chamada do algoritmo genético
[x_opt, d_max] = ga(@(x) -calcular_distancia_maxima(x(1), x(2)), 2, [], [], [], [], lb, ub, [], options);

% Exibe os resultados
disp(['Massa de água ótima: ', num2str(x_opt(1)), ' kg']);
disp(['Ângulo ótimo: ', num2str(x_opt(2)), ' graus']);
disp(['Distância máxima horizontal: ', num2str(-d_max), ' m']);

% Função local calcular_distancia_maxima
function d_max = calcular_distancia_maxima(m_agua, theta)
    global m_garrafa V_garrafa A_bocal A_frontal P_inicial h_rampa rho_agua rho_ar P_atm C_a g
    m_0 = m_agua;

    % Ângulo em radianos
    theta_rad = deg2rad(theta);
    
    % Condições iniciais
    v_x = 0;    % Velocidade horizontal inicial
    v_y = 0;    % Velocidade vertical inicial
    x = 0;      % Posição horizontal inicial
    y = 0;      % Posição vertical inicial (altura da rampa)
    P_inicial = 300000;

    % Intervalo de tempo para integração numérica
    dt = 0.05;
    tempo_max = 10; % Tempo máximo de simulação em segundos
    t = 0;          % Inicializa o contador de tempo
    k = 1;          % Parâmetro de empuxo inicial
    % Loop de simulação de movimento na Rampa
    disp(['Angulo theta ', num2str(theta)]);
    while x <= h_rampa / tan(theta_rad) && t < tempo_max
        % Cálculo do empuxo
        E = k * (P_inicial - P_atm) * A_bocal / 2;
        
        % Cálculo da força de arrasto
        v = sqrt(v_x^2 + v_y^2);
        F_a = 0.5 * C_a * rho_ar * A_frontal * v^2;
        
        % Cálculo das forças
        F_x = E * cos(theta_rad) - F_a * cos(theta_rad) - (m_garrafa + m_agua) * g * sin(theta_rad) * cos(theta_rad);
        F_y = E * sin(theta_rad) - F_a * sin(theta_rad) - (m_garrafa + m_agua) * g * sin(theta_rad) * sin(theta_rad);
        
        % Atualização das velocidades
        v_x = v_x + (F_x / (m_garrafa + m_agua)) * dt;
        v_y = v_y + (F_y / (m_garrafa + m_agua)) * dt;
        
        % Atualização das posições
        x = x + v_x * dt;
        y = y + v_y * dt;

        % Atualização da pressão e velocidade da água
        U = sqrt((2 * (P_inicial - P_atm)) / rho_agua);
        P_inicial = 300000 * (1 - m_0 / (rho_agua * V_garrafa)) / (1 - m_agua / (rho_agua * V_garrafa));

        % Verificação se a água foi ejetada completamente
        if m_agua > 0
            m_agua = m_agua - rho_agua * A_bocal * U * dt;
            if m_agua < 0
                m_agua = 0;
            end
        else
            k = 0;  % Sem empuxo após ejetar toda a água
        end

        % Incrementa o tempo
        t = t + dt;
    end
    
    % Loop de simulação de movimento fora da rampa
    while y > 0 && t < tempo_max
        % Cálculo do empuxo
        E = k * (P_inicial - P_atm) * A_bocal / 2;
        
        % Cálculo da força de arrasto
        v = sqrt(v_x^2 + v_y^2);
        F_a = 0.5 * C_a * rho_ar * A_frontal * v^2;
        
        % Cálculo das forças
        F_x = E * (v_x / v) - F_a * (v_x / v);
        F_y = E * (v_y / v) - F_a * (v_y / v) - (m_garrafa + m_agua) * g;
        
        % Atualização das velocidades
        v_x = v_x + (F_x / (m_garrafa + m_agua)) * dt;
        v_y = v_y + (F_y / (m_garrafa + m_agua)) * dt;
        
        % Atualização das posições
        x = x + v_x * dt;
        y = y + v_y * dt;
        P_inicial = 300000 * (1 - m_0 / (rho_agua * V_garrafa)) / (1 - m_agua / (rho_agua * V_garrafa))
        % Atualização da pressão e velocidade da água
        if m_agua > 0
            U = sqrt((2 * (P_inicial - P_atm)) / rho_agua);
        else
            U = 0;
            k=0;
        end
        
        % Verificação se a água foi ejetada completamente
        if m_agua > 0
            m_agua = m_agua - rho_agua * A_bocal * U * dt;
            if m_agua < 0
                m_agua = 0;
            end
        else
            k = 0;  % Sem empuxo após ejetar toda a água
        end
        
        % Incrementa o tempo
        t = t + dt;
    end
    
    % Distância máxima horizontal percorrida
    d_max = x;
end