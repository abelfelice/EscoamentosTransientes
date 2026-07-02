% Solução Gauss-Seidel do Escoamento de Poiseuille Bidimensional Transiente Laminar - ABEL

clc; clear all; close all;

%% Parâmetros físicos e geométricos
tf = 60;              % Tempo total [s]
h = 0.01;             % Altura [m]
L = h;                % Comprimento [m]
N = 100;              % Nós verticais
M = N;                % Nós horizontais
nu = 1e-6;            % Viscosidade cinemática [m^2/s]
rho = 1000;           % Massa específica [kg/m^3]
gradp = -1e1;         % Gradiente de pressão [Pa/m]
CFL = 1;              % CFL para controle do passo de tempo

%% Malha e parâmetros numéricos
delta = h / (N-1);            
dt = CFL * (delta^2) / nu;    
Nt = ceil(tf / dt);           
y = linspace(0, h, N);        
z = linspace(0, L, M);        
c = (1/rho)* gradp;             

%% Inicialização do campo de velocidades
u = zeros(N, M, Nt);           % u(j,k,n): j = vertical, k = horizontal, n = tempo

%% Coeficientes do método implícito (Gauss-Seidel)
a0 = 1 + 4 * (nu * dt / delta^2);
ae = -nu * dt / delta^2;
ad = -nu * dt / delta^2;
ab = -nu * dt / delta^2;
ac = -nu * dt / delta^2;

%% Parâmetros do Gauss-Seidel
max_iter = 1e4;
tol = 1e-8;

%% Loop principal no tempo
for n = 2:Nt
    % Inicializa próxima solução com o valor anterior
    u(:, :, n) = u(:, :, n-1);  
    
    for iter = 1:max_iter
        u_old = u(:, :, n); % Armazena para verificar convergência
        
        for j = 2:N-1
            for k = 2:M-1
                b = -c * dt + u(j, k, n-1);
                u(j, k, n) = (b - (ae * u(j-1, k, n) + ad * u(j+1, k, n) + ...
                                  ab * u(j, k-1, n) + ac * u(j, k+1, n))) / a0;
            end
        end
        
        % Critério de convergência
        erro = max(max(abs(u(:, :, n) - u_old)));
        if erro < tol
            break;
        end
    end
end

[Z, Y] = meshgrid(z, y);  % Use os vetores z e y completos que incluem fronteiras


%% Solução analítica

% Parâmetros Físicos
mu = rho * nu;    % Viscosidade dinâmica [Pa.s]
G = -gradp;       % Gradiente de pressão G = -dp/dx

% Parâmetros da Série
N_series = 100;   % Número de termos para a série. 

% Inicializa matriz da solução analítica
%u_analytical = zeros(N, M);

% Primeiro termo da equação
term1 = (G / (2*mu)) * Y .* (h - Y);

% Segundo termo (somatório da série infinita)
sum_series = zeros(N, M);
for n = 1:N_series
    beta_n = (2*n - 1) * pi / h;
    

    % No entanto, para valores usuais de N_series, o cálculo direto funciona.
    numerator_sinh = sinh(beta_n * Z) + sinh(beta_n * (L - Z));
    denominator_sinh = sinh(beta_n * L);
    
    sin_term = sin(beta_n * Y);
    
    term_n = (1 / (2*n - 1)^3) * (numerator_sinh / denominator_sinh) .* sin_term;
    sum_series = sum_series + term_n;
end

% Combina os termos para a solução final
u_analytical = term1 - (4*G*h^2 / (mu * pi^3)) * sum_series;

% Garante condições de contorno de não escorregamento 
u_analytical(1, :) = 0;
u_analytical(end, :) = 0;
u_analytical(:, 1) = 0;
u_analytical(:, end) = 0;

%% COMPARAÇÃO E PLOTAGEM DOS RESULTADOS 

u_numerical_ss = u(:, :, end); % Solução numérica em estado estacionário

figure('Name', 'Comparação de Resultados: Numérico vs. Analítico', 'Position', [100, 100, 1500, 500]);

% Subplot 1: Solução Numérica
subplot(1, 2, 1);
surf(Z, Y, u_numerical_ss);
shading interp;
xlabel('z [m]');
ylabel('y [m]');
zlabel('u [m/s]');
title('Solução Numérica (Gauss-Seidel)');
colorbar;
colormap(turbo);
view(45, 30);
axis tight;

% Subplot 2: Solução Analítica
subplot(1, 2, 2);
surf(Z, Y, u_analytical);
shading interp;
xlabel('z [m]');
ylabel('y [m]');
zlabel('u [m/s]');
title('Solução Analítica');
colorbar;
colormap(turbo);
view(45, 30);
axis tight;

% Título geral
sgtitle(sprintf('Comparação de perfis - t = %.2f s (Estado Estacionário)', tf), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Gráficos de Corte - Comparação Numérico vs Analítico

% Criar nova figura para os perfis de velocidade
figure('Name', 'Perfis de Velocidade - Cortes', 'Position', [100, 100, 1400, 600]);

% CORTE 1: Perfil ao longo de y (z fixo no meio do canal) 
k_meio = round(M/2);  % Índice do meio em z
z_corte = z(k_meio);  % Posição z do corte

subplot(1, 2, 1);
plot(y, u_numerical_ss(:, k_meio), 'b-', 'LineWidth', 2, 'DisplayName', 'Numérico');
hold on;
plot(y, u_analytical(:, k_meio), 'r--', 'LineWidth', 2, 'DisplayName', 'Analítico');
grid on;
xlabel('y [m]', 'FontSize', 12);
ylabel('u [m/s]', 'FontSize', 12);
title(sprintf('Perfil de Velocidade em z = 0.005 m (meio do canal)'), 'FontSize', 14);
legend('Location', 'north', 'FontSize', 11);
set(gca, 'FontSize', 11);

% Adicionar anotação do erro neste corte
erro_corte_y = abs(u_numerical_ss(:, k_meio) - u_analytical(:, k_meio));
erro_max_corte_y = max(erro_corte_y);
text(0.5*h, 0.1*max(u_analytical(:, k_meio)), ...
    sprintf('Erro máx: %.2e m/s', erro_max_corte_y), ...
    'FontSize', 10, 'BackgroundColor', 'white');

% CORTE 2: Perfil ao longo de z (y fixo no meio do canal)
j_meio = round(N/2);  % Índice do meio em y
y_corte = y(j_meio);  % Posição y do corte

subplot(1, 2, 2);
plot(z, u_numerical_ss(j_meio, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Numérico');
hold on;
plot(z, u_analytical(j_meio, :), 'r--', 'LineWidth', 2, 'DisplayName', 'Analítico');
grid on;
xlabel('z [m]', 'FontSize', 12);
ylabel('u [m/s]', 'FontSize', 12);
title(sprintf('Perfil de Velocidade em y = 0.005 m (meio do canal)'), 'FontSize', 14);
legend('Location', 'north', 'FontSize', 11);
set(gca, 'FontSize', 11);

% Adicionar anotação do erro neste corte
erro_corte_z = abs(u_numerical_ss(j_meio, :) - u_analytical(j_meio, :));
erro_max_corte_z = max(erro_corte_z);
text(0.5*L, 0.1*max(u_analytical(j_meio, :)), ...
    sprintf('Erro máx: %.2e m/s', erro_max_corte_z), ...
    'FontSize', 10, 'BackgroundColor', 'white');

% Título geral
sgtitle(sprintf('Comparação em corte - t = %.2f s (Estado Estacionário)', tf), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Erro 

% Erro Absoluto
erro_abs = abs(u_numerical_ss - u_analytical);

% Exibição do erro quantitativo no console
max_erro_abs = max(erro_abs(:));
max_vel_analytical = max(u_analytical(:));
erro_relativo_percentual = (max_erro_abs / max_vel_analytical) * 100;

fprintf(' Análise de Verificação \n');
fprintf('Velocidade máxima (Numérica):   %.6f m/s\n', max(u_numerical_ss(:)));
fprintf('Velocidade máxima (Analítica):  %.6f m/s\n', max_vel_analytical);
fprintf('---------------------------------\n');
fprintf('Erro absoluto máximo:           %.4e m/s\n', max_erro_abs);
fprintf('Erro relativo máximo (%%):       %.4f %%\n', erro_relativo_percentual);

maior_linha_u = max(u(:,:,end));
maior_valor_u = max(max(u(:,:,end)));
Re = (maior_valor_u*h)/nu;





%% ================= COMPARAÇÃO E PLOTAGEM DOS RESULTADOS ==============

u_numerical_ss = u(:, :, end); % Solução numérica em estado estacionário

%% =====================================================================
%% Figura 1 — Superfícies 3D
%% =====================================================================

figure('Name', 'Comparação: Numérico vs Contínuo', ...
       'Position', [100, 100, 1700, 700]);

%% ----------------- Solução Numérica -----------------

subplot(1,2,1);

surf(Z, Y, u_numerical_ss);

shading interp;
colormap(turbo);
colorbar;

view(45,30);
axis tight;

xlabel('z [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('y [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

zlabel('u [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title('Solução Numérica (Gauss-Seidel)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

%% ----------------- Solução Analítica -----------------

subplot(1,2,2);

surf(Z, Y, u_analytical);

shading interp;
colormap(turbo);
colorbar;

view(45,30);
axis tight;

xlabel('z [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('y [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

zlabel('u [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title('Solução Analítica', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

%% ----------------- Título Geral -----------------

sgtitle(sprintf('Comparação de Perfis — t = %.2f s (Estado Estacionário)', tf), ...
    'FontSize', 22, ...
    'FontWeight', 'bold');

%% =====================================================================
%% Figura 2 — Perfis de Velocidade em Corte
%% =====================================================================

figure('Name', 'Perfis de Velocidade em Corte', ...
       'Position', [100, 100, 1600, 700]);

%% ----------------- CORTE EM y -----------------

k_meio = round(M/2);
z_corte = z(k_meio);

subplot(1,2,1);

plot(y, u_numerical_ss(:, k_meio), ...
    'b-', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Numérico');

hold on;

plot(y, u_analytical(:, k_meio), ...
    'r--', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Contínuo');

grid on;
box on;

xlabel('y [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('u [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title(sprintf('Perfil de Velocidade em z = 0.005 m', z_corte), ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

legend('Location', 'north', ...
    'FontSize', 13);

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

%% Erro no corte em y

erro_corte_y = abs(u_numerical_ss(:, k_meio) - u_analytical(:, k_meio));
erro_max_corte_y = max(erro_corte_y);

text(0.45*h, 0.08*max(u_analytical(:, k_meio)), ...
    sprintf('Erro máx = %.2e m/s', erro_max_corte_y), ...
    'FontSize', 12, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');

hold off;

%% ----------------- CORTE EM z -----------------

j_meio = round(N/2);
y_corte = y(j_meio);

subplot(1,2,2);

plot(z, u_numerical_ss(j_meio, :), ...
    'b-', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Numérico');

hold on;

plot(z, u_analytical(j_meio, :), ...
    'r--', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Contínuo');

grid on;
box on;

xlabel('z [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('u [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title(sprintf('Perfil de Velocidade em y = %.3f m', y_corte), ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

legend('Location', 'north', ...
    'FontSize', 13);

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

%% Erro no corte em z

erro_corte_z = abs(u_numerical_ss(j_meio, :) - u_analytical(j_meio, :));
erro_max_corte_z = max(erro_corte_z);

text(0.45*L, 0.08*max(u_analytical(j_meio, :)), ...
    sprintf('Erro máx = %.2e m/s', erro_max_corte_z), ...
    'FontSize', 12, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');

hold off;

%% ----------------- Título Geral -----------------

sgtitle(sprintf('Comparação dos Perfis em Corte — t = %.2f s', tf), ...
    'FontSize', 22, ...
    'FontWeight', 'bold');