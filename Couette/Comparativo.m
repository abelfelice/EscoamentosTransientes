% Solução Implícita vs Explícita do escoamento de Couette plano unidimensional ABEL

clc; clear; close all;

%% Criação de estrutura de armazenamento de dados

Nmax = 100;     % Maior número de nós usado neste programa
Ns = 5:1:Nmax;  % Vetor com diferentes números de nós para testar
numNs = length(Ns);

% Estrutura para armazenar os resultados
resultados = struct('N', [], 'y', [], 'norma_L2', [], 'dy', [], 'invref', []);

%% Parâmetros de entrada (comuns)

h  = 0.01;    % Distância entre as duas placas [m]
Nu = 1e-6;    % Viscosidade cinemática da água [m²/s]
CFL = 0.5;

% Condições de contorno

V0 = 0;  % [m/s]
Vf = 1;  % [m/s]
t0 = 0;  % [s]
tf = 60; % [s]

% Instantes de tempo para comparação
t_plots = [0.5, 1.0, 2.5, 5.0, 10.0, 60.0];
n_plots = length(t_plots);

for idx = 1:numNs
    N = Ns(idx);

    % Criação da malha na extensão do domínio

    dy = h/(N-1);
    dt = CFL*(dy^2)/Nu;
    Nt = ceil((tf-t0)/dt);
    y  = linspace(0, h, N);
    alpha = (Nu*dt)/dy^2;
    invref = 1/dy;

    % Condição de estabilidade
    if dt > CFL*(dy^2)/Nu
        error('O programa não atinge convergência numérica.');
    end

    %% Método Implícito

    u_imp = zeros(N, Nt);
    u_imp(1,:) = V0;
    u_imp(N,:) = Vf;
    u_imp(:,1) = t0;

    principal = -(1 + 2*alpha) * ones(N-2, 1);
    superior  =  alpha * ones(N-3, 1);
    inferior  =  alpha * ones(N-3, 1);
    A = diag(principal) + diag(superior, 1) + diag(inferior, -1);

    B = zeros(N-2, 1);
    for n = 1:Nt
        B(end) = B(end) - alpha*Vf;
        x = A \ B;
        B = -x;
        u_imp(:, n) = [V0; x; Vf];
    end

    %% Método Explícito

    u_exp = zeros(N, Nt);
    u_exp(1,:) = V0;
    u_exp(N,:) = Vf;
    u_exp(:,1) = t0;

    for n = 1:Nt-1
        for i = 2:N-1
            u_exp(i,n+1) = u_exp(i,n) + Nu*(dt/dy^2)*(u_exp(i+1,n) - 2*u_exp(i,n) + u_exp(i-1,n));
        end
    end

    %% Extração dos perfis nos instantes desejados

    perfis_imp = zeros(N, n_plots);
    perfis_exp = zeros(N, n_plots);

    for k = 1:n_plots
        col = min(ceil(t_plots(k)/dt), Nt);
        perfis_imp(:, k) = u_imp(:, col);
        perfis_exp(:, k) = u_exp(:, col);
    end

    %% Cálculo da Norma L2 (entre implícito e explícito em t = 0.5 s)

    col_ref = min(ceil(t_plots(1)/dt), Nt);
    diff_vec = u_imp(:, col_ref) - u_exp(:, col_ref);
    norma_L2       = (1/(N-2)) * sqrt(sum(diff_vec.^2));
    norma_Linfinito = max(abs(diff_vec));

    %% Armazenar resultados

    resultados(idx).N        = N;
    resultados(idx).y        = y;
    resultados(idx).norma_L2 = norma_L2;
    resultados(idx).dy       = dy;
    resultados(idx).invref   = invref;

    fprintf('N = %3d  |  Norma L2 = %.6e  |  Norma L∞ = %.6e\n', N, norma_L2, norma_Linfinito);
end

%% Vetores auxiliares para o gráfico de convergência

invref_valores   = [resultados.invref];
norma_L2_valores = [resultados.norma_L2];

%% ========================  PLOTAGEM  =================================

cores = lines(n_plots);

%% Janela 1: Perfis de velocidade
figure('Position', [50 50 850 650]);

hold on;
box on;
grid on;

h_lines = gobjects(2*n_plots, 1);

for k = 1:n_plots

    % Implícito — linha contínua
    h_lines(2*k-1) = plot(perfis_imp(:,k), y, '-', ...
        'Color', cores(k,:), ...
        'LineWidth', 2.5);

    % Explícito — marcadores circulares
    step = max(1, floor(N/40));
    idx_pts = 1:step:N;

    h_lines(2*k) = plot(perfis_exp(idx_pts,k), y(idx_pts), 'o', ...
        'Color', cores(k,:), ...
        'MarkerSize', 6, ...
        'LineWidth', 1.5);
end

%% Legenda

leg_str = cell(2*n_plots, 1);

for k = 1:n_plots
    leg_str{2*k-1} = sprintf('Implícito   t = %.1f s', t_plots(k));
    leg_str{2*k}   = sprintf('Explícito   t = %.1f s', t_plots(k));
end

legend(h_lines, leg_str, ...
    'Location', 'southeast', ...
    'FontSize', 12);

%% Configurações dos eixos

xlabel('Velocidade u [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('Distância y [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title('Perfis de Velocidade — Implícito vs. Explícito', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

axis([0 Vf 0 h]);

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

hold off;

%% =====================================================================
%% Janela 2: Norma L2 vs refinamento
%% =====================================================================

figure('Position', [950 50 850 650]);

plot(invref_valores, norma_L2_valores, ...
    'b-o', ...
    'LineWidth', 2.0, ...
    'MarkerSize', 6);

grid on;
box on;

xlabel('1 / \Deltay  (inverso do refinamento)', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('Norma L_2', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

title('Norma L_2 (Implícito − Explícito) vs. Refinamento', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);
