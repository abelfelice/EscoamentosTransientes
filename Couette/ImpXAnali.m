%% Comparação: Solução Analítica vs. Solução Implícita
%  Escoamento de Couette plano unidimensional transiente - ABEL
%  Plotagem dos perfis numéricos sobrepostos aos perfis analíticos

clc; clear; close all;

%% ========================  PARÂMETROS COMUNS  ========================

h  = 0.01;     % Distância entre as placas [m]
Vf = 1.0;      % Velocidade da placa superior [m/s]
V0 = 0.0;      % Velocidade da placa inferior [m/s]
Nu = 1e-6;     % Viscosidade cinemática (água) [m²/s]
N  = 100;      % Número de nós da malha

t0 = 0;        % Tempo inicial [s]
tf = 60;       % Tempo final [s]

t_plots = [0.5, 1, 2.5, 5, 10, 60]; % Instantes para comparação [s]

%% ====================  SOLUÇÃO IMPLÍCITA (NUMÉRICA)  =================
tic

dy = h / (N - 1);
CFL = 1;
dt = CFL * (dy^2) / Nu;
Nt = ceil((tf - t0) / dt);
y  = linspace(0, h, N)';
alpha = (Nu * dt) / dy^2;

% Pré-alocação da matriz de velocidades
u_num = zeros(N, Nt);
u_num(1, :)  = V0;
u_num(N, :)  = Vf;
u_num(:, 1)  = 0;

% Matriz tridiagonal (sistema implícito)
principal = -(1 + 2*alpha) * ones(N-2, 1);
superior  =  alpha * ones(N-3, 1);
inferior  =  alpha * ones(N-3, 1);
A = diag(principal) + diag(superior, 1) + diag(inferior, -1);

% Marcha temporal
B = zeros(N-2, 1);
for n = 1:Nt
    B(end) = B(end) - alpha * Vf;
    x = A \ B;
    B = -x;
    u_num(:, n) = [V0; x; Vf];
end

tempo_exec = toc;
fprintf('Tempo de execução da solução implícita: %.4f s\n', tempo_exec);

%% ====================  SOLUÇÃO ANALÍTICA  ============================

n_terms = 200; % Termos da série de Fourier

% Pré-alocação: cada coluna armazena o perfil analítico de um instante
u_anl = zeros(N, length(t_plots));

for i = 1:length(t_plots)
    t = t_plots(i);

    % Estado estacionário
    u_steady = Vf * (y / h);

    % Termo transiente (série de Fourier)
    u_transient = zeros(N, 1);
    for k = 1:n_terms
        lambda_k = k * pi / h;
        term_k = (2*Vf) / (k*pi) * (-1)^k * sin(lambda_k * y) ...
                 * exp(-Nu * lambda_k^2 * t);
        u_transient = u_transient + term_k;
    end

    u_anl(:, i) = u_steady + u_transient;
end

%% ====================  PLOTAGEM COMPARATIVA  =========================

% Paleta de cores (uma cor por instante de tempo)
cores = lines(length(t_plots));

figure('Name', 'Contínuo vs. Implícito', ...
       'Position', [100 100 900 700]);

hold on;
box on;
grid on;

legend_entries = cell(1, 2*length(t_plots));

for i = 1:length(t_plots)

    t = t_plots(i);

    % Índice correspondente ao instante desejado
    idx = min(ceil(t / dt), Nt);

    %% Perfil analítico — linha contínua
    plot(u_anl(:, i), y, '-', ...
        'Color', cores(i,:), ...
        'LineWidth', 2.5);

    legend_entries{2*i - 1} = ...
        sprintf('Contínuo   t = %.1f s', t);

    %% Perfil implícito — marcadores circulares

    step = max(1, floor(N/40));   % Aproximadamente 20 marcadores
    idx_pts = 1:step:N;

    plot(u_num(idx_pts, idx), y(idx_pts), 'o', ...
        'Color', cores(i,:), ...
        'MarkerSize', 6, ...
        'LineWidth', 1.5);

    legend_entries{2*i} = ...
        sprintf('Implícito   t = %.1f s', t);

end

%% Configurações do gráfico

title('Escoamento de Couette Transiente — Contínuo vs. Implícito', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');

xlabel('Velocidade u(y,t) [m/s]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('Distância y [m]', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

legend(legend_entries, ...
    'Location', 'southeast', ...
    'FontSize', 12);

axis([0 Vf 0 h]);

set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.2);

hold off;



