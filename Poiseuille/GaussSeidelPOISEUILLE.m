% Solução Gauss-Seidel do Escoamento de Poiseuille Bidimensional Transiente Laminar - ABEL

clc; clear all; close all;

%% Parâmetros físicos e geométricos
tic
tf = 60;              % Tempo total [s]
h = 0.01;             % Altura [m]
L = h;                % Comprimento [m]
N = 50;               % Nós verticais
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

%% Plotagem para ver o resultado final

[Z, Y] = meshgrid(z, y);  % Use os vetores z e y completos que incluem fronteiras

% ANIMAÇÃO FLUIDA DO ESCOAMENTO 

f = figure; % Uma janela única

for n = 1:Nt
    if ~ishandle(f) % Use ishandle(f) para checar se a janela está aberta
        disp('Animação encerrada pelo usuário.');
        break
    end %Encerra a simulação quando clicka na janela da figura para fechar

% Plot da animação do escoamento 
    surf(Z, Y, u(:,:,n), 'EdgeColor', 'none');
    view(45, 30);
    colormap(turbo);
    colorbar;
    title('Evolução da velociade com tempo');
    xlabel('Coordenada z (m)');
    ylabel('Coordenada y (m)');
    zlabel('Velocidade (m/s)');
    zlim([0 max(u(:))]);
    axis([0 L 0 h 0 1.1*max(u(:))]);
    shading interp;
    drawnow;     % Atualiza a tela para cada iteração
end

% %% Plotagem alternativa para um só passo de tempo
% 
% % Geração da malha para plotagem
% [X, Y] = meshgrid(z, y);
% 
% % Gráfico 3D da velocidade no tempo final
% figure;
% surf(X, Y, u(:, :, end));  % u(y,z) no tempo final
% shading interp;
% xlabel('z [m]');
% ylabel('y [m]');
% zlabel('u [m/s]');
% title('Perfil de Velocidade 3D - Escoamento de Poiseuille');
% colorbar;
% colormap("turbo");
% view(45, 30); % Ângulo de visualização
% 
% tempo_execucao = toc; % Calcula o tempo de execução
% fprintf('Tempo real de execução: %.6f segundos\n', tempo_execucao);

%% Cálculo de Parâmetros de Escoamento no Tempo Final
% Extração do campo de velocidade no tempo final
u_final = u(:, :, end);

% Busca do valor máximo global de velocidade na malha
u_max = max(u_final, [], 'all'); % O(N^2)

% O diâmetro hidráulico (Dh) para um duto de seção quadrada (h x L, com h=L) 
% é dado por Dh = 4A/P = 4(h^2)/(4h) = h.
Dh = h;

% Cálculo do Número de Reynolds global do escoamento
Re = (u_max * Dh) / nu;

%% Exibição dos Resultados
fprintf('\n--- Resultados Analíticos ---\n');
fprintf('Velocidade Máxima (u_max): %.5f m/s\n', u_max);
fprintf('Número de Reynolds (Re)  : %.2f\n', Re);
fprintf('Regime de Escoamento     : %s\n', evalc('if Re < 2300, disp(''Laminar''); else, disp(''Turbulento/Transição''); end'));

% =========================================================================
% ANIMAÇÃO FLUIDA DO ESCOAMENTO E EXPORTAÇÃO MP4
% =========================================================================

% Configuração do VideoWriter
nome_video = 'Animacao_Poiseuille_GaussSeidel.mp4';
v = VideoWriter(nome_video, 'MPEG-4');
v.FrameRate = 30; % Ajuste a taxa de quadros conforme a fluidez desejada
open(v);

f = figure('Name', 'Animação 3D - Poiseuille', 'Position', [100 100 800 600]);

% Inicialização do objeto gráfico surface fora do loop (Otimização de memória)
% Pegamos o quadro inicial u(:,:,1)
hSurf = surf(Z, Y, u(:,:,1), 'EdgeColor', 'none');
view(45, 30);
colormap(turbo);
colorbar;
xlabel('Coordenada z [m]');
ylabel('Coordenada y [m]');
zlabel('Velocidade u [m/s]');
shading interp;

% Fixação estrita dos eixos baseada no máximo global para evitar distorção no vídeo
vel_max_global = max(u(:));
zlim([0 vel_max_global]);
axis([0 L 0 h 0 1.1*vel_max_global]);

% Como Nt será muito grande para um dt numérico restritivo, pulamos alguns
% frames para o vídeo não ficar longo e pesado. O passo "step" dita isso.
step = 5; 

for n = 1:step:Nt
    if ~ishandle(f)
        disp('Animação abortada antes da conclusão da gravação.');
        break
    end
    
    % Atualização direta apenas dos dados z (elevação) e c (cor) do surface
    % Isso é computacionalmente muito mais leve do que usar surf() no loop
    set(hSurf, 'ZData', u(:,:,n), 'CData', u(:,:,n));
    
    % Atualiza o título com a métrica temporal real
    title(sprintf('Evolução da velocidade - Tempo: %.4f s', n*dt));
    drawnow;
    
    % Captura e escrita do frame
    frame = getframe(f);
    writeVideo(v, frame);
end

% Fechamento seguro do arquivo de vídeo
if ishandle(f)
    close(v);
    fprintf('Vídeo MP4 salvo com sucesso no diretório atual: %s\n', nome_video);
end
