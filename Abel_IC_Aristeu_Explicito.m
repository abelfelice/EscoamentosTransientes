% Solução Explicita do escoamento de Couette plano unidimensional ABEL

clc; clear all; close all;

%% Modelo computacional Explícito:

% Parâmetros de entrada
tic
h = 0.01;     % Distância entre as duas placas [m]
N = 100;      % Numero de nós da malha
Nu = 1e-6;    % Viscosidade cinemática da água [m^2/s]
CFL = 0.5;

% Condições de contorno e de valor inicial

V0 = 0;  %[m/s]
Vf = 1;  %[m/s]
t0 = 0;  %[s]
tf = 60; %[s]

% Criação da malha na extensão do domínio

dy = h/(N-1);          % Distância de um nó ao outro na malha
dt = CFL*(dy^2)/Nu;    % Cálculo do valor do passo de tempo dt (fórmula)
y = linspace (0,h,N);  % Vetor dos intervalos do espaço
Nt = ceil((tf-t0)/dt); % Número total de passos de tempo

% Condição de estabilidade

if dt > CFL*(dy^2)/Nu
    error('O programa não atinge convergência numérica.');
end

%Implementação do modelo matemático

u = zeros(N,Nt); % Pré inicialização da matriz velocidade
u(1,:) = V0;     % Na linha 1 todos elementos são 0
u(N,:) = Vf;     % Na linha N todos os elementos são 1
u(:,1) = t0;     % Coluna 1 todos elementos são igual a 0


%Looping principal que dá resultados de velocidade a cada posição da matriz
for n = 1:Nt-1
    for i=2:N-1
        u(i,n+1) = u(i,n) + Nu*(dt/dy^2)*(u(i+1,n)-2*u(i,n)+u(i-1,n));
        %Obs: Na formula de cima foi necessária uma adaptação pois os
        %valores dos indices devem ser positivos
    end
end

%Isso é para na hora de plotar, o plot ser com a altura da placa
t1 = u(:,ceil(tf/dt));
t2 = u(:,ceil(30/dt));
t3 = u(:,ceil(20/dt));
t4 = u(:,ceil(15/dt));
t5 = u(:,ceil(10/dt));
t6 = u(:,ceil(7.5/dt));
t7 = u(:,ceil(5/dt));
t8 = u(:,ceil(2.5/dt));
t9 = u(:,ceil(1/dt));
t10 = u(:,ceil(0.5/dt));

% Vetor que armazena uma amostra das velocidades na
% posicao N-1 ao longo dos diferentes tempos
evolve = u(N-1,:);

%Vetor com Nt valores
aux = linspace(0,tf,Nt);

%% Plotar resultado

%Evolução dos perfis de velocidade:
% figure('Position', [100 100 1200 500]);
% subplot(1,2,1);
hold on
plot(t1,y,LineWidth=1.5);
plot(t2,y,LineWidth=1.5);
plot(t3,y,LineWidth=1.5);
plot(t4,y,LineWidth=1.5);
plot(t5,y,LineWidth=1.5);
plot(t6,y,LineWidth=1.5);
plot(t7,y,LineWidth=1.5);
plot(t8,y,LineWidth=1.5);
plot(t9,y,LineWidth=1.5);
plot(t10,y,LineWidth=1.5);
legend("t = Tempo final","t = 30","t = 20","t = 15","t = 10", ...
    "t = 7.5","t = 5","t = 2.5","t = 1","t = 0.5",'Location','best');
xlabel('Velocidade u [m/s]');
ylabel('Distância y [m]');
title('Perfil de Velocidade do Escoamento de Couette');
grid on;
axis([0 1 0 0.01]);
hold off

% % Evolução das velocidades na posição N-1 ao longo do tempo
% subplot(1,2,2);
% plot(aux,evolve,LineWidth=1.5);
% grid on;
% xlabel('Tempo total de simulação [s]');
% ylabel('Velocidade u [m/s]');
% title('Evolução da velocidade em um mesmo ponto ao longo dos tempos');

tempo_execucao = toc; % Calcula o tempo de execução
fprintf('Tempo real de execução: %.6f segundos\n', tempo_execucao);

% Animação do perfil de velocidade ao longo do tempo
fig = figure('Name','Animação do Perfil de Velocidade','Position',[100 100 600 500]);

ax  = axes('Parent',fig);
hLn = plot(ax,u(:,1),y,'b','LineWidth',1.5);   % só cria 1 vez
xlabel(ax,'Velocidade u [m/s]');
ylabel(ax,'Distância y [m]');
title(ax,'Perfil de Velocidade');
axis(ax,[0 1 0 h]);
grid(ax,'on');

pauseTime = 0.02;   % 0,15 s ≈ 6–7 fps (ajuste a gosto)

for n = 1:10:Nt
    if ~isvalid(fig); break;
    end
    set(hLn,'XData',u(:,n));  % atualiza a linha
    title(ax,sprintf('Perfil de Velocidade – Tempo = %.2f s',(n)*dt));
    drawnow
    pause(pauseTime);         % controla a velocidade
end

% Animação do perfil de velocidade ao longo do tempo e Gravação MP4

% Configuração do objeto VideoWriter
nome_video = 'Animacao_Escoamento_Couette.mp4';
v = VideoWriter(nome_video, 'MPEG-4');
v.FrameRate = 30; % Define a taxa de quadros (fps). Ajuste conforme a fluidez desejada.
open(v); % Abre o arquivo de vídeo para gravação

fig = figure('Name','Animação do Perfil de Velocidade','Position',[100 100 600 500]);

ax  = axes('Parent',fig);
hLn = plot(ax,u(:,1),y,'b','LineWidth',1.5);   % Inicializa a curva do perfil
xlabel(ax,'Velocidade u [m/s]');
ylabel(ax,'Distância y [m]');
title(ax,'Perfil de Velocidade');
axis(ax,[0 1 0 h]);
grid(ax,'on');

% Note que a função pause() foi removida. O VideoWriter já dita a velocidade 
% do arquivo final a partir do parâmetro v.FrameRate e do salto no loop (1:10:Nt).

for n = 1:10:Nt
    if ~isvalid(fig)
        break; 
    end
    
    set(hLn,'XData',u(:,n));  % Atualiza os dados vetoriais da linha
    title(ax,sprintf('Perfil de Velocidade – Tempo = %.2f s',(n)*dt));
    drawnow;
    
    % Captura o quadro atual da figura e escreve no objeto de vídeo
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Fecha e salva o arquivo de vídeo de forma segura
if isvalid(fig) % Garante que só fecha se a janela não foi fechada abruptamente
    close(v);
    fprintf('Animação salva com sucesso no arquivo: %s\n', nome_video);
end