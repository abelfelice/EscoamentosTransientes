% Solução Implicita do Escoamento de Poiseuille Bidimensional Transiente Laminar ABEL

clc; clear; close all;

%% Modelo computacional Implícito:

% Parâmetros de entrada
tic
tf = 60;       % Tempo total de simulação [s] 
h = 0.01;      % Distância vertical do tubo [m]
L = h;         % Distância horizontal do tubo [m]
N = 49;        % Número de nós verticais da malha
M = N;         % Número de nós da horizontais da malha
Nu = 1e-6;     % Viscosidade cinemática da água [m^2/s]     
CFL = 1;       % Para o implícito não precisa de CFL
pho = 1000;    % Massa específica da água [kg/m^3];
gradp = -1e1;  % Gradiente de Pressão [Pa] (Coloquei menos por causa do sentido do escoamento)

% Criação da malha na extensão do domínio

delta = h/(N-1);           % Distância de um nó a outro
dt = CFL*(delta^2)/Nu;     % Tamnho do passo de tempo baseado na condição de estabilidade
Nt = ceil(tf/dt);          % Número total de passos de tempo
steppersec = Nt/tf;        % Número de passos de tempo em 1 segundo da vida real
y = linspace (0,h,N);      % Vetor dos intervalos do espaço verticais
z = y;                     % Vetor dos intervalos do espaço horizontais

% Coeficientes da equação discretizada, que entrarão na matriz

c = (1/pho)* gradp;
alpha = (Nu*dt)/delta^2;   
Beta = 1 + ((4 * Nu * dt)/delta^2);

% Inicialização e Condições de contorno

u = ones(N,M,Nt); % Pré alocação da matriz de velocidades

u(:,:,1) = 0; % Condição de contorno: velocidade zero nas paredes
u(1,:,:) = 0; % Condição de contorno: velocidade zero nas paredes
u(:,1,:) = 0; % Condição de contorno: velocidade zero nas paredes
u(N,:,:) = 0; % Condição de contorno: velocidade zero nas paredes
u(:,M,:) = 0; % Condição de contorno: velocidade zero nas paredes

% Criando a matriz pentadiagonal espaçada (com uma coluna de zeros no meio)

Dim = N-2; % O zero secreto vai estar na posição Dim da diagonal
Dim2 = Dim^2; %Isso é a dimensão da matriz

% Valores das diagonais
principal = Beta * ones(Dim2, 1);    % Diagonal principal

superior = -alpha * ones(Dim2-1, 1);   % Diagonal acima
superior(Dim:Dim:end) = 0;
superior3 = -alpha * ones (Dim2-Dim,1);

inferior = -alpha * ones(Dim2-1, 1);   % Diagonal abaixo
inferior(Dim:Dim:end) = 0;
inferior3 = -alpha *ones(Dim2-Dim,1);

% Criando a matriz tridiagonal
A = diag(principal) + diag(superior, 1) + diag(superior3,Dim) + ...
    diag(inferior, -1) + diag(inferior3,-Dim);

B = ones(Dim2,1);
B = B*(-c * dt);

for n = 1:Nt
    x = A\B;
    xmat = reshape(x, Dim, Dim);
    
    u(2:N-1, 2:M-1, n+1) = xmat; % Inclui em u os resultados do meio x 
    
    B = x - c * dt; % Atualiza B para o proximo passo de tempo
end

% Após o loop, u contém a solução completa com condições de contorno

%% Plotagem para ver o resultado final

[Z, Y] = meshgrid(z, y);  % Use os vetores z e y completos que incluem fronteiras

% ANIMAÇÃO FLUIDA DO ESCOAMENTO 

% f = figure; % Uma janela única
% 
% for n = 1:Nt
%     if ~ishandle(f) % Use ishandle(f) para checar se a janela está aberta
%         disp('Animação encerrada pelo usuário.');
%         break
%     end %Encerra a simulação quando clicka na janela da figura para fechar
% 
% % Plot da animação do escoamento 
%     surf(Z, Y, u(:,:,n), 'EdgeColor', 'none');
%     view(45, 30);
%     colormap(turbo);
%     colorbar;
%     title('Evolução da velociade com tempo');
%     xlabel('Coordenada z (m)');
%     ylabel('Coordenada y (m)');
%     zlabel('Velocidade (m/s)');
%     zlim([0 max(u(:))]);
%     axis([0 L 0 h 0 1.1*max(u(:))]);
%     shading interp;
%     drawnow;     % Atualiza a tela para cada iteração
% end

% --- Derivada du/dy no último passo de tempo, perfil em z central ---

coluna_z = round(M/2);            % coluna central do domínio z
u_perfil = u(:, coluna_z, Nt);    % perfil de velocidades ao longo de y no instante final

du_dy = zeros(size(u_perfil));    % inicializa vetor das derivadas

% Paredes
du_dy(1) = (4*u_perfil(2)-3*u_perfil(1)-u_perfil(3))/(2*delta); % Diferença progressiva com 3 pontos
du_dy(N) = (u_perfil(N-2)-4*u_perfil(N-1)+3*u_perfil(N))/(2*delta); % Diferença regresiva com 3 pontos

for i = 2:N-1
    du_dy(i) = ((u_perfil(i+1) - u_perfil(i))/delta + (u_perfil(i) - u_perfil(i-1))/delta)/2;
end

% Plot do gráfico das derivadas
figure;
plot(y, du_dy, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Coordenada y (m)',FontSize=14);
ylabel('Valor da derivada no ponto',FontSize=14);
title('Derivada da velocidade (du/dy) ao longo de y no tempo final',FontSize=14);
grid on;

%% Plotagem alternativa para um só passo de tempo

% Geração da malha para plotagem
[X, Y] = meshgrid(z, y);

% Gráfico 3D da velocidade no tempo final
figure;
surf(X, Y, u(:, :, end));  % u(y,z) no tempo final
shading interp;
xlabel('z [m]');
ylabel('y [m]');
zlabel('u [m/s]');
title('Perfil de Velocidade 3D - Escoamento de Poiseuille');
colorbar;
colormap("turbo")
view(45, 30); % Ângulo de visualização

tempo_execucao = toc; % Calcula o tempo de execução
fprintf('Tempo real de execução: %.6f segundos\n', tempo_execucao);