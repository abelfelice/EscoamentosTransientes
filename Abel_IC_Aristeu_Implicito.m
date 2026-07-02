% Solução Implicita do escoamento de Couette plano unidimensional ABEL

clc; clear all; close all;

%% Modelo computacional Implícito:

% Parâmetros de entrada
tic
h = 0.01;     % Distância entre as duas placas [m]
N = 100;      % Numero de nós da malha
Nu = 1e-6;    % Viscosidade cinemática da água [m^2/s]     
CFL = 0.1;      % Para o implícito não precisa de CFL
Mu = 1e-3;    % Viscosidade dinâmica da água
L = 1;        % Tamanho da placa [m]

% Condições de contorno e de valor inicial

V0 = 0;  %[m/s]
Vf = 1;  %[m/s]
t0 = 0;  %[s]
tf = 60; %[s]


% Criação da malha na extensão do domínio

dy = h/(N-1);          % Distância de um nó ao outro na malha
dt = CFL*(dy^2)/Nu;    % Tamnho do passo de tempo baseado na condição de estabilidade
Nt = ceil((tf-t0)/dt); % Número total de passos de tempo
y = linspace (0,h,N);  % Vetor dos intervalos do espaço
alpha = (Nu*dt)/dy^2;  % Coeficiente da matriz
steppersec = Nt/tf;    % Número de passos de tempo em 1 segundo da vida real

% Implementação do modelo matemático 

u = zeros(N,Nt);  % Pré inicialização da matriz velocidade
u(1,:) = V0;      % Na linha 1 todos elementos são 0
u(N,:) = Vf;      % Na linha N todos os elementos são 1
u(:,1) = t0;      % Coluna 1 todos elementos são igual a 0

% Criando a matriz tridiagonal

principal = -(1 + 2*alpha) * ones(N-2, 1);    % Diagonal principal
superior = alpha * ones(N-3, 1);              % Diagonal acima
inferior = alpha * ones(N-3, 1);              % Diagonal abaixo
A = diag(principal) + diag(superior, 1) + diag(inferior, -1); 
x = zeros(N-2,Nt);

% Impondo a condição de contorno e resolvendo o sistema linear

B = zeros(N-2,1);
for n = 1:Nt
B(end) = B(end) - alpha*Vf;

x(:,n) = A\B;
B = -x(:,n);
u(:,n) = [V0; x(:,n); Vf];
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
% Posicao N-1 ao longo dos diferentes tempos
evolve = u(N-1,:);

%Vetor com Nt valores 
aux = linspace(0,tf,Nt);

%% Plotar resultado

%Evolução dos perfis de velocidade:
figure('Name', 'Perfis de Velocidade','Position', [100 100 1200 500]);
subplot(1,2,1);
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

% Evolução das velocidades na posição N-1 ao longo do tempo
subplot(1,2,2);
plot(aux,evolve,LineWidth=1.5);
grid on; 
xlabel('Tempo total de simulação [s]');
ylabel('Velocidade u [m/s]');
title('Evolução da velocidade em um mesmo ponto ao longo dos tempos');
tempo_execucao = toc; % Calcula o tempo de execução
fprintf('Tempo real de execução: %.6f segundos\n', tempo_execucao);

% % Animação do perfil de velocidade ao longo do tempo
% fig = figure('Name','Animação do Perfil de Velocidade','Position',[100 100 600 500]);
% 
% ax  = axes('Parent',fig);
% hLn = plot(ax,u(:,1),y,'b','LineWidth',1.5);   % só cria 1 vez
% xlabel(ax,'Velocidade u [m/s]');
% ylabel(ax,'Distância y [m]');
% title(ax,'Perfil de Velocidade');
% axis(ax,[0 1 0 h]);
% grid(ax,'on');
% 
% pauseTime = 0.02;   % 0,15 s ≈ 6–7 fps (ajuste a gosto)
% 
% for n = 1:10:Nt
%     if ~isvalid(fig); break;
%     end
%     set(hLn,'XData',u(:,n));  % atualiza a linha
%     title(ax,sprintf('Perfil de Velocidade – Tempo = %.2f s',(n)*dt));
%     drawnow
%     pause(pauseTime);         % controla a velocidade
% end

Re = (Vf*h)/Nu;


%% Seção final - Potências: Discreta x Continua

% % ----- Potência mecânica Discreta (1 ordem) -----------
% 
% % Pré-alocação
% F   = zeros(1,Nt);   % Força viscosa na placa superior [N]
% Pot = zeros(1,Nt);   % Potência mecânica Discreta       [W]
% 
% 
% % Tensão na parede (y = h) por diferença regressiva:
% %   F(n)   = mu * (u_N - u_{N-1})/dy * L^2
% %   Pot(n) = F(n) * Vf            (Vf constante)
% for n = 1:Nt
%     F(n)   = Mu * ((u(N,n) - u(N-1,n))/dy) * L^2;
%     Pot(n) = F(n) * Vf;
% end


% ----- Potência mecânica Discreta (2 ordem) -----------

% Pré-alocação

F    = zeros(1,Nt);
Pot  = zeros(1,Nt);
potv = zeros(1,Nt);

for n = 1:Nt
    % --- Derivada de 2ª ordem ---
    dudy        = zeros(N,1);
    dudy(N)     = ( 3*u(N,n)  - 4*u(N-1,n) + u(N-2,n)) / (2*dy);  % regressiva  (y=h)

    % --- potência mecânica (tensão de parede em y=h) ---
    F(n)   = Mu * dudy(N) * L^2;
    Pot(n) = F(n) * Vf;

end


% ----- Potência mecânica Continua ----------

% Gradiente do perfil de Fourier avaliado em y = h:
%   tau_w(t) = (mu*Vf/h) * [ 1 + 2*sum_k exp(-k^2*pi^2*nu*t/h^2) ]
%   Pmec(t)  = tau_w(t) * Vf * L^2
t      = (1:Nt)*dt;     % tempo físico de cada coluna (coluna n -> t = n*dt)
Nmodes = 1000;          % nº de termos da série
serie  = ones(1,Nt);    % o termo "1 +"
for k = 1:Nmodes
    serie = serie + 2*exp( -(k^2)*(pi^2)*Nu*t / (h^2) );
end
Pot_an = (Mu*Vf^2*L^2/h) .* serie;

% Regime permanente (conferência): Pmec(inf) = mu*Vf^2*L^2/h
Pot_perm = Mu*Vf^2*L^2/h;
fprintf('Potência mecânica em regime permanente: %.4f W\n', Pot_perm);

% ----- Comparação -------------------------------------------------------
figure('Name','Potência mecânica: Discreta x Continua', ...
       'Position',[100 100 850 500]);
plot(t, Pot,    'b-',  'LineWidth', 1.5); hold on
plot(t, Pot_an, 'r--', 'LineWidth', 1.5);
yline(Pot_perm, 'k:',  'LineWidth', 1.2);
hold off
grid on
xlabel('Tempo [s]');
ylabel('Potência mecânica [W]');
title('Potência mecânica fornecida pela placa: Discreta x Continua');
legend('Discreta (diferenças finitas)','Continua (série de Fourier)', ...
       'Regime permanente','Location','northeast');
ylim([0 1]);   % foca na convergência; o pico em t->0 é singular e fica fora de escala


%% Potência viscosa: Discreta x Continua
% (usa t, Nmodes e Pot_perm definidos na seção da potência mecânica)

potv = zeros(1,Nt);   % Potência viscosa (Discreta) [W]

% ----- Potência viscosa Discreta ---------------------------------------
% Phi = mu*(du/dy)^2 integrado na altura, vezes a área L^2.
% Gradiente nas faces: (du/dy)_{i+1/2} = (u_{i+1}-u_i)/dy
% Integral pela regra do ponto médio: soma das contribuições de cada face.
for n = 1:Nt
    potv(n) = 0;                                    % zera o acumulador
    for i = 2:N-1
        dudy2(i)    = ((u(i+1,n) - u(i,n)) / dy).^2;         % gradiente na face i
        potv(n) = Mu * sum(dudy2) * dy * L^2; % ACUMULA cada face
    end
end


% ----- Potência viscosa Continua --------------------------------------
% Integrando (du/dy)^2 do perfil de Fourier (ortogonalidade dos cossenos):
%   integral_0^h (du/dy)^2 dy = (Vf^2/h)*[1 + 2*sum_k exp(-2*k^2*pi^2*nu*t/h^2)]
%   Pdis(t) = mu*L^2 * (integral acima)
% OBS: o expoente tem 2*k^2 -- a dissipação envolve o QUADRADO da série,
%      então cada modo decai com o dobro da taxa da potência mecânica.
serie_v = ones(1,Nt);
for k = 1:Nmodes
    serie_v = serie_v + 2*exp( -2*(k^2)*(pi^2)*Nu*t / (h^2) );
end
potv_an = (Mu*Vf^2*L^2/h) .* serie_v;

% ----- Comparação -------------------------------------------------------
figure('Name','Potência viscosa: Discreta x Continua', ...
       'Position',[100 100 850 500]);
plot(t, potv,    'b-',  'LineWidth', 1.5); hold on
plot(t, potv_an, 'r--', 'LineWidth', 1.5);
yline(Pot_perm, 'k:',  'LineWidth', 1.2);
hold off
grid on
xlabel('Tempo [s]');
ylabel('Potência viscosa transformada [W]');
title('Potência viscosa transformada no fluido: Discreta x Continua');
legend('Discreta (diferenças finitas)','Continua (série de Fourier)', ...
       'Regime permanente','Location','northeast');
ylim([0 1]);   % foco na convergência; pico em t->0 (singular) fora de escala




%% Overlay - Potência mecânica x viscosa (Continua e Discreta)
figure('Name','Mecânica x Viscosa','Position',[100 100 1200 500]);

% ----- Solução Continua -----
subplot(1,2,1);
fill([t fliplr(t)], [Pot_an fliplr(potv_an)], [0.85 0.85 0.6], ...
     'EdgeColor','none','FaceAlpha',0.35); hold on
plot(t, Pot_an,  'b-', 'LineWidth',1.5);
plot(t, potv_an, 'r-', 'LineWidth',1.5);
plot([t(1) t(end)], [Pot_perm Pot_perm], 'k:', 'LineWidth',1.2);
hold off; grid on
xlabel('Tempo [s]'); ylabel('Potência [W]');
title('Continua: mecânica x viscosa');
legend('dE_c/dt (energia cinética)','Mecânica (entrada)', ...
       'Viscosa (transformada)','Regime permanente','Location','northeast');
ylim([0 1]);

% ----- Solução Discreta -----
subplot(1,2,2);
fill([t fliplr(t)], [Pot fliplr(potv)], [0.85 0.85 0.6], ...
     'EdgeColor','none','FaceAlpha',0.35); hold on
plot(t, Pot,  'b-', 'LineWidth',1.5);
plot(t, potv, 'r-', 'LineWidth',1.5);
plot([t(1) t(end)], [Pot_perm Pot_perm], 'k:', 'LineWidth',1.2);
hold off; grid on
xlabel('Tempo [s]'); ylabel('Potência [W]');
title('Discreta: mecânica x viscosa');
legend('dE_c/dt (energia cinética)','Mecânica (entrada)', ...
       'Viscosa (transformada)','Regime permanente','Location','northeast');
ylim([0 1]);

%% Diferença de potências:  Pmec - Pvis  =  dEc/dt
dPot    = Pot    - potv;      % Discreta
dPot_an = Pot_an - potv_an;   % Continua

figure('Name','Pmec - Pvis (taxa de energia cinética)','Position',[100 100 850 500]);
plot(t, dPot,    'b-',  'LineWidth',1.5); hold on
plot(t, dPot_an, 'r--', 'LineWidth',1.5);
yline(0, 'k:', 'LineWidth',1.0);
hold off
grid on
xlabel('Tempo [s]');
ylabel('P_{mec} - P_{vis}  =  dE_c/dt   [W]');
title('Diferença entre potência mecânica e viscosa');
legend('Discreta','Continua','Location','northeast');
ylim([0 0.5]);   % foco na convergência; pico singular em t->0 fica fora de escala


%% Potências contínuas: mecânica x viscosa — área preenchida + escalas ajustadas
% Usa Pot_an, potv_an, Pot_perm e t (já calculados).

figure('Name','Potências contínuas: mecânica x viscosa','Position',[100 100 950 560]);

% --- área entre as curvas (= dEc/dt) ---
fill([t fliplr(t)], [Pot_an fliplr(potv_an)], [0.95 0.85 0.45], ...
     'EdgeColor','none','FaceAlpha',0.5); hold on

% --- curvas ---
plot(t, Pot_an,  'b-', 'LineWidth', 2.0);
plot(t, potv_an, 'r-', 'LineWidth', 2.0);
yline(Pot_perm, 'k:', 'LineWidth', 1.2);
hold off
grid on

xlabel('Tempo [s]', 'FontSize', 16);
ylabel('Potência [W]', 'FontSize', 16);
title('Potência mecânica e viscosa (Solução Contínua)', 'FontSize', 18);
legend('E_c(t) (área entre curvas)','Potência Mecânica ', ...
       ' Potência Viscosa ','Regime permanente', ...
       'Location','northeast', 'FontSize', 13);
set(gca, 'FontSize', 14);

% --- escalas ajustadas (aplicar por último, senão o fill reescala) ---
xlim([0 20]);          % 0 até o tempo final da simulação (60 s)
ylim([0 20]);           % foco na convergência ao platô 0,1 W
yticks(0:5:20);

%% Energia cinética do escoamento
rho = Mu / Nu;                 % densidade [kg/m^3]  (= 1000 para a água)

%--- numérica: integral de (1/2) rho u^2 na altura, vezes L^2 ---
Ec = zeros(1,Nt);
for n = 1:Nt
    Ec(n) = 0.5 * rho * L^2 * trapz(y(:), u(:,n).^2);
end

% Ec = zeros(1,Nt);
% for n = 1:Nt
%     Ec(n) = 0.5 * rho * L^2 * dy * sum( u(2:N-1,n).^2 ); Não deu certo
% end

% --- analítica (mesma série de Fourier) ---
serie1 = zeros(1,Nt);          % sum beta_n / n^2
serie2 = zeros(1,Nt);          % sum beta_n^2 / n^2
for k = 1:Nmodes
    bk     = exp( -(k^2)*(pi^2)*Nu*t / (h^2) );   % beta_k
    serie1 = serie1 + bk    / k^2;
    serie2 = serie2 + bk.^2 / k^2;
end
Ec_an = (rho*Vf^2*L^2*h/2) .* ( 1/3 - (4/pi^2)*serie1 + (2/pi^2)*serie2 );


%% Energia cinética em função do tempo: numérica x analítica
Ec_perm = rho*Vf^2*L^2*h/6;        % valor de regime permanente [J]

figure('Name','Energia cinética','Position',[100 100 850 500]);
plot(t, Ec,    'b-',  'LineWidth',1.5); hold on
plot(t, Ec_an, 'r--', 'LineWidth',1.5);
yline(Ec_perm, 'k:',  'LineWidth',1.2);
hold off
grid on
xlabel('Tempo [s]');
ylabel('Energia cinética E_c [J]');
title('Energia cinética do escoamento: Discreta x Contínua');
legend('Discreta','Contínua','Regime permanente','Location','southeast');
ylim([0 1.1*Ec_perm]);






%% ================================================================
%  Energia cinética analítica — duas rotas independentes
%   Rota 1 (campo):     E_c = (1/2) rho L^2 ∫_0^h u^2 dy
%   Rota 2 (potências): E_c = ∫_0^t [ P_mec(τ) - P_vis(τ) ] dτ
%  Devem coincidir (é o balanço de energia em forma analítica).
%  Reaproveita h, N, Nu, Mu, L, Vf, dt, Nt do seu script.
% ================================================================
rho    = Mu/Nu;            % densidade [kg/m^3]  (viscosidade dinâmica / cinemática)
Nmodes = 1000;             % nº de modos das séries
t      = (1:Nt)*dt;        % vetor de tempo (coluna n  ->  t = n*dt)

% ---------- Rota 1: direto do campo de velocidades ----------
% Prefator com rho (inércia); viscosidade entra só via beta_k (= nu no expoente).
%   E_c = (rho Vf^2 L^2 h /2) [ 1/3 - (4/pi^2) Σ β_k/k^2 + (2/pi^2) Σ β_k^2/k^2 ]
S1 = zeros(1,Nt);   S2 = zeros(1,Nt);
for k = 1:Nmodes
    bk = exp( -(k^2)*(pi^2)*Nu.*t / (h^2) );     % β_k(t)
    S1 = S1 + bk    / k^2;                       % Σ β_k / k^2
    S2 = S2 + bk.^2 / k^2;                       % Σ β_k^2 / k^2
end
Ec_campo = (rho*Vf^2*L^2*h/2) .* ( 1/3 - (4/pi^2)*S1 + (2/pi^2)*S2 );

% ---------- Rota 2: a partir das potências ----------
% Potências analíticas carregam mu (viscosidade DINÂMICA -> tensão = mu*du/dy):
%   P_mec = (mu Vf^2 L^2/h)[1 + 2 Σ β_k]
%   P_vis = (mu Vf^2 L^2/h)[1 + 2 Σ β_k^2]
% Integrando a diferença no tempo (E_c(0)=0), com a_k = k^2 pi^2 nu / h^2:
%   ∫β_k dτ = (1-β_k)/a_k ,   ∫β_k^2 dτ = (1-β_k^2)/(2 a_k)
%   E_c = (2 mu Vf^2 L^2/h) Σ [ (1-β_k)/a_k - (1-β_k^2)/(2 a_k) ]
% O fator (2 mu/h)·(1/a_k) traz h^2/nu, e mu/nu = rho -> reaparece a inércia.
Pmec_an = ones(1,Nt);   Pvis_an = ones(1,Nt);    % os "1" das séries
Sbeta   = zeros(1,Nt);                            % parte dependente de β da integral
for k = 1:Nmodes
    bk = exp( -(k^2)*(pi^2)*Nu.*t / (h^2) );
    ak = (k^2)*(pi^2)*Nu / (h^2);                 % taxa de decaimento do modo k
    Pmec_an = Pmec_an + 2*bk;                     % monta Σ β_k
    Pvis_an = Pvis_an + 2*bk.^2;                  % monta Σ β_k^2
    Sbeta   = Sbeta - bk./ak + (bk.^2)./(2*ak);   % -Σ β_k/a_k + Σ β_k^2/(2 a_k)
end
Pmec_an = (Mu*Vf^2*L^2/h) .* Pmec_an;             % potência mecânica analítica [W]
Pvis_an = (Mu*Vf^2*L^2/h) .* Pvis_an;             % potência viscosa  analítica [W]
const   = h^2/(12*Nu);                            % Σ 1/(2 a_k) em forma fechada (ζ(2)=π²/6)
Ec_potencias = (2*Mu*Vf^2*L^2/h) .* ( const + Sbeta );

% ---------- Comparação ----------
Ec_perm = rho*Vf^2*L^2*h/6;                       % valor de regime permanente
figure('Name','Energia cinética analítica: duas rotas','Position',[100 100 900 560]);

%subplot(3,1,[1 2]);
plot(t, Ec_campo,     'b-',  'LineWidth', 2.0); hold on
plot(t, Ec_potencias, 'r--', 'LineWidth', 1.6);
yline(Ec_perm, 'k:', 'LineWidth', 1.2);
hold off; grid on
xlabel('Tempo [s]'); ylabel('E_c  [J]');
title('Energia cinética: direto do campo  \times  a partir das potências');
legend('Rota 1:  \rho L^2 \int u^2 dy / 2', ...
       'Rota 2:  \int (P_{mec}-P_{vis}) d\tau', ...
       'Regime permanente', 'Location','southeast');


fprintf('E_c final  (rota 1, campo)      = %.6f J\n', Ec_campo(end));
fprintf('E_c final  (rota 2, potências)  = %.6f J\n', Ec_potencias(end));
fprintf('E_c permanente                  = %.6f J\n', Ec_perm);
fprintf('Diferença máxima entre as rotas = %.3e J\n', max(abs(Ec_campo - Ec_potencias)));







































 

