clc; clear;
syms beta alpha c delta_t;  % Define as variáveis simbólicas
N = 5;
M = N;
Dim = N-2; %o zero secreto vai estar na posição Dim da diagonal
Dim2 = Dim^2;

%Obs: A segunda diagonal secundária estará Dim-1 posições deslocada em relação a
%primeira diagonal secundária

%Obs2: A segunda diagonal secundária é maciça

% Define a matriz U simbólica
U = sym('U', [N, M]);  % Matriz 5x5 para incluir os índices de contorno

% Define a matriz U^{n-1} simbólica
U_anterior = sym('U_anterior', [N, M]);

% Cria o vetor de equações simbólicas
equacoes = sym(zeros(Dim2, 1));  % Vetor para armazenar as equações

% Índice para o vetor de equações
indice = 1;

for j = 2:N-1
    for k = 2:M-1
        equacoes(indice) = beta * U(j, k) - alpha * (U(j+1, k) + U(j-1, k) + U(j, k+1) + U(j, k-1)) == U_anterior(j, k) - c * delta_t;
        indice = indice + 1;  % Incrementa o índice
    end
end

% Exibe o vetor de equações
disp('Vetor de equações simbólicas:');
disp(equacoes);

% Cria um vetor de todas as variáveis U(j, k)
variaveis_U = [];
for j = 2:N-1
    for k = 2:M-1
        variaveis_U = [variaveis_U, U(j, k)];
    end
end

% Converte as equações para a matriz A
[A] = equationsToMatrix(equacoes, variaveis_U);
disp('Matriz A:');
disp(A);
