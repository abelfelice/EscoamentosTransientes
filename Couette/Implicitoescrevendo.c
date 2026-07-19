#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Função para resolver o sistema tridiagonal (Método de Thomas)
void solveTridiagonal(double *a, double *b, double *c, double *d, double *x, int n) {
    double *c_prime = (double *)malloc(n * sizeof(double));
    double *d_prime = (double *)malloc(n * sizeof(double));

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / m;
    }

    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    free(c_prime);
    free(d_prime);
}

int main() {

     // Início da medição do tempo de execução 
    clock_t t_start = clock();

    // Parâmetros de entrada
    double h = 0.01;   // Distância entre as placas [m]
    int N = 100;        // Número de nós da malha
    double Nu = 1e-6;  // Viscosidade cinemática [m^2/s]
    double CFL = 1;    // Implicito n precisa de CFL

    // Condições de contorno
    double V0 = 0.0;  // Velocidade da placa inferior [m/s]
    double Vf = 1.0;  // Velocidade da placa superior [m/s]
    double tf = 60.0; // Tempo final [s]

    // Criação da malha no domínio
    double dy = h / (N - 1);
    double dt = CFL * (dy * dy) / Nu; // Passo de tempo
    int Nt = (int)(round(tf / dt)); // Número total de passos de tempo
    double alpha = Nu * dt / (dy * dy);

    // Vetores para y, u e auxiliares
    double y[N];
    // Alocação dinâmica da matriz u (evita estouro de pilha para N ou tf grandes)
    double **u = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        u[i] = (double *)malloc(Nt * sizeof(double));
    }
    double a[N - 2], b[N - 2], c[N - 2], d[N - 2], x[N - 2];

    // Inicializando o vetor y
    for (int i = 0; i < N; i++) {
        y[i] = i * dy;
    }

    // Inicializando a matriz u com zeros
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Nt; j++) {
            u[i][j] = 0.0;
        }
    }

    // Condições de contorno
    for (int n = 0; n < Nt; n++) {
        u[0][n] = V0;      // Velocidade da placa inferior
        u[N - 1][n] = Vf;  // Velocidade da placa superior
    }

    // Configuração da matriz tridiagonal
    for (int i = 0; i < N - 2; i++) {
        a[i] = alpha;            // Elementos da diagonal inferior
        b[i] = -(1 + 2 * alpha); // Elementos da diagonal principal
        c[i] = alpha;            // Elementos da diagonal superior
    }

    // Implementação do método implícito
    for (int n = 0; n < Nt - 1; n++) {
        // Configurando o vetor d para o sistema linear
        for (int i = 0; i < N - 2; i++) {
            d[i] = -u[i + 1][n];
        }

        // Ajustando as condições de contorno no vetor d
        d[0] -= alpha * u[0][n + 1];        // Condição de contorno inferior
        d[N - 3] -= alpha * u[N - 1][n + 1]; // Condição de contorno superior

        // Resolvendo o sistema tridiagonal
        solveTridiagonal(a, b, c, d, x, N - 2);

        // Atualizando os valores da matriz u
        for (int i = 1; i < N - 1; i++) {
            u[i][n + 1] = x[i - 1];
        }
    }

    // Fim da medição
    clock_t t_end = clock();
    double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    printf("Tempo de execucao do laco numerico: %.4f s\n", elapsed);

    // Índices correspondentes a t = 38s, 40s e 55s (batem com os rótulos do script Python)
    int idx1 = (int)(1.0 / dt);
    int idx2 = (int)(5.0 / dt);
    int idx3 = (int)(10.0 / dt);
    int idx4 = (int)(30.0 / dt);
    int idx5 = (int)(60.0 / dt);

    // Proteção contra índices inválidos (caso tf ou dt mudem)
    if (idx1 < 0) idx1 = 0; if (idx1 >= Nt) idx1 = Nt - 1;
    if (idx2 < 0) idx2 = 0; if (idx2 >= Nt) idx2 = Nt - 1;
    if (idx3 < 0) idx3 = 0; if (idx3 >= Nt) idx3 = Nt - 1;
    if (idx4 < 0) idx4 = 0; if (idx4 >= Nt) idx4 = Nt - 1;
    if (idx5 < 0) idx5 = 0; if (idx5 >= Nt) idx5 = Nt - 1;

    // Impressão do perfil de velocidade em t = 55s
    printf("Perfil de velocidade em t = %.2f s:\n", idx5 * dt);
    for (int i = 0; i < N; i++) {
        printf("y = %.5f, u = %.5f\n", y[i], u[i][idx5]);
    }

    // Salvar os resultados no arquivo
    FILE *file = fopen("couette_resultssimp.txt", "w");
    if (file == NULL) {
        printf("Erro ao criar o arquivo.\n");
        for (int i = 0; i < N; i++) free(u[i]);
        free(u);
        return 1;
    }

    fprintf(file, "y u_tempo1 u_tempo2 u_tempo3 u_tempo4 utempo5\n");
    for (int i = 0; i < N; i++) {
        fprintf(file, "%.5f %.5f %.5f %.5f %.5f %.5f\n", y[i], u[i][idx1], u[i][idx2], u[i][idx3], u[i][idx4], u[i][idx5]);
    }

    fclose(file);
    printf("Resultados salvos em couette_resultsimp.txt\n");

    






    // Liberação da memória alocada
    for (int i = 0; i < N; i++) free(u[i]);

    free(u);
    return 0;
}
