#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main() {

    // Início da medição do tempo de execução 
    clock_t t_start = clock();

    // Parâmetros de entrada
    double h = 0.01;  // Distância entre as duas placas [m]
    int N = 100;      // Número de nós da malha
    double Nu = 1e-6; // Viscosidade cinemática [m^2/s]
    double CFL = 0.5;

    // Condições de contorno
    double V0 = 0.0;  // Velocidade da placa inferior [m/s]
    double Vf = 1.0;  // Velocidade da placa superior [m/s]
    double t0 = 0.0;  // Tempo inicial [s]
    double tf = 60.0; // Tempo final [s]

    // Criação da malha no domínio
    double dy = h / (N - 1);
    double dt = CFL * (dy * dy) / Nu;  // Passo de tempo
    int Nt = (int)round((tf - t0) / dt) + 1;  // Número total de passos de tempo

    double y[N];
    // Alocação dinâmica da matriz u (evita estouro de pilha para N ou tf grandes)
    double **u = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        u[i] = (double *)malloc(Nt * sizeof(double));
    }

    // Limpando a memória do vetor y
    for (int i = 0; i < N; i++) {
        y[i] = 0;
    }

    // Limpando a memória da matriz u
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Nt; j++) {
            u[i][j] = 0;
        }
    }

    // Linspace do vetor y
    for (int i = 0; i < N; i++) {
        y[i] = i * dy;
    }

    // Condições de contorno
    for (int n = 0; n < Nt; n++) {
        u[0][n] = V0;     // Velocidade na placa inferior
        u[N - 1][n] = Vf; // Velocidade na placa superior
    }



    // Implementação do método explícito
    for (int n = 0; n < Nt - 1; n++) {
        for (int i = 1; i < N - 1; i++) { // Posição 1 até N-2 (interno)
            u[i][n + 1] = u[i][n] + Nu * (dt / (dy * dy)) * (u[i + 1][n] - 2 * u[i][n] + u[i - 1][n]);
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

    // Salvar os resultados no arquivo
    FILE *file = fopen("couette_resultssexp.txt", "w");
    if (file == NULL) {
        printf("Erro ao criar o arquivo.\n");
        for (int i = 0; i < N; i++) free(u[i]);
        free(u);
        return 1;
    }

    fprintf(file, "y u_tempo1 u_tempo2 u_tempo3 u_tempo4 u_tempo5\n");
    for (int i = 0; i < N; i++) {
        fprintf(file, "%.5f %.5f %.5f %.5f %.5f %.5f\n", y[i], u[i][idx1], u[i][idx2], u[i][idx3], u[i][idx4], u[i][idx5]);
    }

    fclose(file);
    printf("Resultados salvos em couette_resultssexp.txt\n");

    // Liberação da memória alocada
    for (int i = 0; i < N; i++) free(u[i]);
    free(u);

    return 0;
}
