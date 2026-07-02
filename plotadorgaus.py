import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def plotar_simulacao(caminho_arquivo='poiseuille_fortran_final.dat'):
    """
    Lê e plota os dados do último passo de tempo de uma simulação.

        caminho_arquivo (str): O nome do arquivo de dados a ser lido.
    """
    print(f"Lendo dados do passo de tempo final do arquivo: {caminho_arquivo}")
    
    try:
        # Carrega os dados do arquivo. Como o arquivo só contém o estado final (t=Nt),
        # a matriz 'dados' representará o último passo de tempo.
        dados = np.loadtxt(caminho_arquivo)
    except FileNotFoundError:
        print(f"ERRO: Arquivo '{caminho_arquivo}' não encontrado.")
        print("Certifique-se de que o script Python está no mesmo diretório que o arquivo de dados, ou forneça o caminho correto.")
        sys.exit(1)
    except Exception as e:
        print(f"Ocorreu um erro ao ler o arquivo: {e}")
        sys.exit(1)

    print("Arquivo lido com sucesso. Processando dados do estado final...")

    # Extrai as colunas de coordenadas e velocidade
    y_flat = dados[:, 0]
    z_flat = dados[:, 1]
    u_flat = dados[:, 2]

    # Determina as dimensões da malha (N x M)
    N = len(np.unique(y_flat))
    M = len(np.unique(z_flat))

    # Remodela os vetores para matrizes 2D para a plotagem
    Y = y_flat.reshape(M, N)
    Z = z_flat.reshape(M, N)
    U = u_flat.reshape(M, N)
    
    print(f"Dados do último passo de tempo remodelados para uma grade de {N}x{M}.")

    # --- Criação do Gráfico 3D ---
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(Y, Z, U, cmap='turbo', edgecolor='none')

    ax.set_title('Campo de Velocidades Final - Simulação Fortran')
    ax.set_xlabel('Coordenada Y (m)')
    ax.set_ylabel('Coordenada Z (m)')
    ax.set_zlabel('Velocidade U (m/s)')
    fig.colorbar(surf, shrink=0.6, aspect=10, label='Velocidade (m/s)')
    ax.autoscale(tight=True)
    ax.view_init(elev=30, azim=45) 

    print("Gráfico do último passo de tempo gerado. Exibindo...")
    plt.show()


if __name__ == "__main__":
    plotar_simulacao()