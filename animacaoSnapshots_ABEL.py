import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import imageio.v2 as imageio

pasta = os.path.dirname(os.path.abspath(__file__))

os.makedirs(os.path.join(pasta, "frames"), exist_ok=True)

snap_files = sorted(glob.glob(os.path.join(pasta, "snapshots", "snap_*.dat")))

if not snap_files:
    print("Nenhum snapshot encontrado em snapshots/")
    exit()

print(f"{len(snap_files)} snapshots encontrados.")

# Pré-calcula limite global da norma para escala consistente entre frames
print("Calculando limites globais...")
mag_max = -np.inf

for snap_file in snap_files:
    dados = np.loadtxt(snap_file, skiprows=1)
    mag = np.sqrt(dados[:, 2]**2 + dados[:, 3]**2)
    mag_max = max(mag_max, mag.max())

frame_paths = []

for idx, snap_file in enumerate(snap_files):
    with open(snap_file, 'r') as f:
        header = f.readline()
    tempo = float(header.split("=")[1].strip())

    dados = np.loadtxt(snap_file, skiprows=1)

    x = dados[:, 0]
    y = dados[:, 1]
    u = dados[:, 2]
    v = dados[:, 3]

    x_unicos = np.unique(x)
    y_unicos = np.unique(y)
    Nx = len(x_unicos)
    Ny = len(y_unicos)

    X = x.reshape(Ny, Nx)
    Y = y.reshape(Ny, Nx)
    U = u.reshape(Ny, Nx)
    V = v.reshape(Ny, Nx)
    Vel_Mag = np.sqrt(U**2 + V**2)

    # Subamostrar para streamplot (evita lentidão em malhas grandes)
    step = max(1, Nx // 50)
    Xs = X[::step, ::step]
    Ys = Y[::step, ::step]
    Us = U[::step, ::step]
    Vs = V[::step, ::step]

    fig, ax = plt.subplots(figsize=(7, 6))
    fig.suptitle(f"t = {tempo:.4f} s  —  frame {idx+1}/{len(snap_files)}", fontsize=13)

    cf = ax.contourf(X, Y, Vel_Mag, levels=100, cmap='turbo', vmin=0, vmax=mag_max)
    plt.colorbar(cf, ax=ax, label='|V| [m/s]')
    ax.streamplot(Xs.T, Ys.T, Us.T, Vs.T, color='w', linewidth=0.5, density=2.0, arrowsize=0.8)
    ax.set_title("Magnitude de Velocidade + Linhas de Corrente")
    ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    plt.tight_layout()

    frame_path = os.path.join(pasta, "frames", f"frame_{idx+1:05d}.png")
    plt.savefig(frame_path, dpi=80)
    plt.close(fig)
    frame_paths.append(frame_path)
    print(f"  [{idx+1}/{len(snap_files)}] {snap_file}  (t={tempo:.4f}s)")

print("Gerando GIF...")
frames = [imageio.imread(fp) for fp in frame_paths]
imageio.mimsave(os.path.join(pasta, "animacao_acpv.gif"), frames, fps=15, loop=0)
print(f"GIF salvo: animacao_acpv.gif  ({len(frames)} frames, 15 fps)")
