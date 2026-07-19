import numpy as np
import matplotlib.pyplot as plt
import os

pasta = os.path.dirname(os.path.abspath(__file__))

dados = np.loadtxt(os.path.join(pasta, "resultadosacpv.dat"), skiprows=3)

x = dados[:, 0]
y = dados[:, 1]
u = dados[:, 2]
v = dados[:, 3]
p = dados[:, 4]

x_unicos = np.unique(x)
y_unicos = np.unique(y)

Nx = len(x_unicos)
Ny = len(y_unicos)

X = x.reshape(Ny,Nx)
Y = y.reshape(Ny,Nx)
U = u.reshape(Ny,Nx)
V = v.reshape(Ny,Nx)
P = p.reshape(Ny,Nx)



print("Umin/Umax =", U.min(), U.max())
print("Vmin/Vmax =", V.min(), V.max())

print(f"Malha: {Nx} x {Ny}")

# Plotagem de U

plt.figure(figsize=(8, 6))
cf = plt.contourf(X, Y, U, cmap='jet',levels=200, vmin=U.min(), vmax=U.max())
plt.title("Campo de velocidades Horizontais")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.colorbar(cf, label="u [m/s]")
plt.show()
plt.close()

# Plotagem de V

plt.figure(figsize=(8, 6))
cf = plt.contourf(X, Y, V, cmap='jet',levels=200, vmin=V.min(), vmax=V.max())
plt.title("Campo de velocidades Verticais")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.colorbar(cf, label="v [m/s]")
plt.show()
plt.close()

# Plotagem de P

plt.figure(figsize=(8, 6))
cf = plt.contourf(X, Y, P, cmap='jet',levels=200)
plt.title("Campo de Pressao")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.colorbar(cf, label="P [Pa]")
plt.show()
plt.close()

# Plotagem do campo misto

Vel_Mag = np.sqrt(U**2 + V**2)

plt.figure(figsize=(8, 6))
cf = plt.contourf(X, Y, Vel_Mag, levels=200, cmap='turbo')
cb = plt.colorbar(cf)
cb.set_label('Magnitude de Velocidade [m/s]')
plt.streamplot(X.T, Y.T, U.T, V.T, color='k', linewidth=0.7, density=2.5, arrowsize=1.0)
plt.title(f"Magnitude de Velocidade e Linhas de Corrente\n(Malha {Nx}x{Ny})")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.show()
plt.close()


# Plotagem das velocidades Ghia x Abel 

Umeioabel = np.loadtxt(os.path.join(pasta, "resultadosmeiosU.dat"))
Vmeioabel = np.loadtxt(os.path.join(pasta, "resultadosmeiosV.dat"))

Tamanho_y_abel = len(Vmeioabel)
y_abel = np.linspace(0.0, 1.0, Tamanho_y_abel)

Tamanho_x_abel = len(Umeioabel)
x_abel = np.linspace(0.0, 1.0, Tamanho_x_abel)


# Coordenada Y (Ordem Crescente)
y_ghia = np.array([
    0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000,
    0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000
])

# Velocidades u (Re=100)
ghia_u_dados_100 =  np.array([
        0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581,
        -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.00000])


# Coordenada X (Ordem Crescente)
x_ghia = np.array([
    0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000,
    0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000
])

# Velocidades v (Re=100)
ghia_v_dados_100 = np.array([

    0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454,
    -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.00000
    ])

# Velocidades u (Re=400)
ghia_u_dados_400 = np.array([
    0.00000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726, -0.17119, -0.11477,
    0.02135, 0.16256, 0.29093, 0.55892, 0.61756, 0.68439, 0.75837, 1.00000
])

# Velocidades u (Re=1000)
ghia_u_dados_1000 = np.array([
    0.00000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805, -0.10648, -0.06080,
    0.05702, 0.18719, 0.33304, 0.46604, 0.51117, 0.57492, 0.65928, 1.00000
])

# Velocidades v (Re=400)
ghia_v_dados_400 = np.array([
    0.00000, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186,
    -0.38598, -0.44993, -0.23827, -0.22847, -0.19254, -0.15663, -0.12146, 0.00000
])

# Velocidades v (Re=1000)
ghia_v_dados_1000 = np.array([
    0.00000, 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075, 0.32235, 0.02526,
    -0.31966, -0.42665, -0.51550, -0.39188, -0.33714, -0.27669, -0.21388, 0.00000
])


# Velocidades u (Re=10000) 
ghia_u_dados_10000 = np.array([
    0.00000, -0.33093, -0.33066, -0.32709, -0.30045, -0.26084, -0.28857, -0.05437, -0.03039,
    -0.03816, 0.05186, 0.12317, 0.46036, 0.46120, 0.47244, 0.47221, 1.00000
])


# Velocidades v (Re=10000)
ghia_v_dados_10000 = np.array([
    0.00000, 0.34707, 0.38038, 0.40879, 0.44976, 0.40742, 0.26442, 0.24523, 0.03418,
    -0.28842, -0.40398, -0.51868, -0.52839, -0.45785, -0.37180, -0.28017, 0.00000
])

# Meus dados 

plt.plot(y_ghia,ghia_u_dados_10000, "o", label = "Ghia")
plt.plot(x_abel,Umeioabel,label="Abel")
plt.title("Velocidades centrais (u) comparadas")
plt.xlabel("Y [m]")
plt.ylabel("u [m/s]")
plt.legend()
plt.show()
plt.close()
plt.plot(x_ghia,ghia_v_dados_10000, "o", label = "Ghia")
plt.plot(y_abel,Vmeioabel, label = "Abel")
plt.title("Velocidades centrais (v) comparadas")
plt.xlabel("X [m]")
plt.ylabel("v [m/s]")
plt.legend()
plt.show()
plt.close()



