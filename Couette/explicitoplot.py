import numpy as np
import matplotlib.pyplot as plt

# Carregar os dados do arquivo
data = np.loadtxt("couette_resultssexp.txt", skiprows=1)

y = data[:, 0]  # Primeira coluna: y
u_tempo1 = data[:, 1]  # Segunda coluna: u no tempo 1
u_tempo2 = data[:, 2]  # Terceira coluna: u no tempo 2
u_tempo3 = data[:, 3]  # Quarta coluna: u no tempo 3
u_tempo4 = data[:, 4]  # Quarta coluna: u no tempo 4
u_tempo5 = data[:, 5]  # Quarta coluna: u no tempo 5

# Plotar os perfis de velocidade
plt.plot(u_tempo1, y, label="t = 1 s") 
plt.plot(u_tempo2,y, label="t = 5 s") 
plt.plot(u_tempo3,y, label="t = 10 s") 
plt.plot(u_tempo4,y, label="t = 30 s") 
plt.plot(u_tempo5,y, label="t = 60 s") 
plt.ylabel("Distância y (m)")
plt.xlabel("Velocidade u (m/s)")
plt.title("Perfis de Velocidade do Escoamento de Couette")
plt.legend()
plt.grid()
plt.show()
