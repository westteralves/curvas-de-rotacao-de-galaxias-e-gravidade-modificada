import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define semente "fixa" para garantir a reprodutibilidade dos resultados
np.random.seed(42)


# Bojo central

num_pontos_bojo = 4000
# Distribuição exponencial para simular o decaimento da densidade estelar
raio_bojo = np.random.exponential(scale=0.7, size=num_pontos_bojo)
phi_bojo = np.random.uniform(0, 2 * np.pi, num_pontos_bojo)
cos_theta_bojo = np.random.uniform(-1, 1, num_pontos_bojo)
theta_bojo = np.arccos(cos_theta_bojo)

# Coordenadas com achatamento vertical (bojo oblato)
x_bojo = raio_bojo * np.sin(theta_bojo) * np.cos(phi_bojo) * 0.9
y_bojo = raio_bojo * np.sin(theta_bojo) * np.sin(phi_bojo) * 0.9
z_bojo = raio_bojo * np.cos(theta_bojo) * 0.45

# 2. Braços espirais

num_pontos_disco = 30000
num_braços = 4
x_braco, y_braco, z_braco, raio_braco = [], [], [], []

for i in range(num_braços):
    pontos_por_braco = num_pontos_disco // num_braços
    angulo_offset = i * (np.pi / 2)          # defasagem de 90 graus entre os braços
    r = np.linspace(0.6, 9.8, pontos_por_braco) + np.random.normal(0, 0.25, pontos_por_braco)
    r = np.clip(r, 0.4, 11.0)
    # Perfil logarítmico para produzir o enrolamento característico
    theta = 2.5 * np.log(r + 0.2) + angulo_offset + np.random.normal(0, 0.18, pontos_por_braco)
    x_braco.extend(r * np.cos(theta))
    y_braco.extend(r * np.sin(theta))
    # Espessura vertical cresce suavemente com o raio
    z_braco.extend(np.random.normal(0, 0.08 + 0.02 * r, pontos_por_braco))
    raio_braco.extend(r)

x_braco = np.array(x_braco)
y_braco = np.array(y_braco)
z_braco = np.array(z_braco)
raio_braco = np.array(raio_braco)


# Halo de matéria escura

# Malha para a cúpula transparente que delimita o halo
u_cupula = np.linspace(0, 2 * np.pi, 80)
v_cupula = np.linspace(0, np.pi, 80)
raio_halo = 13.5
x_cupula = raio_halo * np.outer(np.cos(u_cupula), np.sin(v_cupula))
y_cupula = raio_halo * np.outer(np.sin(u_cupula), np.sin(v_cupula))
z_cupula = raio_halo * 0.8 * np.outer(np.ones(np.size(u_cupula)), np.cos(v_cupula))

# Partículas do halo com distribuição de densidade em lei de potência
num_particulas_halo = 4000
r_halo_pts = np.random.power(a=1.8, size=num_particulas_halo) * raio_halo
phi_halo = np.random.uniform(0, 2 * np.pi, num_particulas_halo)
cos_theta_halo = np.random.uniform(-1, 1, num_particulas_halo)
theta_halo = np.arccos(cos_theta_halo)

x_halo = r_halo_pts * np.sin(theta_halo) * np.cos(phi_halo)
y_halo = r_halo_pts * np.sin(theta_halo) * np.sin(phi_halo)
z_halo = r_halo_pts * 0.8 * np.cos(theta_halo)


# Setas indicativas da dispersão de velocidades

num_latitudes = 10
num_longitudes = 16
raio_setas = raio_halo * 0.85

pos_x, pos_y, pos_z = [], [], []
dir_x, dir_y, dir_z = [], [], []

for lat in range(1, num_latitudes):
    theta_lat = lat * np.pi / (num_latitudes + 1)
    for lon in range(num_longitudes):
        phi_lon = lon * 2 * np.pi / num_longitudes
        x = raio_setas * np.sin(theta_lat) * np.cos(phi_lon)
        y = raio_setas * np.sin(theta_lat) * np.sin(phi_lon)
        z = raio_setas * 0.8 * np.cos(theta_lat)
        proj_xy = np.sqrt(x**2 + y**2)
        if proj_xy > 0.1:
            # Direção tangencial à rotação do halo
            tx = -y / proj_xy
            ty = x / proj_xy
            tz = 0
            norma = np.sqrt(tx**2 + ty**2 + tz**2)
            tx, ty, tz = tx / norma, ty / norma, tz / norma
            # Módulo da velocidade decai com a distância ao centro
            dist_centro = np.sqrt(x**2 + y**2 + (z / 0.8)**2)
            dist_norm = dist_centro / raio_halo
            modulo_vel = 1.2 * (1 - np.exp(-2.5 * dist_norm))
            dispersao = np.random.normal(0, 0.08, 1)[0]
            modulo_vel = max(0.15, modulo_vel + dispersao)
            tx *= modulo_vel
            ty *= modulo_vel
            tz *= modulo_vel
        else:
            tx, ty, tz = 0, 0, 0
        pos_x.append(x)
        pos_y.append(y)
        pos_z.append(z)
        dir_x.append(tx)
        dir_y.append(ty)
        dir_z.append(tz)


# Figura

fig = plt.figure(figsize=(12, 10), facecolor='#06060c', dpi=300)
ax = fig.add_subplot(111, projection='3d', facecolor='#06060c')

# Halo com partículas difusas e superfície quase transparente
ax.scatter(x_halo, y_halo, z_halo, c='#4a7dbd', s=1.8, alpha=0.12, edgecolors='none')
ax.plot_surface(x_cupula, y_cupula, z_cupula, color='#4b7bec', alpha=0.05,
                rstride=2, cstride=2, linewidth=0, edgecolor='none', shade=False)

# Vetores da dispersão de velocidades
ax.quiver(pos_x, pos_y, pos_z, dir_x, dir_y, dir_z,
          length=0.8, normalize=False, color='#88ccff', linewidth=1.2,
          arrow_length_ratio=0.35, alpha=0.75)

# Disco com braços espirais com gradiente de cor e tamanho variável
tamanhos = 1.8 + 3.0 * (raio_braco / raio_braco.max())
ax.scatter(x_braco, y_braco, z_braco, c=raio_braco, cmap='cool', s=tamanhos,
           alpha=0.75, edgecolors='none')

# Bojo com pontos mais brilhantes e concentrados
ax.scatter(x_bojo, y_bojo, z_bojo, c='#fed330', s=10.0, alpha=0.9,
           edgecolors='#fa8231', linewidth=0.3)

# Anéis concêntricos de referência no plano do disco
angulo_anel = np.linspace(0, 2 * np.pi, 300)
for raio_anel in [3.0, 5.5, 8.0]:
    ax.plot(raio_anel * np.cos(angulo_anel), raio_anel * np.sin(angulo_anel),
            np.zeros_like(angulo_anel), color='#fed330', linestyle='--',
            linewidth=0.6, alpha=0.4)

# Centro da galáxia
ax.scatter([0], [0], [0], c='white', s=40, alpha=0.95)

# Ajustes dos limites e orientação da cena
limite = raio_halo * 1.02
ax.set_xlim([-limite, limite])
ax.set_ylim([-limite, limite])
ax.set_zlim([-limite * 0.8, limite * 0.8])
ax.view_init(elev=26, azim=55)
ax.axis('off')

# Legenda
elementos_legenda = [
    Line2D([0], [0], marker='o', color='w', label='Halo de Matéria Escura (CDM)',
            markerfacecolor='#4a7dbd', markersize=8, alpha=0.5),
    Line2D([0], [0], marker='o', color='w', label='Disco e Bojo Bariônico',
            markerfacecolor='#fed330', markersize=8),
    Line2D([0], [0], color='#88ccff', lw=2, label='Dispersão de Velocidades (setas)'),
]
ax.legend(handles=elementos_legenda, loc='lower left', frameon=True,
          facecolor='#0a0a1a', edgecolor='#4b7bec', fontsize=12, labelcolor='white')

plt.title('Galáxia Espiral com Halo de Matéria Escura e Dispersão de Velocidades',
          color='white', fontsize=16, pad=20, fontweight='medium')

# Figuras em png e pdf
plt.savefig('galaxia_halo_dispersao.png', dpi=600,
            bbox_inches='tight', facecolor=fig.get_facecolor())
plt.savefig('galaxia_halo_dispersao.pdf',
            bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()