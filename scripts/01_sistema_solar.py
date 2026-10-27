import numpy as np
import matplotlib.pyplot as plt

# Dados dos planetas extraidos da NASA Fact Sheets
distancias_ua = np.array([0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.191, 30.069])
velocidades_obs = np.array([47.87, 35.02, 29.78, 24.07, 13.07, 9.69, 6.81, 5.43])

# Cores dos planetas
cores = {
    'Mercúrio': '#b0b0b0',
    'Vênus': '#e6b800',
    'Terra': '#2e8b57',
    'Marte': '#c1440e',
    'Júpiter': '#d4a373',
    'Saturno': '#d4af37',
    'Urano': '#4fc3f7',
    'Netuno': '#1a237e'
}
nomes = list(cores.keys())

# Constantes
G = 6.67430e-11
M_sol = 1.989e30
UA = 1.496e11   # 1 UA em metros

# Para os gráficos, os pontos estão entre 0.1 e 35 UA, assim a curva fica visivelmente suave
raios_ua = np.linspace(0.1, 35, 1000)
raios_metros = raios_ua * UA   # converte para metros

# Velocidade kepleriana discutida no TCC: v(r) = sqrt(GM/r)
# O resultado é em m/s, dividimos por 1000 para ter km/s
vel_teorica = np.sqrt(G * M_sol / raios_metros) / 1000

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

plt.figure(figsize=(9, 5), dpi=300)

# Curva teórica
plt.plot(raios_ua, vel_teorica, color='#1a73e8', linewidth=2.5,
         label=r'$V(r) = \sqrt{GM_{\odot}/r}$')

# Atribuição da cor ao planeta
for i, planeta in enumerate(nomes):
    plt.scatter(distancias_ua[i], velocidades_obs[i], color=cores[planeta], s=90,
                edgecolor='black', linewidth=0.8, zorder=5,
                label=planeta)

# Rótulos e grade
plt.xlabel('Distância ao Sol (UA)', fontsize=13)
plt.ylabel('Velocidade Orbital (km/s)', fontsize=13)
plt.title('Curva de Rotação do Sistema Solar', fontsize=15, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.ylim(0, 52)

# Legenda do gráfico
plt.legend(fontsize=9, frameon=True, edgecolor='gray', ncol=2)

plt.savefig('curva_rotacao_solar.pdf', bbox_inches='tight')
plt.savefig('curva_rotacao_solar.png', dpi=600, bbox_inches='tight')
plt.show()