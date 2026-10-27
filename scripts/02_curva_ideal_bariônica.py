import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os


def ler_dados_sparc(arquivo):
    colunas = ['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    caminhos = [os.path.join('dados', arquivo),
                os.path.join('data', arquivo),
                arquivo]
    for caminho in caminhos:
        if os.path.exists(caminho):
            print(f"Lendo arquivo: {caminho}")
            return pd.read_csv(caminho, sep=r'\s+', comment='#', names=colunas)
    raise FileNotFoundError(f"Arquivo {arquivo} não encontrado.")


def v_bario_total(df, Y_disco, Y_bojo):
    return (df['Vgas']**2
            + Y_disco * df['Vdisk']**2
            + Y_bojo * df['Vbul']**2)


# Leitura dos dados observacionais
dados = ler_dados_sparc("UGC02885_rotmod.dat")
r = dados['Rad'].values
v_obs = dados['Vobs'].values
erro = dados['errV'].values

# Grade suave para as curvas contínuas
r_suave = np.linspace(0.01, 78, 800)

# Curva ideal: interpolação cúbica dos dados observados
interp_ideal = interp1d(r, v_obs, kind='cubic', fill_value='extrapolate')
v_ideal = interp_ideal(r_suave)

# Curva bariônica com razões massa-luminosidade unitárias
v_bar = interp1d(r, np.sqrt(v_bario_total(dados, 1.0, 1.0)),
                 kind='cubic', fill_value='extrapolate')(r_suave)

# Configuração de estilo
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Gráfico
plt.figure(figsize=(9, 5), dpi=600)

plt.errorbar(r, v_obs, yerr=erro, fmt='o', color='#f57c00',
             capsize=3, capthick=1.5, elinewidth=1.5,
             markersize=6, alpha=0.8, label='Dados SPARC')

plt.plot(r_suave, v_ideal, color='#1a237e', linewidth=1.8,
         label='Curva ideal (Interpolada)')

plt.plot(r_suave, v_bar, color='#d32f2f', linewidth=1.5,
         linestyle='--', label='Contribuição bariônica')

plt.xlabel('Raio (kpc)', fontsize=13)
plt.ylabel('Velocidade (km/s)', fontsize=13)
plt.title('Curva Ideal e Contribuição Bariônica para a Galáxia de Rubin (UGC 2885)',
          fontsize=14, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.xlim(0, 78)
plt.ylim(175, 315)
plt.legend(fontsize=10, frameon=True, edgecolor='gray', loc='lower right')

plt.savefig('grafico_curva_ideal_e_massa_barionica_ugc2885.pdf', bbox_inches='tight')
plt.savefig('grafico_curva_ideal_e_massa_barionica_ugc2885.pdf', dpi=600, bbox_inches='tight')
plt.show()