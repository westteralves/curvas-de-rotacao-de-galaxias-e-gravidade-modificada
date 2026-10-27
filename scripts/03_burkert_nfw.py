import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os

G = 4.302e-6


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


def massa_burkert(r, rho0, rc):
    x = r / rc
    return 2 * np.pi * rho0 * rc**3 * (
        np.log(1 + x) + 0.5 * np.log(1 + x**2) - np.arctan(x))


def velocidade_burkert(r, df, Y_disco, Y_bojo, rho0, rc):
    v2_bar = interp1d(df['Rad'].values,
                      v_bario_total(df, Y_disco, Y_bojo),
                      kind='cubic', fill_value='extrapolate')(r)
    v2_halo = np.maximum(G * massa_burkert(r, rho0, rc) / r, 0.0)
    return np.sqrt(v2_bar + v2_halo)


def massa_nfw(r, rho0, rs):
    x = r / rs
    return 4 * np.pi * rho0 * rs**3 * (np.log(1 + x) - x / (1 + x))


def velocidade_nfw(r, df, Y_disco, Y_bojo, rho0, rs):
    v2_bar = interp1d(df['Rad'].values,
                      v_bario_total(df, Y_disco, Y_bojo),
                      kind='cubic', fill_value='extrapolate')(r)
    v2_halo = np.maximum(G * massa_nfw(r, rho0, rs) / r, 0.0)
    return np.sqrt(v2_bar + v2_halo)


# Leitura dos dados observacionais
dados = ler_dados_sparc("UGC02885_rotmod.dat")
r = dados['Rad'].values
v_obs = dados['Vobs'].values
erro = dados['errV'].values


def ajuste_burkert(r, Y_disco, Y_bojo, rho0, rc):
    return velocidade_burkert(r, dados, Y_disco, Y_bojo, rho0, rc)


def ajuste_nfw(r, Y_disco, Y_bojo, rho0, rs):
    return velocidade_nfw(r, dados, Y_disco, Y_bojo, rho0, rs)


limites = ([0.1, 0.1, 1e5, 0.1], [1.2, 1.2, 1e10, 100.0])

popt_burk, _ = curve_fit(ajuste_burkert, r, v_obs, sigma=erro,
                         p0=[0.5, 0.5, 5e6, 8.0],
                         bounds=limites, maxfev=20000)

popt_nfw, _ = curve_fit(ajuste_nfw, r, v_obs, sigma=erro,
                        p0=[0.5, 0.5, 5e6, 15.0],
                        bounds=limites, maxfev=20000)

print("\nBurkert")
print(f"  Y_disco = {popt_burk[0]:.3f}")
print(f"  Y_bojo  = {popt_burk[1]:.3f}")
print(f"  rho0    = {popt_burk[2]:.3e} M_sun/kpc^3")
print(f"  r_c     = {popt_burk[3]:.3f} kpc")

print("\nNFW")
print(f"  Y_disco = {popt_nfw[0]:.3f}")
print(f"  Y_bojo  = {popt_nfw[1]:.3f}")
print(f"  rho0    = {popt_nfw[2]:.3e} M_sun/kpc^3")
print(f"  r_s     = {popt_nfw[3]:.3f} kpc")


# Grade suave para as curvas contínuas
r_suave = np.linspace(0.01, 78, 800)

v_bar = interp1d(r, np.sqrt(v_bario_total(dados, 1.0, 1.0)),
                 kind='cubic', fill_value='extrapolate')(r_suave)

v_burk = velocidade_burkert(r_suave, dados, *popt_burk)
v_nfw = velocidade_nfw(r_suave, dados, *popt_nfw)

# Configuração de estilo
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Gráfico
plt.figure(figsize=(9, 5), dpi=600)

plt.errorbar(r, v_obs, yerr=erro, fmt='o', color='#f57c00',
             capsize=3, capthick=1.5, elinewidth=1.5,
             markersize=6, alpha=0.8, label='Dados SPARC')

plt.plot(r_suave, v_bar, color='#d32f2f', linewidth=1.5,
         linestyle='--', label='Bariônica')

plt.plot(r_suave, v_burk, color='#ffeb3b', linewidth=1.8,
         label='Burkert + bariônica')

plt.plot(r_suave, v_nfw, color='#757575', linewidth=1.8,
         label='NFW + bariônica')

plt.xlabel('Raio (kpc)', fontsize=13)
plt.ylabel('Velocidade (km/s)', fontsize=13)
plt.title('Modelos de Matéria Escura para a Galáxia de Rubin (UGC 2885)',
          fontsize=14, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.xlim(0, 78)
plt.ylim(175, 315)
plt.legend(fontsize=10, frameon=True, edgecolor='gray', loc='lower right')

plt.savefig('grafico_modelos_burkert_e_nfw_ugc2885.pdf', bbox_inches='tight')
plt.savefig('grafico_modelos_burkert_e_nfw_ugc2885.png', dpi=600, bbox_inches='tight')
plt.show()