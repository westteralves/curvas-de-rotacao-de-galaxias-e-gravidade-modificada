import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os

# Constante de Milgrom: a0 = 1.2e-10 m/s^2 convertido para (km/s)^2 / kpc
# 1.2e-10 * (3.086e13) = 3703.2
a0 = 3703.2


def ler_dados_sparc(arquivo):
    """Tenta ler o arquivo .dat do SPARC. Procura primeiro na pasta 'dados', depois na pasta atual."""
    colunas = ['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    caminhos = [
        os.path.join('dados', arquivo),
        os.path.join('data', arquivo),
        arquivo
    ]
    for caminho in caminhos:
        if os.path.exists(caminho):
            print(f"Lendo arquivo: {caminho}")
            return pd.read_csv(caminho, sep=r'\s+', comment='#', names=colunas)
    raise FileNotFoundError(f"Arquivo {arquivo} não encontrado.")


def v_bario_total(df, Y_disco, Y_bojo):
    """Soma quadrática das contribuições: gás + disco (Y_disco) + bojo (Y_bojo)."""
    v2_gas = df['Vgas']**2
    v2_disco = Y_disco * (df['Vdisk']**2)
    v2_bojo = Y_bojo * (df['Vbul']**2)
    return v2_gas + v2_disco + v2_bojo


def v_mond(r, df, Y_disco, Y_bojo):
    """Velocidade MOND. Interpola a massa bariônica para raios suaves (apenas para visualização)."""
    v2_bar_orig = v_bario_total(df, Y_disco, Y_bojo)
    r_orig = df['Rad'].values

    # Interpola a velocidade bariônica para os novos raios (curva contínua)
    interp_func = interp1d(r_orig, v2_bar_orig, kind='cubic', fill_value='extrapolate')
    v2_bar = interp_func(r)

    # Solução da equação \mu(g/a0) * g = g_newton com \mu(x) = x/(1+x)
    g_newton = v2_bar / r
    g_efetivo = (g_newton + np.sqrt(g_newton**2 + 4*g_newton*a0)) / 2.0
    return np.sqrt(g_efetivo * r)


def v_yukawa(r, df, Y_disco, Y_bojo, alpha, lambd):
    """Velocidade prevista pelo potencial de Yukawa (correção à gravidade newtoniana)."""
    v2_bar_orig = v_bario_total(df, Y_disco, Y_bojo)
    r_orig = df['Rad'].values

    interp_func = interp1d(r_orig, v2_bar_orig, kind='cubic', fill_value='extrapolate')
    v2_bar = interp_func(r)

    # Fator de correção do potencial de Yukawa
    fator = 1 + alpha * (1 + r/lambd) * np.exp(-r/lambd)
    fator = np.maximum(fator, 0.0)  # evita valores negativos (não físicos)
    return np.sqrt(v2_bar * fator)


# Carrega os dados
dados = ler_dados_sparc("UGC02885_rotmod.dat")
r = dados['Rad'].values
v_obs = dados['Vobs'].values
erro = dados['errV'].values

# Ajuste MOND
def mond_fit(r, Y_disco, Y_bojo):
    return v_mond(r, dados, Y_disco, Y_bojo)

# Limites físicos para as razões massa-luminosidade (0.1 a 1.2)
bounds_mond = ([0.1, 0.1], [1.2, 1.2])
popt_mond, _ = curve_fit(mond_fit, r, v_obs, sigma=erro,
                         p0=[0.5, 0.5], bounds=bounds_mond)

# Ajuste Yukawa
def yukawa_fit(r, Y_disco, Y_bojo, alpha, lambd):
    return v_yukawa(r, dados, Y_disco, Y_bojo, alpha, lambd)

bounds_yuk = ([0.1, 0.1, 0.0, 0.1], [1.2, 1.2, 5.0, 1000.0])
popt_yuk, _ = curve_fit(yukawa_fit, r, v_obs, sigma=erro,
                        p0=[0.5, 0.5, 0.6, 50.0], bounds=bounds_yuk)

# Cálculo do chi^2 reduzido (n_params = número de parâmetros livres)
def chi2_red(v_mod, n_params):
    chi2 = np.sum(((v_obs - v_mod) / erro)**2)
    return chi2 / (len(v_obs) - n_params)

v_mond_ajustada = mond_fit(r, *popt_mond)
v_yuk_ajustada = yukawa_fit(r, *popt_yuk)

print("\nRESULTADOS MOND")
print(f"Y_disco = {popt_mond[0]:.3f}, Y_bojo = {popt_mond[1]:.3f}")
print(f"χ²_red = {chi2_red(v_mond_ajustada, 2):.3f}")

print("\nRESULTADOS YUKAWA")
print(f"Y_disco = {popt_yuk[0]:.3f}, Y_bojo = {popt_yuk[1]:.3f}")
print(f"α = {popt_yuk[2]:.3f}, λ = {popt_yuk[3]:.1f} kpc")
print(f"χ²_red = {chi2_red(v_yuk_ajustada, 4):.3f}")

# Gera uma curva contínua com muitos pontos para visualização
r_suave = np.linspace(0.5, 78, 500)

# Interpolação cúbica para suavizar a curva (apenas visual, não afeta os parâmetros)
v_mond_suave = v_mond(r_suave, dados, *popt_mond)
v_yuk_suave = v_yukawa(r_suave, dados, *popt_yuk)

# Curva bariônica suave (Y=1)
v_bar_orig = np.sqrt(v_bario_total(dados, 1.0, 1.0))
interp_bar = interp1d(r, v_bar_orig, kind='cubic', fill_value='extrapolate')
v_bar_suave = interp_bar(r_suave)

# Configuração de fontes para os gráficos
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Gráfico 1: ajuste MOND
plt.figure(figsize=(9, 5), dpi=600)

plt.errorbar(r, v_obs, yerr=erro, fmt='o', color='#f57c00',
             capsize=3, capthick=1.5, elinewidth=1.5,
             markersize=6, alpha=0.8, label='Dados SPARC')

plt.plot(r_suave, v_mond_suave, color='#1a73e8', linewidth=2.5,
         label=f'MOND ($\chi^2_{{\mathrm{{red}}}} = {chi2_red(v_mond_ajustada, 2):.3f}$)')

plt.xlabel('Raio (kpc)', fontsize=13)
plt.ylabel('Velocidade (km/s)', fontsize=13)
plt.title('Ajuste MOND para a Galáxia de Rubin (UGC 2885)', fontsize=14, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.xlim(0, 78)
plt.ylim(240, 315)
plt.legend(fontsize=10, frameon=True, edgecolor='gray', loc='lower right')

plt.savefig('ajuste_mond_ugc2885.pdf', bbox_inches='tight')
plt.savefig('ajuste_mond_ugc2885.png', dpi=600, bbox_inches='tight')
plt.show()

# Gráfico 2: ajuste Yukawa
plt.figure(figsize=(9, 5), dpi=600)

plt.errorbar(r, v_obs, yerr=erro, fmt='o', color='#f57c00',
             capsize=3, capthick=1.5, elinewidth=1.5,
             markersize=6, alpha=0.8, label='Dados SPARC')

plt.plot(r_suave, v_yuk_suave, color='#2e7d32', linewidth=2.5,
         label=f'Yukawa ($\chi^2_{{\mathrm{{red}}}} = {chi2_red(v_yuk_ajustada, 4):.3f}$)')

plt.xlabel('Raio (kpc)', fontsize=13)
plt.ylabel('Velocidade (km/s)', fontsize=13)
plt.title('Ajuste do Potencial de Yukawa para a Galáxia de Rubin (UGC 2885)', fontsize=14, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.xlim(0, 78)
plt.ylim(240, 315)
plt.legend(fontsize=10, frameon=True, edgecolor='gray', loc='lower right')

plt.savefig('ajuste_yukawa_ugc2885.pdf', bbox_inches='tight')
plt.savefig('ajuste_yukawa_ugc2885.png', dpi=600, bbox_inches='tight')
plt.show()

# Gráfico 3: comparação dos modelos com a curva bariônica
plt.figure(figsize=(9, 5), dpi=600)

plt.errorbar(r, v_obs, yerr=erro, fmt='o', color='#f57c00',
             capsize=3, capthick=1.5, elinewidth=1.5,
             markersize=6, alpha=0.8, label='Dados SPARC')

plt.plot(r_suave, v_bar_suave, color='#d32f2f', linewidth=1.8, linestyle='--',
         label='Bariônica')

plt.plot(r_suave, v_mond_suave, color='#1a73e8', linewidth=2.0,
         label=f'MOND ($\chi^2_{{\mathrm{{red}}}} = {chi2_red(v_mond_ajustada, 2):.3f}$)')

plt.plot(r_suave, v_yuk_suave, color='#2e7d32', linewidth=2.0,
         label=f'Yukawa ($\chi^2_{{\mathrm{{red}}}} = {chi2_red(v_yuk_ajustada, 4):.3f}$)')

plt.xlabel('Raio (kpc)', fontsize=13)
plt.ylabel('Velocidade (km/s)', fontsize=13)
plt.title('Comparação dos Ajustes para a Galáxia de Rubin (UGC 2885)', fontsize=14, fontweight='medium')
plt.grid(True, alpha=0.2, linestyle='--')
plt.xlim(0, 78)
plt.ylim(175, 315)
plt.legend(fontsize=10, frameon=True, edgecolor='gray', loc='lower right')

plt.savefig('comparacao_ajustes_ugc2885.pdf', bbox_inches='tight')
plt.savefig('comparacao_ajustes_ugc2885.png', dpi=600, bbox_inches='tight')
plt.show()