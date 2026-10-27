import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os

# Aceleração de Milgrom em (km/s)^2 / kpc
a0 = 3703.2

# Y_disco e Y_bojo correspondem às razões massa-luminosidade Upsilon_disco e Upsilon_bojo do texto.


def ler_arquivo_sparc(nome_arquivo):
    """Lê um arquivo .dat do catálogo SPARC e retorna um DataFrame com os dados."""
    colunas = ['Raio', 'Vobs', 'errV', 'Vgas', 'Vdisco', 'Vbojo', 'SBdisco', 'SBbojo']
    try:
        dados = pd.read_csv(nome_arquivo, sep=r'\s+', comment='#', names=colunas)
        dados = dados.replace([np.inf, -np.inf], np.nan).dropna()
        return dados
    except Exception as e:
        print(f"  Erro ao ler {nome_arquivo}: {e}")
        return None


def quadrado_velocidade_barionica(dados, Y_disco, Y_bojo):
    """Retorna o quadrado da velocidade circular bariônica total."""
    v2_gas = dados['Vgas'].values**2
    v2_disco = Y_disco * (dados['Vdisco'].values**2)
    v2_bojo = Y_bojo * (dados['Vbojo'].values**2)
    return v2_gas + v2_disco + v2_bojo


def velocidade_mond(raio, dados, Y_disco, Y_bojo):
    """Velocidade circular prevista pela MOND com μ(x) = x/(1+x)."""
    v2_bar = quadrado_velocidade_barionica(dados, Y_disco, Y_bojo)
    g_newton = v2_bar / raio
    g_efetivo = (g_newton + np.sqrt(g_newton**2 + 4 * g_newton * a0)) / 2.0
    return np.sqrt(g_efetivo * raio)


def velocidade_yukawa(raio, dados, Y_disco, Y_bojo, alpha, lambd):
    """Velocidade circular prevista pelo potencial de Yukawa."""
    v2_bar = quadrado_velocidade_barionica(dados, Y_disco, Y_bojo)
    fator = 1 + alpha * (1 + raio / lambd) * np.exp(-raio / lambd)
    fator = np.maximum(fator, 0.0)
    return np.sqrt(v2_bar * fator)


def ajustar_galaxia(dados, nome_galaxia):
    """Realiza os ajustes MOND e Yukawa para uma galáxia."""
    raio = dados['Raio'].values
    v_obs = dados['Vobs'].values
    erro = dados['errV'].values

    if len(raio) < 4:
        return None

    def mond_ajuste(raio, Yd, Yb):
        return velocidade_mond(raio, dados, Yd, Yb)

    try:
        params_mond, _ = curve_fit(mond_ajuste, raio, v_obs, sigma=erro,
                                   p0=[0.5, 0.5],
                                   bounds=([0.1, 0.1], [1.2, 1.2]),
                                   maxfev=5000)
        v_mond = mond_ajuste(raio, *params_mond)
        chi2_mond = np.sum(((v_obs - v_mond) / erro)**2)
        chi2_red_mond = chi2_mond / (len(raio) - 2) if len(raio) > 2 else np.nan
    except Exception:
        params_mond = [np.nan, np.nan]
        chi2_red_mond = np.nan

    def yukawa_ajuste(raio, Yd, Yb, alpha, lambd):
        return velocidade_yukawa(raio, dados, Yd, Yb, alpha, lambd)

    try:
        params_yuk, _ = curve_fit(yukawa_ajuste, raio, v_obs, sigma=erro,
                                  p0=[0.5, 0.5, 0.6, 50.0],
                                  bounds=([0.1, 0.1, 0.0, 0.1],
                                          [1.2, 1.2, 5.0, 1000.0]),
                                  maxfev=5000)
        v_yuk = yukawa_ajuste(raio, *params_yuk)
        chi2_yuk = np.sum(((v_obs - v_yuk) / erro)**2)
        chi2_red_yuk = chi2_yuk / (len(raio) - 4) if len(raio) > 4 else np.nan
    except Exception:
        params_yuk = [np.nan, np.nan, np.nan, np.nan]
        chi2_red_yuk = np.nan

    return {
        'galaxia': nome_galaxia,
        'Y_disco_mond': params_mond[0],
        'Y_bojo_mond': params_mond[1],
        'chi2_red_mond': chi2_red_mond,
        'Y_disco_yuk': params_yuk[0],
        'Y_bojo_yuk': params_yuk[1],
        'alpha': params_yuk[2],
        'lambda': params_yuk[3],
        'chi2_red_yuk': chi2_red_yuk,
        'n_pontos': len(raio)
    }


def processar_todas_galaxias(pasta_dados='dados'):
    """Processa todos os arquivos .dat do SPARC em uma pasta."""
    if not os.path.exists(pasta_dados):
        print(f"ERRO: A pasta '{pasta_dados}' não existe!")
        return pd.DataFrame()

    arquivos = glob.glob(os.path.join(pasta_dados, '*.dat'))
    if not arquivos:
        arquivos = glob.glob(os.path.join(pasta_dados, '*.txt'))
    if not arquivos:
        arquivos = [f for f in glob.glob(os.path.join(pasta_dados, '*')) if os.path.isfile(f)]

    print(f"Encontrados {len(arquivos)} arquivos.")
    if len(arquivos) == 0:
        print("Nenhum arquivo encontrado.")
        return pd.DataFrame()

    resultados = []
    for i, arq in enumerate(arquivos):
        nome = os.path.basename(arq).replace('_rotmod.dat', '').replace('_rotmod.txt', '')
        if nome == os.path.basename(arq):
            nome = os.path.splitext(os.path.basename(arq))[0]

        print(f"[{i+1}/{len(arquivos)}] Processando {nome}...")
        dados = ler_arquivo_sparc(arq)
        if dados is None or len(dados) < 4:
            print(f"  Pulando {nome} (dados insuficientes).")
            continue

        res = ajustar_galaxia(dados, nome)
        if res is not None:
            resultados.append(res)

        if (i + 1) % 10 == 0:
            pd.DataFrame(resultados).to_csv('todos_resultados_parciais.csv', index=False)

    df_resultados = pd.DataFrame(resultados)
    df_resultados.to_csv('todos_resultados.csv', index=False)
    print(f"Processamento concluído. {len(df_resultados)} galáxias analisadas.")
    return df_resultados


def exibir_estatisticas(df):
    """Exibe estatísticas detalhadas dos ajustes e exporta chi2_red_todos.csv."""
    print("\n" + "=" * 60)
    print("ESTATÍSTICAS GERAIS")
    print("=" * 60)
    print(f"Número total de galáxias processadas: {len(df)}")
    mond_validos = df['chi2_red_mond'].notna().sum()
    yuk_validos = df['chi2_red_yuk'].notna().sum()
    print(f"Galáxias com chi2 MOND válido: {mond_validos}")
    print(f"Galáxias com chi2 Yukawa válido: {yuk_validos}")

    print("\n--- MOND ---")
    mond = df['chi2_red_mond'].dropna()
    print(f"  Média de chi2_red = {mond.mean():.3f}")
    print(f"  Mediana de chi2_red = {mond.median():.3f}")
    print(f"  Desvio padrão = {mond.std():.3f}")
    print(f"  Fração com chi2_red < 2 = {(mond < 2).mean() * 100:.1f}%")
    print(f"  Fração com chi2_red < 5 = {(mond < 5).mean() * 100:.1f}%")
    print(f"  Fração com chi2_red < 10 = {(mond < 10).mean() * 100:.1f}%")

    print("\n--- Yukawa (alpha livre) ---")
    yuk = df['chi2_red_yuk'].dropna()
    print(f"  Média de chi2_red = {yuk.mean():.3f}")
    print(f"  Mediana de chi2_red = {yuk.median():.3f}")
    print(f"  Desvio padrão = {yuk.std():.3f}")
    print(f"  Fração com chi2_red < 2 = {(yuk < 2).mean() * 100:.1f}%")
    print(f"  Fração com chi2_red < 5 = {(yuk < 5).mean() * 100:.1f}%")
    print(f"  Fração com chi2_red < 10 = {(yuk < 10).mean() * 100:.1f}%")

    print("\n" + "=" * 60)
    print("ESTATÍSTICAS COM FILTRO (chi2_red < 10)")
    print("=" * 60)
    mask_mond = df['chi2_red_mond'] < 10
    mask_yuk = df['chi2_red_yuk'] < 10
    print(f"Galáxias com MOND chi2 < 10: {mask_mond.sum()}/{len(df)}")
    print(f"Galáxias com Yukawa chi2 < 10: {mask_yuk.sum()}/{len(df)}")
    if mask_mond.sum() > 0:
        mond_filt = df.loc[mask_mond, 'chi2_red_mond']
        print(f"  MOND (filtrado): média={mond_filt.mean():.3f}, mediana={mond_filt.median():.3f}")
    if mask_yuk.sum() > 0:
        yuk_filt = df.loc[mask_yuk, 'chi2_red_yuk']
        print(f"  Yukawa (filtrado): média={yuk_filt.mean():.3f}, mediana={yuk_filt.median():.3f}")

    print("\n" + "=" * 60)
    print("PARÂMETROS MOND")
    print("=" * 60)
    Yd = df['Y_disco_mond'].dropna()
    Yb = df['Y_bojo_mond'].dropna()
    print(f"  Y_disco: média={Yd.mean():.3f}, mediana={Yd.median():.3f}, std={Yd.std():.3f}")
    print(f"  Y_bojo:  média={Yb.mean():.3f}, mediana={Yb.median():.3f}, std={Yb.std():.3f}")

    print("\n--- Parâmetros Yukawa (válidos) ---")
    alpha = df['alpha'].dropna()
    lambd = df['lambda'].dropna()
    print(f"  alpha: média={alpha.mean():.3f}, mediana={alpha.median():.3f}, std={alpha.std():.3f}")
    print(f"  lambda: média={lambd.mean():.1f}, mediana={lambd.median():.1f}, std={lambd.std():.1f}")

    chi2_df = df[['galaxia', 'chi2_red_mond', 'chi2_red_yuk']].copy()
    chi2_df.to_csv('chi2_red_todos.csv', index=False)
    print("\nArquivo 'chi2_red_todos.csv' gerado com todos os chi2_red.")


def plotar_resultados(df):
    """Gera os três gráficos da análise estatística."""
    if df.empty:
        print("Sem dados para plotar.")
        return

    df_limpo = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['chi2_red_mond', 'chi2_red_yuk'])

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 13
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.dpi'] = 600

    mond_chi2 = df_limpo['chi2_red_mond'].dropna()
    yuk_chi2 = df_limpo['chi2_red_yuk'].dropna()

    med_mond = mond_chi2.median()
    med_yuk = yuk_chi2.median()
    mean_mond = mond_chi2.mean()
    mean_yuk = yuk_chi2.mean()

    # Histograma
    fig1, ax1 = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 10, 50)
    ax1.hist(mond_chi2, bins=bins, alpha=0.7, color='#1a73e8',
             edgecolor='black', linewidth=0.5,
             label=f'MOND (mediana = {med_mond:.2f})')
    ax1.hist(yuk_chi2, bins=bins, alpha=0.7, color='#2e7d32',
             edgecolor='black', linewidth=0.5,
             label=f'Yukawa (mediana = {med_yuk:.2f})')
    ax1.set_xlabel('$\\chi^2_{\\mathrm{red}}$', fontsize=13)
    ax1.set_ylabel('Número de galáxias', fontsize=13)
    ax1.set_title('Comparação da qualidade dos ajustes ($\\chi^2_{\\mathrm{red}} \\leq 10$)', fontsize=14, fontweight='medium')
    ax1.grid(True, alpha=0.2, linestyle='--')
    ax1.legend(frameon=True, edgecolor='gray', loc='upper right')
    ax1.set_xlim(0, 10)
    ax1.set_xticks(np.arange(0, 10.5, 0.5))
    fig1.tight_layout()
    fig1.savefig('histograma_qualidade.pdf', bbox_inches='tight')
    fig1.savefig('histograma_qualidade.png', bbox_inches='tight')
    plt.show()

    # Pirâmide
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 20, 30)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    hist_mond, _ = np.histogram(mond_chi2, bins=bins)
    hist_yuk, _ = np.histogram(yuk_chi2, bins=bins)
    ax2.barh(bin_centers, -hist_yuk, height=bins[1]-bins[0],
             color='#2e7d32', alpha=0.7, edgecolor='black', linewidth=0.5,
             label='Yukawa')
    ax2.barh(bin_centers, hist_mond, height=bins[1]-bins[0],
             color='#1a73e8', alpha=0.7, edgecolor='black', linewidth=0.5,
             label='MOND')
    ax2.axvline(x=0, color='black', linewidth=1.5, linestyle='-')
    ax2.axvline(x=mean_mond, color='#1a73e8', linestyle='--', linewidth=1.5, alpha=0.8,
                label=f'Média MOND = {mean_mond:.2f}')
    ax2.axvline(x=-mean_yuk, color='#2e7d32', linestyle='--', linewidth=1.5, alpha=0.8,
                label=f'Média Yukawa = {mean_yuk:.2f}')
    ax2.set_xlabel('Número de galáxias', fontsize=13)
    ax2.set_ylabel('$\\chi^2_{\\mathrm{red}}$', fontsize=13)
    ax2.set_title('Distribuição de $\\chi^2_{\\mathrm{red}}$ para os dois modelos', fontsize=14, fontweight='medium')
    ax2.grid(True, alpha=0.2, linestyle='--', axis='x')
    max_count = max(hist_mond.max(), hist_yuk.max())
    ax2.set_xlim(-max_count * 1.15, max_count * 1.15)
    ax2.set_ylim(0, 20)
    ax2.set_yticks(np.arange(0, 21, 1))
    passo_x = 4
    x_ticks = np.arange(-max_count, max_count + passo_x, passo_x)
    if 0 not in x_ticks:
        x_ticks = np.concatenate([x_ticks, [0]])
        x_ticks = np.sort(x_ticks)
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels([str(abs(int(tick))) for tick in x_ticks])
    ax2.legend(frameon=True, edgecolor='gray', loc='upper right')
    fig2.tight_layout()
    fig2.savefig('piramide_distribuicao.pdf', bbox_inches='tight')
    fig2.savefig('piramide_distribuicao.png', bbox_inches='tight')
    plt.show()

    # Correlação
    fig3, ax3 = plt.subplots(figsize=(9, 6.5))
    x = df_limpo['Y_disco_mond']
    y = df_limpo['chi2_red_mond']
    ax3.scatter(x, y, alpha=0.6, s=30, color='#1a73e8', edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('$\\Upsilon_{\\mathrm{disco}}$ (MOND)', fontsize=13)
    ax3.set_ylabel('$\\chi^2_{\\mathrm{red}}$', fontsize=13)
    ax3.set_title('Correlação entre $\\Upsilon_{\\mathrm{disco}}$ e qualidade do ajuste MOND', fontsize=14, fontweight='medium')
    ax3.grid(True, alpha=0.2, linestyle='--')
    y_max = y.max()
    y_lim = y_max + 0.05 * y_max if y_max > 0 else 1.0
    if y_lim < 22:
        y_lim = 22
    ax3.set_xlim(0, max(x) + 0.1)
    ax3.set_ylim(0, y_lim)
    fig3.tight_layout()
    fig3.savefig('correlacao_yd_chi2.pdf', bbox_inches='tight')
    fig3.savefig('correlacao_yd_chi2.png', bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    caminho = 'dados'
    print(f"Procurando arquivos em: {os.path.abspath(caminho)}")
    df_resultados = processar_todas_galaxias(pasta_dados=caminho)

    if df_resultados.empty:
        print("Nenhuma galáxia foi processada. Verifique os arquivos.")
    else:
        exibir_estatisticas(df_resultados)
        plotar_resultados(df_resultados)
        print("\nArquivos gerados:")
        print("  - todos_resultados.csv")
        print("  - chi2_red_todos.csv")
        print("  - histograma_qualidade.pdf/png")
        print("  - piramide_distribuicao.pdf/png")
        print("  - correlacao_yd_chi2.pdf/png")