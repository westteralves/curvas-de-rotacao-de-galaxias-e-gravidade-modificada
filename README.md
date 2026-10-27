Curvas de Rotação de Galáxias e Gravidade Modificada

![DOI](https://zenodo.org/badge/10.5281/zenodo.23681694.svg)

Códigos em Python desenvolvidos para o meu Trabalho de Conclusão de Curso (TCC) em Bacharelado em Física na Universidade Federal do Oeste da Bahia (UFOB).

Autor: Westter Alves Vieira
Orientador: Prof. Dr. Elias Brito Alves Junior
Ano: 2026
Sobre o projeto

O trabalho compara dois modelos de gravidade modificada na descrição de curvas de rotação de galáxias espirais:

MOND (Dinâmica Newtoniana Modificada): utiliza a função de interpolação μ(x)=x/(1+x) e a escala a0=1,2×10−10 m/s2.
Potencial de Yukawa: inspirado no limite de campo fraco das teorias f(R), com parâmetros de acoplamento α e alcance λ.

Os ajustes foram aplicados à galáxia UGC 2885 (Galáxia de Rubin) e estendidos para as 175 galáxias do catálogo SPARC.

Requisitos e instalação

É necessário ter o Python 3.10 ou superior instalado.

As dependências do projeto (numpy, pandas, scipy e matplotlib) podem ser instaladas com:

bash
pip install -r requirements.txt

Como reproduzir os resultados
1. Baixar e extrair os dados do SPARC

Execute o script para baixar e extrair as curvas de rotação oficiais na pasta data/:

bash
python data/download_sparc.py

2. Executar os scripts

Os códigos estão na pasta scripts/ e correspondem à ordem do Apêndice B do TCC.

Arquivo	O que o código faz	Seção do TCC
01_sistema_solar.py	Gera a curva de rotação do Sistema Solar	B.1
02_curva_ideal_bariônica.py	Plota a curva ideal e a contribuição bariônica da UGC 2885	B.2
03_burkert_nfw.py	Ajusta os perfis de Burkert e NFW à curva da UGC 2885	B.3
04_galaxia_halo.py Gera a representação da galáxia com halo de matéria escura	B.4
05_ugc2885_mond_yukawa.py	Ajusta os modelos MOND e Yukawa para a UGC 2885	B.5
06_sparc_analise_completa.py	Processa as 175 galáxias do SPARC e gera gráficos estatísticos	B.6
Exemplo de execução

Para executar a análise MOND e Yukawa da galáxia UGC 2885:

bash
python scripts/05_ugc2885_mond_yukawa.py

Citação

Se este código for útil para sua pesquisa, cite o TCC:

bibtex
@misc{vieira2026tcc,
  author       = {Vieira, Westter Alves},
  title        = {Curvas de Rotação de Galáxias e Gravidade Modificada},
  year         = {2026},
  howpublished = {Trabalho de Conclusão de Curso, Universidade Federal do Oeste da Bahia},
  doi          = {10.5281/zenodo.23681694},
  url          = {https://doi.org/10.5281/zenodo.23681694}
}

Créditos e dados de terceiros

Os dados utilizados são provenientes do catálogo SPARC e pertencem aos seus autores originais.

Ao utilizar esses dados, cite o artigo original:

Lelli, F.; McGaugh, S. S.; Schombert, J. M. SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. The Astronomical Journal, v. 152, n. 6, p. 157, 2016.

Licença

Este projeto está sob a licença MIT. Consulte o arquivo LICENSE para mais detalhes.
