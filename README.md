# Curvas de Rotação de Galáxias e Gravidade Modificada

Códigos em Python desenvolvidos para o meu Trabalho de Conclusão de Curso (TCC) em Bacharelado em Física na Universidade Federal do Oeste da Bahia (UFOB).

* **Autor:** Westter Alves Vieira
* **Orientador:** Prof. Dr. Elias Brito Alves Junior
* **Ano:** 2026

## Sobre o projeto

O trabalho compara dois modelos de gravidade modificada na descrição de curvas de rotação de galáxias espirais:

* **MOND** (Dinâmica Newtoniana Modificada), utilizando a função de interpolação $\mu(x) = x/(1+x)$ e a escala $a_0 = 1{,}2 \times 10^{-10}\,\text{m/s}^2$.
* **Potencial de Yukawa**, inspirado no limite de campo fraco das teorias $f(R)$, com parâmetros de acoplamento $\alpha$ e alcance $\lambda$.

Os ajustes foram aplicados à galáxia UGC 2885 (Galáxia de Rubin) e estendidos para as 175 galáxias do catálogo SPARC.

## Requisitos e Instalação

É necessário ter o Python 3.10 ou superior instalado. As dependências do projeto (`numpy`, `pandas`, `scipy`, `matplotlib`) podem ser instaladas com:

```bash
pip install -r requirements.txt

