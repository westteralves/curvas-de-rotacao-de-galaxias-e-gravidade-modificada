import os
import zipfile

NOME_ZIP = "Rotmod_LTG.zip"

def extrair_local():
    pasta_data = os.path.dirname(os.path.abspath(__file__))
    caminho_zip = os.path.join(pasta_data, NOME_ZIP)

    if not os.path.exists(caminho_zip):
        print(f"ERRO: O arquivo {NOME_ZIP} não foi encontrado na pasta 'data/'.")
        print("Baixe-o manualmente em:")
        print("https://github.com/arwells-research/compact-fiber-galactic-rotation/blob/main/data/sparc/Rotmod_LTG.zip")
        print("e coloque-o dentro da pasta 'data/'.")
        return

    print(f"Extraindo {NOME_ZIP}...")
    try:
        with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
            zip_ref.extractall(pasta_data)
        print("Extração concluída. Os arquivos .dat estão na pasta 'data/'.")
    except Exception as e:
        print(f"Erro ao extrair: {e}")

if __name__ == "__main__":
    extrair_local()
