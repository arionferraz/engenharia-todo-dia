# Tipos de elementos finitos: shear locking, volumetric locking e convergência de malha

Código de apoio do artigo **[Tipos de elementos finitos 1D, 2D e 3D: formulações, travamento numérico e critérios de escolha no ANSYS](https://engenhariatododia.com.br/tipos-de-elementos-finitos-1d-2d-3d-locking-ansys/)**, do blog Engenharia Todo Dia.

Um código de elementos finitos em estado plano, escrito do zero em Python/NumPy, que reproduz os números e as animações do artigo.

## Arquivos

| Arquivo | O que faz |
| --- | --- |
| `locking.py` | Viga em balanço (L/h = 10) com quatro formulações de quadrilátero: **Q4** (bilinear, 2×2), **QM6** (modos incompatíveis com condensação estática), **Q4S** (integração seletiva) e **Q8** (serendipity, 3×3). Mede o shear locking em 5 malhas e o volumetric locking com ν de 0,3 a 0,4999. |
| `convergencia.py` | Estudo de convergência com extrapolação de Richardson e **Grid Convergence Index (GCI)** de Roache, aplicado à mesma viga em três malhas (r = 2). |
| `saida_esperada.txt` | Saída dos dois scripts, para conferir a sua execução. |

## Como rodar

```bash
pip install -r requirements.txt
python locking.py
python convergencia.py   # importa locking.py, mantenha os dois na mesma pasta
```

Requer Python 3.8+ e NumPy. Cada execução leva poucos segundos (a malha mais fina tem 81 × 9 nós).

## Resultado principal

Flecha calculada ÷ referência de Timoshenko (4,0312), estado plano de tensão, ν = 0,3:

| Malha | Q4 | QM6 | Q8 |
| --- | --- | --- | --- |
| 5×1 | 0,379 | 0,989 | 0,981 |
| 10×1 | 0,674 | 0,996 | 0,991 |
| 20×2 | 0,888 | 0,995 | 0,997 |
| 40×4 | 0,968 | 0,997 | 0,998 |

Com ν = 0,4999 em estado plano de deformação (malha 20×4), o Q4 com integração completa cai para **0,021**, enquanto o Q4S e o Q8 ficam em torno de 0,95.

## Limitações

- Código didático: matriz de rigidez densa e solver direto do NumPy. Não é feito para malhas grandes.
- Carga na extremidade aplicada como cisalhamento uniforme com cargas nodais consistentes; engaste com todos os nós da face x = 0 presos.
- A referência é a solução de viga de Timoshenko (κ = 5/6), não uma solução de elasticidade 2D exata, o que explica razões levemente diferentes de 1 mesmo nas malhas finas.

## Licença

MIT. Se usar em trabalho ou aula, cite o artigo.
