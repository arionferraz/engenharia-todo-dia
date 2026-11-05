"""
convergencia.py - Engenharia Todo Dia
Estudo de convergência de malha com extrapolação de Richardson e
Grid Convergence Index (GCI, Roache). Serve para qualquer grandeza escalar
(deslocamento, tensão em um ponto, frequência natural) obtida em 3 malhas
com razão de refino constante r.
"""
import numpy as np

def richardson_gci(f1, f2, f3, r, Fs=1.25):
    """f1 = malha mais fina, f3 = mais grossa; r = h_grossa / h_fina (constante).
    Retorna ordem observada p, valor extrapolado e GCI da malha fina [%]."""
    e32, e21 = f3 - f2, f2 - f1
    if e32 * e21 <= 0:
        raise ValueError("Convergência oscilatória ou divergente: GCI não se aplica.")
    p = np.log(e32 / e21) / np.log(r)                 # ordem de convergência observada
    f_ext = f1 + (f1 - f2) / (r**p - 1)               # extrapolação de Richardson
    gci = Fs * abs((f1 - f2) / f1) / (r**p - 1) * 100 # banda de incerteza da malha fina
    return p, f_ext, gci

if __name__ == "__main__":
    from locking import balanco, referencia
    ref = referencia()
    malhas = [(10, 1), (20, 2), (40, 4)]              # refino uniforme, r = 2
    for tipo in ("Q4", "Q8"):
        f3, f2, f1 = [balanco(nx, ny, tipo) for nx, ny in malhas]
        p, fx, gci = richardson_gci(f1, f2, f3, r=2.0)
        print(f"{tipo}: f(10x1)={f3:.4f} f(20x2)={f2:.4f} f(40x4)={f1:.4f} | "
              f"p={p:.2f} | extrapolado={fx:.4f} ({fx/ref:.3f} da ref.) | GCI={gci:.2f}%")
