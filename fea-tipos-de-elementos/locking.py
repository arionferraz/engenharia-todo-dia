"""
locking.py - Engenharia Todo Dia
Demonstração numérica de shear locking e volumetric locking em elementos
planos (estado plano), comparando:
  Q4  - quadrilátero bilinear, integração completa 2x2
  QM6 - Q4 com modos incompatíveis (Wilson/Taylor), condensação estática
  Q4S - Q4 com integração seletiva reduzida (parcela volumétrica com 1 ponto)
  Q8  - quadrilátero quadrático serendipity, integração 3x3
Problema: viga em balanço L x h, engastada em x=0, carga P na extremidade.
"""
import numpy as np

# ---------------------------------------------------------------- material
def matriz_D(E, nu, estado):
    """Matriz constitutiva elástica (estado plano de tensão ou deformação)."""
    if estado == "EPT":   # estado plano de tensão
        c = E / (1 - nu**2)
        return c * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    c = E / ((1 + nu) * (1 - 2 * nu))   # estado plano de deformação
    return c * np.array([[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2]])

def separa_D(E, nu):
    """Parte desviadora e volumétrica de D (EPD) para integração seletiva."""
    G = E / (2 * (1 + nu)); K = E / (3 * (1 - 2 * nu))
    m = np.array([1.0, 1.0, 0.0])
    Dvol = K * np.outer(m, m)
    Ddev = matriz_D(E, nu, "EPD") - Dvol
    return Ddev, Dvol

# ------------------------------------------------------- funções de forma
def forma_Q4(xi, eta):
    N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
    dN = 0.25 * np.array([[-(1-eta), (1-eta), (1+eta), -(1+eta)],
                          [-(1-xi), -(1+xi), (1+xi), (1-xi)]])
    return N, dN

def forma_Q8(xi, eta):
    s = [(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0)]
    N = np.zeros(8); dN = np.zeros((2, 8))
    for i, (a, b) in enumerate(s):
        if i < 4:   # nós de canto
            N[i] = 0.25*(1+a*xi)*(1+b*eta)*(a*xi+b*eta-1)
            dN[0,i] = 0.25*a*(1+b*eta)*(2*a*xi+b*eta)
            dN[1,i] = 0.25*b*(1+a*xi)*(a*xi+2*b*eta)
        elif a == 0:  # nós de meio de lado horizontais
            N[i] = 0.5*(1-xi**2)*(1+b*eta)
            dN[0,i] = -xi*(1+b*eta); dN[1,i] = 0.5*b*(1-xi**2)
        else:         # nós de meio de lado verticais
            N[i] = 0.5*(1+a*xi)*(1-eta**2)
            dN[0,i] = 0.5*a*(1-eta**2); dN[1,i] = -eta*(1+a*xi)
    return N, dN

def matriz_B(dNdx):
    n = dNdx.shape[1]; B = np.zeros((3, 2*n))
    B[0, 0::2] = dNdx[0]; B[1, 1::2] = dNdx[1]
    B[2, 0::2] = dNdx[1]; B[2, 1::2] = dNdx[0]
    return B

def gauss(n):
    p, w = np.polynomial.legendre.leggauss(n)
    return [(a, b, wa*wb) for a, wa in zip(p, w) for b, wb in zip(p, w)]

# ------------------------------------------------------ rigidez elementar
def k_elem(xy, tipo, D, Dsplit=None):
    if tipo == "Q8":
        K = np.zeros((16, 16))
        for xi, eta, w in gauss(3):
            _, dN = forma_Q8(xi, eta); J = dN @ xy
            B = matriz_B(np.linalg.solve(J, dN)); K += B.T @ D @ B * np.linalg.det(J) * w
        return K
    if tipo == "Q4S":   # integração seletiva: desviadora 2x2, volumétrica 1x1
        Ddev, Dvol = Dsplit; K = np.zeros((8, 8))
        for xi, eta, w in gauss(2):
            _, dN = forma_Q4(xi, eta); J = dN @ xy
            B = matriz_B(np.linalg.solve(J, dN)); K += B.T @ Ddev @ B * np.linalg.det(J) * w
        _, dN = forma_Q4(0, 0); J = dN @ xy
        B = matriz_B(np.linalg.solve(J, dN)); K += B.T @ Dvol @ B * np.linalg.det(J) * 4
        return K
    K = np.zeros((8, 8))
    if tipo == "Q4":
        for xi, eta, w in gauss(2):
            _, dN = forma_Q4(xi, eta); J = dN @ xy
            B = matriz_B(np.linalg.solve(J, dN)); K += B.T @ D @ B * np.linalg.det(J) * w
        return K
    # QM6: modos incompatíveis (1-xi²) e (1-eta²), Jacobiano no centro (Taylor)
    _, dN0 = forma_Q4(0, 0); J0 = dN0 @ xy; detJ0 = np.linalg.det(J0)
    Kuu = np.zeros((8, 8)); Kua = np.zeros((8, 4)); Kaa = np.zeros((4, 4))
    for xi, eta, w in gauss(2):
        _, dN = forma_Q4(xi, eta); J = dN @ xy; dJ = np.linalg.det(J)
        Bu = matriz_B(np.linalg.solve(J, dN))
        dP = np.array([[-2*xi, 0.0], [0.0, -2*eta]])        # d/dxi, d/deta dos modos
        Ba = matriz_B(np.linalg.solve(J0, dP)) * detJ0 / dJ
        Kuu += Bu.T @ D @ Bu * dJ * w; Kua += Bu.T @ D @ Ba * dJ * w
        Kaa += Ba.T @ D @ Ba * dJ * w
    return Kuu - Kua @ np.linalg.solve(Kaa, Kua.T)        # condensação estática

# ------------------------------------------------------- montagem e solução
def balanco(nx, ny, tipo, L=10.0, h=1.0, E=1000.0, nu=0.3, P=1.0, estado="EPT"):
    """Retorna o deslocamento vertical médio da extremidade livre."""
    D = matriz_D(E, nu, estado); Ds = separa_D(E, nu) if tipo == "Q4S" else None
    quad = tipo == "Q8"
    mx, my = (2*nx+1, 2*ny+1) if quad else (nx+1, ny+1)
    X, Y = np.meshgrid(np.linspace(0, L, mx), np.linspace(-h/2, h/2, my), indexing="ij")
    nid = lambda i, j: i*my + j
    ndof = 2*mx*my; K = np.zeros((ndof, ndof))
    for e in range(nx):
        for f in range(ny):
            if quad:
                i, j = 2*e, 2*f
                nos = [nid(i,j), nid(i+2,j), nid(i+2,j+2), nid(i,j+2),
                       nid(i+1,j), nid(i+2,j+1), nid(i+1,j+2), nid(i,j+1)]
            else:
                nos = [nid(e,f), nid(e+1,f), nid(e+1,f+1), nid(e,f+1)]
            xy = np.column_stack([X.ravel()[nos], Y.ravel()[nos]])
            ke = k_elem(xy, tipo, D, Ds)
            gl = np.ravel([[2*n, 2*n+1] for n in nos])
            K[np.ix_(gl, gl)] += ke
    # carga de extremidade: cisalhamento uniforme na face livre, cargas nodais
    # consistentes (1/2-1/2 por aresta linear; 1/6-4/6-1/6 por aresta quadrática)
    F = np.zeros(ndof); borda = [nid(mx-1, j) for j in range(my)]
    q = P / h                                  # tração por unidade de altura
    for f in range(ny):
        if quad:
            nos, pesos = borda[2*f:2*f+3], [1/6, 4/6, 1/6]
        else:
            nos, pesos = borda[f:f+2], [1/2, 1/2]
        for n, p in zip(nos, pesos):
            F[2*n+1] -= q * (h/ny) * p
    presos = np.ravel([[2*nid(0,j), 2*nid(0,j)+1] for j in range(my)])
    usados = np.flatnonzero(np.diag(K) > 0)   # Q8: o nó central da grade não pertence a nenhum elemento
    livres = np.setdiff1d(usados, presos)
    u = np.zeros(ndof); u[livres] = np.linalg.solve(K[np.ix_(livres, livres)], F[livres])
    return -np.mean(u[2*np.array([nid(mx-1, j) for j in range(my)])+1])

def referencia(L=10.0, h=1.0, E=1000.0, nu=0.3, P=1.0, estado="EPT"):
    """Flecha de Timoshenko: flexão + cisalhamento (kappa = 5/6)."""
    Ee = E if estado == "EPT" else E/(1 - nu**2)
    I = h**3/12; G = E/(2*(1+nu))
    return P*L**3/(3*Ee*I) + P*L/((5/6)*G*h)

if __name__ == "__main__":
    ref = referencia()
    print(f"Referência (Timoshenko): {ref:.4f}\n")
    print("Shear locking - estado plano de tensão, nu = 0,3, L/h = 10")
    print(f"{'malha':>8} {'Q4':>7} {'QM6':>7} {'Q8':>7}   (flecha / referência)")
    for nx, ny in [(5,1),(10,1),(10,2),(20,2),(40,4)]:
        r = [balanco(nx, ny, t)/ref for t in ("Q4","QM6","Q8")]
        print(f"{nx:>4}x{ny:<3} " + " ".join(f"{v:7.3f}" for v in r))
    print("\nVolumetric locking - estado plano de deformação, malha 20x4")
    print(f"{'nu':>7} {'Q4':>7} {'Q4S':>7} {'Q8':>7}")
    for nu in (0.3, 0.45, 0.49, 0.499, 0.4999):
        rf = referencia(nu=nu, estado="EPD")
        r = [balanco(20, 4, t, nu=nu, estado="EPD")/rf for t in ("Q4","Q4S","Q8")]
        print(f"{nu:7.4f} " + " ".join(f"{v:7.3f}" for v in r))
