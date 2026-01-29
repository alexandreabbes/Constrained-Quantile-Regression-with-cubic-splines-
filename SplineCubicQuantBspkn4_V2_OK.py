import numpy as np
import cvxpy as cp
from scipy.interpolate import BSpline, PPoly
import warnings
import matplotlib.pyplot as plt

def rhotau(u, tau):
    """Fonction de perte quantile"""
    return cp.sum(cp.maximum(tau * u, (tau - 1) * u))

def build_bsplines_and_deriv(knots, degree=4):
    """
    Construit les B-splines de degré 4 et leurs dérivées
    """
    kn = len(knots) - 1
    N = kn + degree  # Nombre de fonctions de base
    
    l_end = np.min(knots)
    r_end = np.max(knots)
    
    # Extension de la séquence de nœuds
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    List_Bsplines = []
    # Pour les contraintes de convexité (dérivée seconde de degré 2)
    DPi2 = np.zeros((N, 3, kn))  # Pour dérivée seconde (quadratique)
    # Pour les contraintes de monotonie (dérivée première de degré 3)
    DPi1 = np.zeros((N, 4, kn))  # Pour dérivée première (cubique)
    
    # Matrice pour contraintes de convexité aux nœuds
    con_array = np.zeros((kn + 1, N))
    
    # Création des B-splines de base
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1
        spline_j = BSpline(s, coefs, degree)
        List_Bsplines.append(spline_j)
        
        # Dérivée seconde pour contraintes de convexité aux nœuds
        spline_j_der2 = spline_j.derivative(2)
        con_array[:, j] = spline_j_der2(knots)
    
    # Calcul des coefficients pour chaque intervalle
    for i in range(kn):  # i indexe les intervalles entre nœuds
        t_k, t_k1 = knots[i], knots[i + 1]
        h = t_k1 - t_k
        
        # Pour chaque fonction de base active sur cet intervalle
        for l in range(degree + 1):
            j = i + l  # Index de la fonction de base
            
            if j < N:
                # Obtenir la B-spline et ses dérivées
                spline_j = List_Bsplines[j]
                
                # Pour contraintes de monotonie: dérivée première (cubique)
                spline_der1 = spline_j.derivative(1)
                
                # Évaluer la dérivée première à 4 points pour obtenir coefficients cubiques
                u_vals = np.linspace(0, 1, 4)
                x_vals = t_k + h * u_vals
                y_vals = spline_der1(x_vals)
                
                # Ajuster un polynôme cubique: a*u^3 + b*u^2 + c*u + d
                A = np.vstack([u_vals**3, u_vals**2, u_vals, np.ones(4)]).T
                coeffs_cubic, _, _, _ = np.linalg.lstsq(A, y_vals, rcond=None)
                DPi1[j, :, i] = coeffs_cubic
                
                # Pour contraintes de convexité: dérivée seconde (quadratique)
                spline_der2 = spline_j.derivative(2)
                
                # Évaluer la dérivée seconde à 3 points pour obtenir coefficients quadratiques
                u_vals_quad = np.linspace(0, 1, 3)
                x_vals_quad = t_k + h * u_vals_quad
                y_vals_quad = spline_der2(x_vals_quad)
                
                # Ajuster un polynôme quadratique: a*u^2 + b*u + c
                A_quad = np.vstack([u_vals_quad**2, u_vals_quad, np.ones(3)]).T
                coeffs_quad, _, _, _ = np.linalg.lstsq(A_quad, y_vals_quad, rcond=None)
                DPi2[j, :, i] = coeffs_quad
    
    return List_Bsplines, DPi1, DPi2, con_array

def apply_karlin_constraints_cubic(coeffs_cubic, monot_sign=1):
    """
    Applique les contraintes de Karlin pour un polynôme cubique (dérivée première)
    selon le théorème 3 de l'article
    """
    # coeffs_cubic = [a, b, c, d] pour a*u^3 + b*u^2 + c*u + d
    a, b, c, d = coeffs_cubic
    
    # Variables auxiliaires selon le théorème 3
    y0 = cp.Variable()
    y1 = cp.Variable()
    y2 = cp.Variable()
    x0 = cp.Variable()
    x1 = cp.Variable()
    x2 = cp.Variable()
    
    constraints = []
    
    # Équations (6a)-(6d) adaptées pour un polynôme cubique
    constraints.append(d == y0)
    constraints.append(c == 2*y1 + x0 - y0)
    constraints.append(b == y2 + 2*x1 - 2*y1)
    constraints.append(a == x2 - y2)
    
    # Contraintes de cône second-order (6e)-(6f)
    constraints.append(cp.SOC(x0 + x2, cp.hstack([x0 - x2, 2*x1])))
    constraints.append(cp.SOC(y0 + y2, cp.hstack([y0 - y2, 2*y1])))
    
    # Pour la monotonie: si monot_sign > 0, p(u) >= 0 pour u ∈ [0,1]
    # si monot_sign < 0, -p(u) >= 0 pour u ∈ [0,1]
    if monot_sign > 0:
        # p(0) = d >= 0 et p(1) = a + b + c + d >= 0
        constraints.append(d >= 0)
        constraints.append(a + b + c + d >= 0)
    elif monot_sign < 0:
        # -p(0) = -d >= 0 et -p(1) = -(a + b + c + d) >= 0
        constraints.append(-d >= 0)
        constraints.append(-(a + b + c + d) >= 0)
    
    return constraints

def apply_karlin_constraints_quadratic(coeffs_quad, cv_sign=1):
    """
    Applique les contraintes de Karlin pour un polynôme quadratique (dérivée seconde)
    selon la Proposition 5 du document PDF
    """
    # coeffs_quad = [a, b, c] pour a*u^2 + b*u + c
    # Note: selon la notation du PDF: p(x) = p0 + p1*x + p2*x^2
    # Donc: p0 = c, p1 = b, p2 = a
    p0, p1, p2 = coeffs_quad[2], coeffs_quad[1], coeffs_quad[0]
    
    constraints = []
    
    if cv_sign != 0:
        # Variable auxiliaire z0 >= 0
        z0 = cp.Variable()
        constraints.append(z0 >= 0)
        
        # Définir les transformations linéaires K1 et K2
        # K1: (p0, p1, p2, z0) -> (p0 + p2 + z0, p0 - p2 - z0)
        # K2: (p0, p1, p2, z0) -> p1 - z0
        
        # Pour cv_sign > 0: p(x) >= 0 sur [0,1]
        if cv_sign > 0:
            # Selon Proposition 5: (p0 + p2 + z0, p0 - p2 - z0, p1 - z0) ∈ Q3
            x1 = p0 + p2 + z0
            x2 = p0 - p2 - z0
            x3 = p1 - z0
            
            # Contrainte SOC: x1 >= ||(x2, x3)||_2
            constraints.append(cp.SOC(x1, cp.hstack([x2, x3])))
            
        # Pour cv_sign < 0: -p(x) >= 0 sur [0,1]
        elif cv_sign < 0:
            # Appliquer la même contrainte à -p(x)
            x1 = -p0 - p2 + z0  # -(p0 + p2) + z0
            x2 = -p0 + p2 - z0  # -(p0 - p2) - z0
            x3 = -p1 - z0       # -p1 - z0
            
            constraints.append(cp.SOC(x1, cp.hstack([x2, x3])))
            constraints.append(z0 >= 0)
    
    return constraints

def SplineQuantBspkn4(xtab, ytab, knots, tau, monot, cv, weight=None):
    """
    Régression quantile avec B-splines de degré 4 et contraintes de forme
    
    Parameters:
    -----------
    xtab, ytab : array-like
        Données x et y
    knots : int or list
        Nombre de nœuds ou liste des nœuds
    tau : float
        Paramètre quantile
    monot : int or list
        Contrainte de monotonie (+1 croissant, -1 décroissant, 0 aucune)
    cv : int or list  
        Contrainte de convexité (+1 convexe, -1 concave, 0 aucune)
    weight : array-like, optional
        Poids des observations
        
    Returns:
    --------
    polyn : BSpline object
        Fonction spline résultante
    """
    
    if weight is None:
        weight = np.ones(len(xtab))
    
    # Tri des données
    sort_idx = np.argsort(xtab)
    xtab = np.array(xtab)[sort_idx]
    ytab = np.array(ytab)[sort_idx]
    weight = np.array(weight)[sort_idx]
    
    n = len(xtab)
    
    if isinstance(knots, int):
        kn = knots - 1
        knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    kn = len(knots) - 1
    degree = 4
    N = kn + degree  # Nombre de fonctions de base
    
    # Gestion des contraintes
    if isinstance(monot, int):
        monot = monot * np.ones(kn)
    else:
        monot = np.array(monot)
    
    if np.isscalar(cv):
        cv_array = cv * np.ones(kn)  # Contraintes par intervalle
        cv_knots = cv * np.ones(kn + 1)  # Contraintes aux nœuds
    else:
        cv_array = np.array(cv)
        # Si cv a kn+1 éléments, les premiers kn sont pour les intervalles
        # et le dernier pour le nœud final
        if len(cv_array) == kn + 1:
            cv_knots = cv_array
            cv_array = cv_array[:-1]  # Les kn premiers pour les intervalles
        else:
            cv_array = cv_array
            cv_knots = np.concatenate([cv_array, [0]])
    
    print(f"Degré: {degree}, Nœuds: {kn}, Fonctions de base: {N}")
    print(f"Contraintes monotonie: {monot}")
    print(f"Contraintes convexité (intervalles): {cv_array}")
    print(f"Contraintes convexité (nœuds): {cv_knots}")
    
    # Construction des B-splines
    List_Bsplines, DPi1, DPi2, con_array = build_bsplines_and_deriv(knots, degree)
    
    # Matrice de design
    env_array = np.zeros((n, N))
    for j in range(N):
        env_array[:, j] = List_Bsplines[j](xtab)
    
    # Variables d'optimisation
    alpha = cp.Variable(N)  # Coefficients des B-splines
    
    # Fonction objectif
    residuals = ytab - env_array @ alpha
    objective = cp.Minimize(rhotau(residuals, tau))
    
    constraints = []
    
    # Contraintes de monotonie (sur la dérivée première, cubique)
    for i in range(kn):
        if monot[i] != 0:
            # Coefficients de la dérivée première sur cet intervalle
            # Somme pondérée des coefficients de base
            coeffs_sum = alpha @ DPi1[:, :, i]
            
            # Appliquer contraintes de Karlin pour cubique
            karlin_constraints = apply_karlin_constraints_cubic(coeffs_sum, monot[i])
            constraints.extend(karlin_constraints)
    
    # Contraintes de convexité (sur la dérivée seconde, quadratique)
    for i in range(kn):
        if cv_array[i] != 0:
            # Coefficients de la dérivée seconde sur cet intervalle
            coeffs_sum = alpha @ DPi2[:, :, i]
            
            # Appliquer contraintes de Karlin pour quadratique
            quad_constraints = apply_karlin_constraints_quadratic(coeffs_sum, cv_array[i])
            constraints.extend(quad_constraints)
    
    # Contraintes de convexité aux nœuds
    for j in range(kn + 1):
        if cv_knots[j] != 0:
            constraints.append(cv_knots[j] * (con_array[j, :] @ alpha) >= 0)
    
    # Résolution
    prob = cp.Problem(objective, constraints)
    
    # Essayer plusieurs solveurs
    try:
        prob.solve(verbose=True, solver='GUROBI')
    except:
        try:
            prob.solve(verbose=True, solver='ECOS')
        except:
            try:
                prob.solve(verbose=True, solver='SCS')
            except:
                prob.solve(verbose=True)
    
    if alpha.value is None:
        warnings.warn("L'optimisation n'a pas convergé")
        return None
    
    print(f"Statut: {prob.status}, Valeur objectif: {prob.value:.4f}")
    
    # Construction de la spline résultante
    l_end = np.min(knots)
    r_end = np.max(knots)
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    polyn = BSpline(s, alpha.value, degree)
    
    return polyn

def compute_derivatives(polyn, x_eval):
    """Calcule les dérivées première et seconde pour vérification"""
    # Dérivée première
    polyn_der1 = polyn.derivative(1)
    y_der1 = polyn_der1(x_eval)
    
    # Dérivée seconde  
    polyn_der2 = polyn.derivative(2)
    y_der2 = polyn_der2(x_eval)
    
    return y_der1, y_der2

def test_spline_degree4():
    """Test avec des splines de degré 4"""
    np.random.seed(42)
    n_points = 200
    
    # Données de test
    xtab = np.linspace(0, 1, n_points)
    ytab= (xtab*(1-xtab)*(xtab-0.3))+ 0.01*np.random.randn(n_points)
    #ytab = 2*xtab + 0.2*np.sin(10*np.pi*xtab) + 0.05*np.random.randn(n_points)
    
    kn = 8
    knots = kn + 1
    
    if isinstance(knots, int):
        kn = knots - 1
        knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    # Contraintes de test - différentes combinaisons
    test_cases = [
        {"name": "Monotone convexe", "monot": [1]*kn, "cv": [1]*kn},
        {"name": "Monotone concave", "monot": [1]*kn, "cv": [-1]*kn},
        {"name": "Non-monotone convexe", "monot": [0]*kn, "cv": [1]*kn},
        {"name": "Non-monotone concave", "monot": [0]*kn, "cv": [-1]*kn},
    ]
    
    fig, axes = plt.subplots(3, len(test_cases), figsize=(16, 12))
    
    for case_idx, test_case in enumerate(test_cases):
        tau = 0.1
        monot = test_case["monot"]
        cv = test_case["cv"]
        
        print(f"\n{'='*50}")
        print(f"Test: {test_case['name']}")
        print(f"{'='*50}")
        
        spline_result = SplineQuantBspkn4(xtab, ytab, knots, tau, monot, cv)
        
        if spline_result is not None:
            x_eval = np.linspace(0, 1, 500)
            y_eval = spline_result(x_eval)
            y_der1, y_der2 = compute_derivatives(spline_result, x_eval)
            
            # Plot 1: Fonction
            axes[0, case_idx].scatter(xtab, ytab, alpha=0.3, s=10, label='Données')
            axes[0, case_idx].plot(x_eval, y_eval, 'r-', linewidth=2, label=f'Spline (τ={tau})')
            axes[0, case_idx].plot(knots, np.ones_like(knots)*max(ytab), 'k|', markersize=10, label='Nœuds')
            axes[0, case_idx].set_xlabel('x')
            axes[0, case_idx].set_ylabel('y')
            axes[0, case_idx].set_title(f'{test_case["name"]}')
            axes[0, case_idx].legend(fontsize='small')
            axes[0, case_idx].grid(True, alpha=0.3)
            
            # Plot 2: Dérivée première
            axes[1, case_idx].plot(x_eval, y_der1, 'b-', linewidth=2, label="f'(x)")
            axes[1, case_idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            
            # Vérifier la monotonie
            if monot[0] != 0:
                if monot[0] > 0:
                    axes[1, case_idx].fill_between(x_eval, 0, y_der1, where=(y_der1>=0), 
                                                   alpha=0.3, color='green', label='f\' ≥ 0')
                    axes[1, case_idx].set_ylim(min(y_der1)*1.1, max(y_der1)*1.1)
                else:
                    axes[1, case_idx].fill_between(x_eval, y_der1, 0, where=(y_der1<=0), 
                                                   alpha=0.3, color='red', label='f\' ≤ 0')
                    axes[1, case_idx].set_ylim(min(y_der1)*1.1, max(y_der1)*1.1)
            
            axes[1, case_idx].set_xlabel('x')
            axes[1, case_idx].set_ylabel("f'(x)")
            axes[1, case_idx].set_title('Dérivée première')
            axes[1, case_idx].legend(fontsize='small')
            axes[1, case_idx].grid(True, alpha=0.3)
            
            # Plot 3: Dérivée seconde
            axes[2, case_idx].plot(x_eval, y_der2, 'g-', linewidth=2, label="f''(x)")
            axes[2, case_idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            
            # Vérifier la convexité
            if cv[0] != 0:
                if cv[0] > 0:
                    axes[2, case_idx].fill_between(x_eval, 0, y_der2, where=(y_der2>=0), 
                                                   alpha=0.3, color='green', label='f\'\' ≥ 0')
                else:
                    axes[2, case_idx].fill_between(x_eval, y_der2, 0, where=(y_der2<=0), 
                                                   alpha=0.3, color='red', label='f\'\' ≤ 0')
            
            axes[2, case_idx].set_xlabel('x')
            axes[2, case_idx].set_ylabel("f''(x)")
            axes[2, case_idx].set_title('Dérivée seconde')
            axes[2, case_idx].legend(fontsize='small')
            axes[2, case_idx].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return spline_result

def test2_spline_degree4():
    """Test avec des splines de degré 4"""
    np.random.seed(42)
    n_points = 200
    
    # Données de test
    xtab = np.linspace(0, 1, n_points)
    ytab = 2*xtab + 0.2*np.sin(10*np.pi*xtab) + 0.05*np.random.randn(n_points)
    
    kn = 10
    knots = kn + 1
    
    if isinstance(knots, int):
        kn = knots - 1
        knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    # Contraintes de test
    monot = [1] * (kn-5)+[0]*5  # Croissant partout
    
    cv = [0] * (kn-5)+[0]*5     # Pas de contrainte de convexité
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    for idx, tau in enumerate([0.1, 0.5, 0.9]):
        spline_result = SplineQuantBspkn4(xtab, ytab, knots, tau, monot, cv)
        
        if spline_result is not None:
            x_eval = np.linspace(0, 1, 500)
            y_eval = spline_result(x_eval)
            
            row = idx // 2
            col = idx % 2
            
            axes[row, col].scatter(xtab, ytab, alpha=0.3, s=10, label='Données')
            axes[row, col].plot(x_eval, y_eval, 'r-', linewidth=2, label=f'τ={tau}')
            axes[row, col].plot(knots, np.ones_like(knots)*max(ytab), 'k|', markersize=10, label='Nœuds')
            axes[row, col].set_xlabel('x')
            axes[row, col].set_ylabel('y')
            axes[row, col].set_title(f'Régression quantile (τ={tau}) avec B-spline degré 4')
            axes[row, col].legend()
            axes[row, col].grid(True, alpha=0.3)
    
    # Afficher les fonctions de base
    axes[1, 1].clear()
    degree = 4
    N = len(knots) - 1 + degree
    l_end = np.min(knots)
    r_end = np.max(knots)
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    x_eval = np.linspace(0, 1, 500)
    for j in range(min(N)):  # Afficher les 10 premières fonctions
        coefs = np.zeros(N)
        coefs[j] = 1
        basis_j = BSpline(s, coefs, degree)
        axes[1, 1].plot(x_eval, basis_j(x_eval), label=f'B{j}')
    
    axes[1, 1].set_title('Fonctions de base B-spline (degré 4)')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('Amplitude')
    axes[1, 1].legend(fontsize='small')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return spline_result


if __name__ == "__main__":
    test_spline_degree4()
