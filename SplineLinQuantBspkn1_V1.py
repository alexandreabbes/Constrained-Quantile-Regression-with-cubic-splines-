import numpy as np
import cvxpy as cp
from scipy.interpolate import BSpline, PPoly
import warnings

def rhotau(u, tau):
    """Fonction de perte quantile"""
    return cp.sum(cp.maximum(tau * u, (tau - 1) * u))

def build_linear_bsplines(knots):
    """
    Construit les B-splines linéaires (degré 1) et leurs dérivées
    
    Parameters:
    -----------
    knots : array
        Nœuds de la spline (knots[0], ..., knots[kn])
    
    Returns:
    --------
    List_Bsplines : list
        Liste des fonctions de base B-spline
    DPi1 : array (N, 2, kn)
        Coefficients des fonctions elles-mêmes (linéaires) sur chaque intervalle
        Utile pour les contraintes de monotonie (dérivée = constante)
    const_deriv : array (N, kn)
        Dérivée (constante) sur chaque intervalle
    """
    kn = len(knots) - 1
    degree = 1
    N = kn + degree  # Nombre de fonctions de base pour degré 1 = kn + 1
    
    l_end = np.min(knots)
    r_end = np.max(knots)
    
    # Extension de la séquence de nœuds pour les B-splines linéaires
    # Pour degré 1, on ajoute 1 nœud à chaque extrémité
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    List_Bsplines = []
    # Coefficients des fonctions (linéaires) sur chaque intervalle
    # Format: [a, b] pour a*u + b sur l'intervalle i (u dans [0,1])
    DPi1 = np.zeros((N, 2, kn))
    # Dérivée constante sur chaque intervalle
    const_deriv = np.zeros((N, kn))
    
    # Création des B-splines de base
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1
        spline_j = BSpline(s, coefs, degree)
        List_Bsplines.append(spline_j)
    
    # Calcul des coefficients sur chaque intervalle
    for i in range(kn):  # i indexe les intervalles entre nœuds
        t_k, t_k1 = knots[i], knots[i + 1]
        h = t_k1 - t_k
        
        # Pour une spline linéaire, 2 fonctions de base actives par intervalle
        for l in range(2):  # degree + 1 = 2
            j = i + l  # Index de la fonction de base
            
            if j < N:
                # Obtenir la B-spline
                spline_j = List_Bsplines[j]
                
                # La fonction est linéaire sur l'intervalle: a*u + b
                # Évaluer à 2 points pour obtenir les coefficients
                u_vals = np.array([0, 1])
                x_vals = t_k + h * u_vals
                y_vals = spline_j(x_vals)
                
                # Ajuster un polynôme linéaire: a*u + b
                # Système: [0, 1] * [a, b]^T = y(0)
                #          [1, 1] * [a, b]^T = y(1)
                A = np.vstack([u_vals, np.ones(2)]).T
                coeffs_lin, _, _, _ = np.linalg.lstsq(A, y_vals, rcond=None)
                DPi1[j, :, i] = coeffs_lin  # [a, b] pour a*u + b
                
                # La dérivée est constante = a/h (car u = (x - t_k)/h)
                # f'(x) = a/h
                const_deriv[j, i] = coeffs_lin[0] / h
    
    return List_Bsplines, DPi1, const_deriv

def apply_linear_function_constraints(coeffs_lin, xmin, xmax, sign=1):
    """
    Applique les contraintes de signe pour une fonction linéaire f(x) = a*x + b
    sur l'intervalle [xmin, xmax]
    
    Parameters:
    -----------
    coeffs_lin : array [a, b]
        Coefficients de la fonction linéaire dans la variable x originale
    xmin, xmax : float
        Bornes de l'intervalle
    sign : int
        +1 pour f(x) >= 0, -1 pour f(x) <= 0
    
    Returns:
    --------
    constraints : list
        Liste des contraintes CVXPY
    """
    a, b = coeffs_lin
    constraints = []
    
    if sign > 0:
        # f(x) >= 0 sur [xmin, xmax] ↔ min(f(xmin), f(xmax)) >= 0
        constraints.append(a * xmin + b >= 0)
        constraints.append(a * xmax + b >= 0)
    elif sign < 0:
        # f(x) <= 0 sur [xmin, xmax] ↔ max(f(xmin), f(xmax)) <= 0
        constraints.append(a * xmin + b <= 0)
        constraints.append(a * xmax + b <= 0)
    
    return constraints

def apply_constant_derivative_constraints(deriv, sign=1):
    """
    Applique les contraintes de monotonie via la dérivée constante
    
    Parameters:
    -----------
    deriv : float or cp.Variable
        Valeur de la dérivée (constante sur l'intervalle)
    sign : int
        +1 pour dérivée >= 0 (croissant), -1 pour dérivée <= 0 (décroissant)
    
    Returns:
    --------
    constraints : list
        Liste des contraintes CVXPY
    """
    constraints = []
    
    if sign > 0:
        constraints.append(deriv >= 0)
    elif sign < 0:
        constraints.append(deriv <= 0)
    
    return constraints

def SplineLinearQuantile(xtab, ytab, knots, tau, monot=0, solver='GUROBI', weight=None):
    """
    Régression quantile avec B-splines de degré 1 (affines par morceaux)
    et contraintes de monotonie
    
    Parameters:
    -----------
    xtab, ytab : array-like
        Données x et y
    knots : int or list
        Nombre de nœuds ou liste des nœuds
    tau : float
        Paramètre quantile (entre 0 et 1)
    monot : int or list
        Contrainte de monotonie (+1 croissant, -1 décroissant, 0 aucune)
        Pour une spline linéaire, la monotonie s'applique à la dérivée (constante)
    solver : str
        Solveur CVXPY à utiliser
    weight : array-like, optional
        Poids des observations
        
    Returns:
    --------
    polyn : BSpline object
        Fonction spline résultante (degré 1)
    """
    
    if weight is None:
        weight = np.ones(len(xtab))
    
    # Tri des données
    sort_idx = np.argsort(xtab)
    xtab = np.array(xtab)[sort_idx]
    ytab = np.array(ytab)[sort_idx]
    weight = np.array(weight)[sort_idx]
    
    n = len(xtab)
    
    # Gestion des nœuds
    if isinstance(knots, int):
        kn = knots - 1
        knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    kn = len(knots) - 1
    degree = 1
    N = kn + degree  # Nombre de fonctions de base = kn + 1
    
    # Gestion des contraintes de monotonie
    if isinstance(monot, int):
        monot = monot * np.ones(kn)  # Une contrainte par intervalle
    else:
        monot = np.array(monot)
    
    print(f"=== Régression quantile avec splines linéaires (degré 1) ===")
    print(f"Nœuds: {kn}, Fonctions de base: {N}")
    print(f"Contraintes monotonie (dérivée): {monot}")
    
    # Construction des B-splines linéaires
    List_Bsplines, DPi1, const_deriv = build_linear_bsplines(knots)
    
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
    
    # Contraintes de monotonie via la dérivée (constante sur chaque intervalle)
    for i in range(kn):
        if monot[i] != 0:
            # Dérivée sur l'intervalle i = somme pondérée des dérivées constantes
            deriv_sum = alpha @ const_deriv[:, i]
            
            # Appliquer contrainte sur la constante
            deriv_constraints = apply_constant_derivative_constraints(deriv_sum, monot[i])
            constraints.extend(deriv_constraints)
    
    # Alternative: contraindre directement la fonction (plus fort)
    # Mais les contraintes sur la dérivée sont plus naturelles pour la monotonie
    
    # Résolution
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False, solver=solver)
    
    if alpha.value is None:
        warnings.warn("L'optimisation n'a pas convergé, essayer un autre solveur")
        return None
    
    print(f"Statut: {prob.status}, Valeur objectif: {prob.value:.4f}")
    
    # Construction de la spline résultante
    l_end = np.min(knots)
    r_end = np.max(knots)
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    polyn = BSpline(s, alpha.value, degree)
    
    return polyn

def compute_linear_derivative(polyn, x_eval):
    """
    Calcule la dérivée d'une spline linéaire (constante par morceaux)
    """
    polyn_der = polyn.derivative(1)
    y_der = polyn_der(x_eval)
    return y_der

def test_linear_spline():
    """
    Fonction de test pour les splines linéaires
    """
    import matplotlib.pyplot as plt
    
    np.random.seed(42)
    n_points = 200
    
    # Données de test
    x = np.linspace(0, 1, n_points)
    y = 2*x + 0.2*np.sin(10*np.pi*x) + 0.05*np.random.randn(n_points)
    
    kn = 15  # Nombre d'intervalles
    knots = np.quantile(x, np.linspace(0, 1, kn + 1))
    
    # Test avec différentes contraintes
    test_cases = [
        {"name": "Sans contrainte", "monot": 0},
        {"name": "Croissante partout", "monot": 1},
        {"name": "Décroissante partout", "monot": -1},
        {"name": "Croissante puis décroissante", "monot": [1]*8 + [-1]*7},
        {"name": "Décroissante puis croissante", "monot": [-1]*7 + [1]*8},
        {"name": "Mixte", "monot": [1,1,1,0,0,-1,-1,-1,0,0,1,1,1,0,0]},
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    
    for idx, test in enumerate(test_cases[:6]):
        tau = 0.5
        spline = SplineLinearQuantile(x, y, knots, tau, 
                                     monot=test["monot"],
                                     solver='CLARABEL')
        
        if spline is not None:
            x_eval = np.linspace(0, 1, 500)
            y_eval = spline(x_eval)
            y_der = compute_linear_derivative(spline, x_eval)
            
            # Fonction
            axes[idx].scatter(x, y, alpha=0.3, s=10, label='Données')
            axes[idx].plot(x_eval, y_eval, 'b-', linewidth=2, label='Spline linéaire')
            axes[idx].plot(knots, np.ones_like(knots)*max(y), 'r|', markersize=10, label='Nœuds')
            axes[idx].set_title(test["name"])
            axes[idx].set_xlabel('x')
            axes[idx].set_ylabel('y')
            axes[idx].legend(fontsize='small')
            axes[idx].grid(True, alpha=0.3)
            
            # Ajouter les contraintes sur le graphique
            info_text = f"τ={tau}\n"
            if isinstance(test["monot"], list):
                info_text += f"M: {test['monot'][0]}...{test['monot'][-1]}"
            else:
                info_text += f"M={test['monot']}"
            axes[idx].text(0.05, 0.95, info_text, transform=axes[idx].transAxes,
                          verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()
    
    # Figure pour les dérivées
    fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8))
    axes2 = axes2.flatten()
    
    for idx, test in enumerate(test_cases[:6]):
        tau = 0.5
        spline = SplineLinearQuantile(x, y, knots, tau, 
                                     monot=test["monot"],
                                     solver='CLARABEL')
        
        if spline is not None:
            x_eval = np.linspace(0, 1, 500)
            y_der = compute_linear_derivative(spline, x_eval)
            
            # Dérivée
            axes2[idx].plot(x_eval, y_der, 'g-', linewidth=2, label="f'(x)")
            axes2[idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            
            # Colorier selon la contrainte
            if isinstance(test["monot"], int):
                if test["monot"] > 0:
                    axes2[idx].fill_between(x_eval, 0, y_der, where=(y_der>=0), 
                                           alpha=0.3, color='green', label='f\' ≥ 0')
                elif test["monot"] < 0:
                    axes2[idx].fill_between(x_eval, y_der, 0, where=(y_der<=0), 
                                           alpha=0.3, color='red', label='f\' ≤ 0')
            
            axes2[idx].set_title(f"{test['name']} - Dérivée")
            axes2[idx].set_xlabel('x')
            axes2[idx].set_ylabel("f'(x)")
            axes2[idx].legend(fontsize='small')
            axes2[idx].grid(True, alpha=0.3)
            
            # Ajouter les nœuds
            for k in knots:
                axes2[idx].axvline(x=k, color='gray', linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.show()
    
    return spline

def compare_linear_with_higher_degree():
    """
    Compare les splines linéaires avec les splines de degré supérieur
    """
    import matplotlib.pyplot as plt
    
    np.random.seed(42)
    n_points = 200
    
    # Données de test avec différentes caractéristiques
    x = np.linspace(0, 1, n_points)
    
    # Données avec changements de pente brutaux
    y = np.zeros_like(x)
    y[x < 0.3] = 1.5 * x[x < 0.3]
    y[(x >= 0.3) & (x < 0.6)] = 0.45 + 0.5 * (x[(x >= 0.3) & (x < 0.6)] - 0.3)
    y[x >= 0.6] = 0.6 + 2 * (x[x >= 0.6] - 0.6)
    y += 0.03 * np.random.randn(n_points)
    
    kn = 15
    knots = np.quantile(x, np.linspace(0, 1, kn + 1))
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    
    # Essayer différents τ
    taus = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
    
    for idx, tau in enumerate(taus):
        row, col = idx // 3, idx % 3
        
        spline_linear = SplineLinearQuantile(x, y, knots, tau, monot=0, solver='CLARABEL')
        
        if spline_linear is not None:
            x_eval = np.linspace(0, 1, 500)
            y_linear = spline_linear(x_eval)
            
            axes[row, col].scatter(x, y, alpha=0.3, s=10, label='Données')
            axes[row, col].plot(x_eval, y_linear, 'b-', linewidth=2, label=f'Linéaire τ={tau}')
            axes[row, col].plot(knots, np.ones_like(knots)*max(y), 'r|', markersize=8, label='Nœuds')
            axes[row, col].set_title(f'Régression quantile τ={tau}')
            axes[row, col].set_xlabel('x')
            axes[row, col].set_ylabel('y')
            axes[row, col].legend(fontsize='small')
            axes[row, col].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Visualisation des fonctions de base
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 5))
    
    degree = 1
    N = len(knots) - 1 + degree
    l_end = np.min(knots)
    r_end = np.max(knots)
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    x_eval = np.linspace(0, 1, 500)
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1
        basis_j = BSpline(s, coefs, degree)
        ax2.plot(x_eval, basis_j(x_eval), label=f'B{j}')
    
    ax2.set_title('Fonctions de base B-spline linéaires (degré 1)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('Amplitude')
    ax2.legend(fontsize='small', ncol=3)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return spline_linear

if __name__ == "__main__":
    print("Test des splines linéaires avec différentes contraintes")
    test_linear_spline()
    
    print("\nComparaison avec différents quantiles")
    compare_linear_with_higher_degree()
