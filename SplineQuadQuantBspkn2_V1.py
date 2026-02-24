#!/usr/bin/env python

__author__ = "Alexandre Abbes"
__copyright__ = "Copyright 2026, Alexandre Adel Abbes"
__credits__ = ["Alexandre Abbes"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Alexandre Abbes"
__email__ = "alexandre.abbes@proton.me"


import numpy as np
import cvxpy as cp
from scipy.interpolate import BSpline, PPoly
import warnings

def rhotau(u, tau):
    """Fonction de perte quantile"""
    return cp.sum(cp.maximum(tau * u, (tau - 1) * u))

def build_quadratic_bsplines(knots):
    """
    Construit les B-splines quadratiques (degré 2) et leurs dérivées
    
    Parameters:
    -----------
    knots : array
        Nœuds de la spline (knots[0], ..., knots[kn])
    
    Returns:
    --------
    List_Bsplines : list
        Liste des fonctions de base B-spline
    DPi1 : array (N, 2, kn)
        Coefficients des dérivées premières (linéaires) sur chaque intervalle
        DPi1[j, :, i] = [a, b] pour a*u + b sur l'intervalle i
    con_array : array (kn+1, N)
        Valeurs des dérivées secondes aux nœuds pour contraintes de convexité
    """
    kn = len(knots) - 1
    degree = 2
    N = kn + degree  # Nombre de fonctions de base
    
    l_end = np.min(knots)
    r_end = np.max(knots)
    
    # Extension de la séquence de nœuds pour les B-splines quadratiques
    # Pour degré 2, on ajoute 2 nœuds à chaque extrémité
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    List_Bsplines = []
    # Pour les contraintes de monotonie (dérivée première, linéaire)
    DPi1 = np.zeros((N, 2, kn))  # [a, b] pour a*u + b
    # Pour les contraintes de convexité aux nœuds (dérivée seconde, constante)
    con_array = np.zeros((kn + 1, N))
    
    # Création des B-splines de base
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1
        spline_j = BSpline(s, coefs, degree)
        List_Bsplines.append(spline_j)
        
        # Dérivée seconde (constante) pour contraintes de convexité aux nœuds
        spline_j_der2 = spline_j.derivative(2)
        con_array[:, j] = spline_j_der2(knots)
    
    # Calcul des coefficients des dérivées premières sur chaque intervalle
    for i in range(kn):  # i indexe les intervalles entre nœuds
        t_k, t_k1 = knots[i], knots[i + 1]
        h = t_k1 - t_k
        
        # Pour une spline quadratique, 3 fonctions de base actives par intervalle
        for l in range(3):  # degree + 1 = 3
            j = i + l  # Index de la fonction de base
            
            if j < N:
                # Obtenir la B-spline et sa dérivée première
                spline_j = List_Bsplines[j]
                spline_der1 = spline_j.derivative(1)
                
                # La dérivée première d'une spline quadratique est linéaire
                # Évaluer à 2 points pour obtenir coefficients linéaires: a*u + b
                u_vals = np.array([0, 1])
                x_vals = t_k + h * u_vals
                y_vals = spline_der1(x_vals)
                
                # Ajuster un polynôme linéaire: a*u + b
                # Système: [0, 1] * [a, b]^T = y(0)
                #          [1, 1] * [a, b]^T = y(1)
                A = np.vstack([u_vals, np.ones(2)]).T
                coeffs_lin, _, _, _ = np.linalg.lstsq(A, y_vals, rcond=None)
                DPi1[j, :, i] = coeffs_lin  # [a, b] pour a*u + b
    
    return List_Bsplines, DPi1, con_array

def apply_linear_constraints(coeffs_lin, monot_sign=1):
    """
    Applique les contraintes pour un polynôme linéaire (dérivée première)
    p(u) = a*u + b >= 0 (ou <= 0) pour u ∈ [0,1]
    
    Parameters:
    -----------
    coeffs_lin : array [a, b]
        Coefficients du polynôme linéaire
    monot_sign : int
        +1 pour p(u) >= 0, -1 pour p(u) <= 0
    
    Returns:
    --------
    constraints : list
        Liste des contraintes CVXPY
    """
    a, b = coeffs_lin
    constraints = []
    
    if monot_sign > 0:
        # p(u) >= 0 sur [0,1] ↔ min(p(0), p(1)) >= 0
        constraints.append(b >= 0)      # p(0) >= 0
        constraints.append(a + b >= 0)  # p(1) >= 0
    elif monot_sign < 0:
        # p(u) <= 0 sur [0,1] ↔ max(p(0), p(1)) <= 0
        constraints.append(b <= 0)      # p(0) <= 0
        constraints.append(a + b <= 0)  # p(1) <= 0
    
    return constraints

def apply_constant_constraints(constant, cv_sign=1):
    """
    Applique les contraintes pour une constante (dérivée seconde)
    
    Parameters:
    -----------
    constant : float or cp.Variable
        Valeur constante
    cv_sign : int
        +1 pour constante >= 0, -1 pour constante <= 0
    
    Returns:
    --------
    constraints : list
        Liste des contraintes CVXPY
    """
    constraints = []
    
    if cv_sign > 0:
        constraints.append(constant >= 0)
    elif cv_sign < 0:
        constraints.append(constant <= 0)
    
    return constraints

def SplineQuadraticQuantile(xtab, ytab, knots, tau, monot=0, cv=0, solver='GUROBI', weight=None):
    """
    Régression quantile avec B-splines de degré 2 et contraintes de forme
    
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
        Pour une spline quadratique, la monotonie s'applique à la dérivée première
    cv : int or list  
        Contrainte de convexité (+1 convexe, -1 concave, 0 aucune)
        Pour une spline quadratique, la convexité s'applique à la dérivée seconde
    solver : str
        Solveur CVXPY à utiliser
    weight : array-like, optional
        Poids des observations
        
    Returns:
    --------
    polyn : BSpline object
        Fonction spline résultante (degré 2)
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
    degree = 2
    N = kn + degree  # Nombre de fonctions de base
    
    # Gestion des contraintes
    if isinstance(monot, int):
        monot = monot * np.ones(kn)
    else:
        monot = np.array(monot)
    
    if np.isscalar(cv):
        cv_array = cv * np.ones(kn)      # Contraintes par intervalle pour la convexité
        cv_knots = cv * np.ones(kn + 1)  # Contraintes aux nœuds
    else:
        cv_array = np.array(cv)
        if len(cv_array) == kn + 1:
            cv_knots = cv_array
            cv_array = cv_array[:-1]
        else:
            cv_array = cv_array
            cv_knots = np.concatenate([cv_array, [0]])
    
    print(f"=== Régression quantile avec splines quadratiques (degré 2) ===")
    print(f"Nœuds: {kn}, Fonctions de base: {N}")
    print(f"Contraintes monotonie (dérivée 1): {monot}")
    print(f"Contraintes convexité (dérivée 2) - intervalles: {cv_array}")
    print(f"Contraintes convexité (dérivée 2) - nœuds: {cv_knots}")
    
    # Construction des B-splines quadratiques
    List_Bsplines, DPi1, con_array = build_quadratic_bsplines(knots)
    
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
    
    # Contraintes de monotonie (sur la dérivée première, linéaire)
    for i in range(kn):
        if monot[i] != 0:
            # Coefficients de la dérivée première sur cet intervalle
            coeffs_sum = alpha @ DPi1[:, :, i]  # [a, b] pour a*u + b
            
            # Appliquer contraintes linéaires
            lin_constraints = apply_linear_constraints(coeffs_sum, monot[i])
            constraints.extend(lin_constraints)
    
    # Contraintes de convexité (sur la dérivée seconde, constante)
    for i in range(kn):
        if cv_array[i] != 0:
            # La dérivée seconde est constante = 2*a (où a est le coeff de u^2)
            # Mais on peut utiliser directement les valeurs aux nœuds
            # Pour chaque intervalle, on peut aussi contraindre les dérivées secondes aux nœuds
            
            # Option 1: Contraindre aux nœuds (plus simple)
            # On le fait dans la boucle suivante
            pass
    
    # Contraintes de convexité aux nœuds
    for j in range(kn + 1):
        if cv_knots[j] != 0:
            # La dérivée seconde au nœud j
            const_constraints = apply_constant_constraints(con_array[j, :] @ alpha, cv_knots[j])
            constraints.extend(const_constraints)
    
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

def compute_quadratic_derivatives(polyn, x_eval):
    """
    Calcule les dérivées première et seconde d'une spline quadratique
    """
    # Dérivée première (linéaire)
    polyn_der1 = polyn.derivative(1)
    y_der1 = polyn_der1(x_eval)
    
    # Dérivée seconde (constante par morceaux)
    polyn_der2 = polyn.derivative(2)
    y_der2 = polyn_der2(x_eval)
    
    return y_der1, y_der2

def test_quadratic_spline():
    """
    Fonction de test pour les splines quadratiques
    """
    import matplotlib.pyplot as plt
    
    np.random.seed(42)
    n_points = 200
    
    # Données de test
    x = np.linspace(0, 1, n_points)
    y = 2*x + 0.2*np.sin(10*np.pi*x) + 0.05*np.random.randn(n_points)
    
    kn = 10
    knots = np.quantile(x, np.linspace(0, 1, kn + 1))
    
    # Test avec différentes contraintes
    test_cases = [
        {"name": "Sans contrainte", "monot": 0, "cv": 0},
        {"name": "Croissante", "monot": 1, "cv": 0},
        {"name": "Décroissante", "monot": -1, "cv": 0},
        {"name": "Convexe", "monot": 0, "cv": 1},
        {"name": "Concave", "monot": 0, "cv": -1},
        {"name": "Croissante & Convexe", "monot": 1, "cv": 1},
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    
    for idx, test in enumerate(test_cases):
        tau = 0.5
        spline = SplineQuadraticQuantile(x, y, knots, tau, 
                                        monot=test["monot"], 
                                        cv=test["cv"],
                                        solver='CLARABEL')
        
        if spline is not None:
            x_eval = np.linspace(0, 1, 500)
            y_eval = spline(x_eval)
            y_der1, y_der2 = compute_quadratic_derivatives(spline, x_eval)
            
            # Fonction
            axes[idx].scatter(x, y, alpha=0.3, s=10, label='Données')
            axes[idx].plot(x_eval, y_eval, 'b-', linewidth=2, label='Spline')
            axes[idx].plot(knots, np.ones_like(knots)*max(y), 'r|', markersize=10, label='Nœuds')
            axes[idx].set_title(test["name"])
            axes[idx].set_xlabel('x')
            axes[idx].set_ylabel('y')
            axes[idx].legend(fontsize='small')
            axes[idx].grid(True, alpha=0.3)
            
            # Petit texte pour les contraintes
            info_text = f"τ={tau}\nM={test['monot']} C={test['cv']}"
            axes[idx].text(0.05, 0.95, info_text, transform=axes[idx].transAxes,
                          verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()
    
    # Figure supplémentaire pour les dérivées
    fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8))
    axes2 = axes2.flatten()
    
    for idx, test in enumerate(test_cases[:6]):
        tau = 0.5
        spline = SplineQuadraticQuantile(x, y, knots, tau, 
                                        monot=test["monot"], 
                                        cv=test["cv"],
                                        solver='CLARABEL')
        
        if spline is not None:
            x_eval = np.linspace(0, 1, 500)
            y_der1, y_der2 = compute_quadratic_derivatives(spline, x_eval)
            
            # Dérivée première
            axes2[idx].plot(x_eval, y_der1, 'g-', linewidth=2, label="f'(x)")
            axes2[idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            
            # Colorier selon la contrainte
            if test["monot"] > 0:
                axes2[idx].fill_between(x_eval, 0, y_der1, where=(y_der1>=0), 
                                       alpha=0.3, color='green', label='f\' ≥ 0')
            elif test["monot"] < 0:
                axes2[idx].fill_between(x_eval, y_der1, 0, where=(y_der1<=0), 
                                       alpha=0.3, color='red', label='f\' ≤ 0')
            
            axes2[idx].set_title(f"{test['name']} - Dérivée première")
            axes2[idx].set_xlabel('x')
            axes2[idx].set_ylabel("f'(x)")
            axes2[idx].legend(fontsize='small')
            axes2[idx].grid(True, alpha=0.3)
            
            # Dérivée seconde sur la même figure en petit
            ax2 = axes2[idx].twinx()
            ax2.plot(x_eval, y_der2, 'm-', linewidth=1, alpha=0.7, label="f''(x)")
            ax2.set_ylabel("f''(x)", color='m')
            ax2.tick_params(axis='y', labelcolor='m')
    
    plt.tight_layout()
    plt.show()
    
    return spline

if __name__ == "__main__":
    test_quadratic_spline()
