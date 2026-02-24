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
import matplotlib.pyplot as plt
import mosek

def rhotau(u, tau):
    """Fonction de perte quantile"""
    return cp.sum(cp.maximum(tau * u, (tau - 1) * u))

def build_bsplines_and_deriv_full(knots, degree=4):
    """
    Construit les B-splines de degré 4 et toutes leurs dérivées
    Retourne aussi les coefficients de la dérivée troisième (linéaire)
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
    # Pour les contraintes de degré 3 (dérivée troisième, linéaire)
    DPi3 = np.zeros((N, 2, kn))  # Pour dérivée troisième (linéaire: a*u + b)
    
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
                
                # Pour contraintes de degré 3: dérivée troisième (linéaire)
                spline_der3 = spline_j.derivative(3)
                
                # Évaluer la dérivée troisième à 2 points pour obtenir coefficients linéaires
                u_vals_lin = np.linspace(0, 1, 2)
                x_vals_lin = t_k + h * u_vals_lin
                y_vals_lin = spline_der3(x_vals_lin)
                
                # Ajuster un polynôme linéaire: a*u + b
                A_lin = np.vstack([u_vals_lin, np.ones(2)]).T
                coeffs_lin, _, _, _ = np.linalg.lstsq(A_lin, y_vals_lin, rcond=None)
                DPi3[j, :, i] = coeffs_lin  # [a, b] pour a*u + b
    
    return List_Bsplines, DPi1, DPi2, DPi3, con_array

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
            x1 = -p0 - p2 + z0
            x2 = -p0 + p2 - z0
            x3 = -p1 - z0
            
            constraints.append(cp.SOC(x1, cp.hstack([x2, x3])))
    
    return constraints

def apply_linear_constraints(coeffs_lin, d3_sign=1):
    """
    Applique les contraintes pour un polynôme linéaire (dérivée troisième)
    p(u) = a*u + b >= 0 (ou <= 0) pour u ∈ [0,1]
    """
    # coeffs_lin = [a, b] pour a*u + b
    a, b = coeffs_lin
    
    constraints = []
    
    if d3_sign > 0:
        # p(u) >= 0 sur [0,1] ↔ min(p(0), p(1)) >= 0
        constraints.append(b >= 0)      # p(0) >= 0
        constraints.append(a + b >= 0)  # p(1) >= 0
    elif d3_sign < 0:
        # p(u) <= 0 sur [0,1] ↔ max(p(0), p(1)) <= 0
        constraints.append(b <= 0)      # p(0) <= 0
        constraints.append(a + b <= 0)  # p(1) <= 0
    
    return constraints

def SplineQuantBspkn4_with_der3(xtab, ytab, knots, tau, monot, cv, d3=None, solver='GUROBI', weight=None):
    """
    Régression quantile avec B-splines de degré 4 et contraintes de forme
    incluant la dérivée troisième
    
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
    d3 : int or list, optional
        Contrainte sur la dérivée troisième (+1 croissante, -1 décroissante, 0 aucune)
    solver : str
        Solveur CVXPY à utiliser
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
        if len(cv_array) == kn + 1:
            cv_knots = cv_array
            cv_array = cv_array[:-1]
        else:
            cv_array = cv_array
            cv_knots = np.concatenate([cv_array, [0]])
    
    # Gestion des contraintes de dérivée troisième
    if d3 is None:
        d3 = np.zeros(kn)
    elif np.isscalar(d3):
        d3 = d3 * np.ones(kn)
    else:
        d3 = np.array(d3)
    
    print(f"Degré: {degree}, Nœuds: {kn}, Fonctions de base: {N}")
    print(f"Contraintes monotonie: {monot}")
    print(f"Contraintes convexité (intervalles): {cv_array}")
    print(f"Contraintes convexité (nœuds): {cv_knots}")
    print(f"Contraintes dérivée 3e: {d3}")
    
    # Construction des B-splines avec dérivée troisième
    List_Bsplines, DPi1, DPi2, DPi3, con_array = build_bsplines_and_deriv_full(knots, degree)
    
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
    
    # Contraintes sur la dérivée troisième (linéaire par morceau)
    for i in range(kn):
        if d3[i] != 0:
            # Coefficients de la dérivée troisième sur cet intervalle
            coeffs_sum = alpha @ DPi3[:, :, i]
            
            # Appliquer contraintes linéaires
            lin_constraints = apply_linear_constraints(coeffs_sum, d3[i])
            constraints.extend(lin_constraints)
    
    # Contraintes de convexité aux nœuds
    for j in range(kn + 1):
        if cv_knots[j] != 0:
            constraints.append(cv_knots[j] * (con_array[j, :] @ alpha) >= 0)
    
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

def compute_all_derivatives(polyn, x_eval):
    """Calcule toutes les dérivées jusqu'au 3ème ordre"""
    y = polyn(x_eval)
    y_der1 = polyn.derivative(1)(x_eval)
    y_der2 = polyn.derivative(2)(x_eval)
    y_der3 = polyn.derivative(3)(x_eval)
    return y, y_der1, y_der2, y_der3

def test_spline_degree4_with_der3():
    """Test avec des splines de degré 4 et contraintes sur la dérivée troisième"""
    np.random.seed(42)
    n_points = 200
    
    # Données de test
    xtab = np.linspace(0, 1, n_points)
    #ytab = 2*xtab + 0.2*np.sin(10*np.pi*xtab) + 0.05*np.random.randn(n_points)
    ytab=xtab**3-xtab**4
    kn = 8
    knots = kn + 1
    knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    # Différents tests
    test_cases = [
        {"name": "Sans contrainte", "monot": [0]*kn, "cv": [0]*kn, "d3": [0]*kn},
        {"name": "Monotone croissante", "monot": [1]*kn, "cv": [0]*kn, "d3": [0]*kn},
        {"name": "Convexe", "monot": [0]*kn, "cv": [1]*kn, "d3": [0]*kn},
        {"name": "Dérivée 3 positive", "monot": [0]*kn, "cv": [0]*kn, "d3": [1]*kn},
        {"name": "Dérivée 3 négative", "monot": [0]*kn, "cv": [0]*kn, "d3": [-1]*kn},
        {"name": "Toutes contraintes", "monot": [1]*kn, "cv": [1]*kn, "d3": [1]*kn},
    ]
    
    fig, axes = plt.subplots(4, len(test_cases), figsize=(20, 16))
    
    for case_idx, test_case in enumerate(test_cases):
        tau = 0.5
        print(f"\n{'='*50}")
        print(f"Test: {test_case['name']}")
        print(f"{'='*50}")
        
        spline_result = SplineQuantBspkn4_with_der3(
            xtab, ytab, knots, tau, 
            test_case["monot"], 
            test_case["cv"],
            test_case["d3"],
            solver='CLARABEL'
        )
        
        if spline_result is not None:
            x_eval = np.linspace(0, 1, 500)
            y, y1, y2, y3 = compute_all_derivatives(spline_result, x_eval)
            
            # Fonction
            axes[0, case_idx].scatter(xtab, ytab, alpha=0.3, s=10)
            axes[0, case_idx].plot(x_eval, y, 'b-', linewidth=2)
            axes[0, case_idx].plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=10)
            axes[0, case_idx].set_title(test_case['name'])
            axes[0, case_idx].grid(True, alpha=0.3)
            
            # Dérivée 1
            axes[1, case_idx].plot(x_eval, y1, 'g-', linewidth=2)
            axes[1, case_idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            axes[1, case_idx].grid(True, alpha=0.3)
            if case_idx == 0:
                axes[1, case_idx].set_ylabel("f'(x)")
            
            # Dérivée 2
            axes[2, case_idx].plot(x_eval, y2, 'r-', linewidth=2)
            axes[2, case_idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            axes[2, case_idx].grid(True, alpha=0.3)
            if case_idx == 0:
                axes[2, case_idx].set_ylabel("f''(x)")
            
            # Dérivée 3
            axes[3, case_idx].plot(x_eval, y3, 'm-', linewidth=2)
            axes[3, case_idx].axhline(y=0, color='k', linestyle='--', alpha=0.5)
            axes[3, case_idx].fill_between(x_eval, 0, y3, where=(y3>=0), alpha=0.3, color='green')
            axes[3, case_idx].fill_between(x_eval, y3, 0, where=(y3<=0), alpha=0.3, color='red')
            axes[3, case_idx].grid(True, alpha=0.3)
            if case_idx == 0:
                axes[3, case_idx].set_ylabel("f'''(x)")
            
            axes[3, case_idx].set_xlabel('x')
    
    plt.tight_layout()
    plt.show()
    
    return spline_result

if __name__ == "__main__":
    test_spline_degree4_with_der3()
