#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script de test pour les régressions quantile avec splines
degrés 1 à 4 avec contraintes de forme.
"""

import sys
import os

# Ajouter le chemin src au PYTHONPATH pour pouvoir importer le package
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
import matplotlib.pyplot as plt

# Import absolu depuis le package
from PysplineQuantReg import (
    SplineLinearQuant,
    SplineQuadraticQuant,
    SplineCubicQuant,
    SplineQuarticQuant
)

# ---------- Fonctions de génération de données (inchangées) ----------

def generate_test_data(n_points=200, noise=0.05, seed=42):
    """Génère des données de test avec différentes caractéristiques."""
    np.random.seed(seed)
    xtab = np.linspace(0, 1, n_points)
    ytab = 3 * xtab + 2 * np.sin(4 * np.pi * xtab) + noise * np.random.randn(n_points)
    return xtab, ytab

def generate_piecewise_data(n_points=200, noise=0.05, seed=42):
    """Génère des données avec des changements de pente."""
    np.random.seed(seed)
    xtab = np.linspace(0, 1, n_points)
    ytab = np.zeros(n_points)
    
    idx1 = xtab < 0.3
    ytab[idx1] = 2 * xtab[idx1]
    
    idx2 = (xtab >= 0.3) & (xtab < 0.6)
    ytab[idx2] = 0.6 + 0.2 * (xtab[idx2] - 0.3)
    
    idx3 = xtab >= 0.6
    ytab[idx3] = 0.7 + 3 * (xtab[idx3] - 0.6)
    
    ytab += noise * np.random.randn(n_points)
    return xtab, ytab

def generate_convex_data(n_points=200, noise=0.05, seed=42):
    """Génère des données convexes (parabole)."""
    np.random.seed(seed)
    xtab = np.linspace(0, 1, n_points)
    ytab = 5 * (xtab - 0.5)**2 + noise * np.random.randn(n_points)
    return xtab, ytab

def generate_concave_data(n_points=200, noise=0.05, seed=42):
    """Génère des données concaves."""
    np.random.seed(seed)
    xtab = np.linspace(0, 1, n_points)
    ytab = -5 * (xtab - 0.5)**2 + 2 + noise * np.random.randn(n_points)
    return xtab, ytab

# ---------- Fonctions de test ----------

def test_all_degrees():
    """Teste toutes les fonctions de régression avec différentes données."""
    
    print("=" * 70)
    print("TEST DES RÉGRESSIONS QUANTILE AVEC SPLINES")
    print("=" * 70)
    
    xtab, ytab = generate_test_data(n_points=200)
    kn = 12
    knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    tau = 0.5
    solver = 'CLARABEL'
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    x_eval = np.linspace(0, 1, 300)
    
    # === 1. Spline linéaire (degré 1) ===
    print("\n--- Spline linéaire (degré 1) ---")
    res_linear1 = SplineLinearQuant(xtab, ytab, knots, tau, monot=0, solver=solver)
    res_linear2 = SplineLinearQuant(xtab, ytab, knots, tau, monot=1, solver=solver)
    res_linear3 = SplineLinearQuant(xtab, ytab, knots, tau, monot=-1, solver=solver)
    
    ax = axes[0, 0]
    ax.scatter(xtab, ytab, alpha=0.3, s=10, color='gray', label='Données')
    if res_linear1 is not None:
        ax.plot(x_eval, res_linear1(x_eval), 'b-', linewidth=2, label='Sans contrainte')
    if res_linear2 is not None:
        ax.plot(x_eval, res_linear2(x_eval), 'g-', linewidth=2, label='Croissante')
    if res_linear3 is not None:
        ax.plot(x_eval, res_linear3(x_eval), 'r-', linewidth=2, label='Décroissante')
    ax.plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=10, label='Nœuds')
    ax.set_title('Spline linéaire (degré 1)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize='small')
    ax.grid(True, alpha=0.3)
    
    # === 2. Spline quadratique (degré 2) ===
    print("\n--- Spline quadratique (degré 2) ---")
    res_quad1 = SplineQuadraticQuant(xtab, ytab, knots, tau, monot=0, cv=0, solver=solver)
    res_quad2 = SplineQuadraticQuant(xtab, ytab, knots, tau, monot=1, cv=0, solver=solver)
    res_quad3 = SplineQuadraticQuant(xtab, ytab, knots, tau, monot=0, cv=1, solver=solver)
    
    ax = axes[0, 1]
    ax.scatter(xtab, ytab, alpha=0.3, s=10, color='gray', label='Données')
    if res_quad1 is not None:
        ax.plot(x_eval, res_quad1(x_eval), 'b-', linewidth=2, label='Sans contrainte')
    if res_quad2 is not None:
        ax.plot(x_eval, res_quad2(x_eval), 'g-', linewidth=2, label='Croissante')
    if res_quad3 is not None:
        ax.plot(x_eval, res_quad3(x_eval), 'orange', linewidth=2, label='Convexe')
    ax.plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=10, label='Nœuds')
    ax.set_title('Spline quadratique (degré 2)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize='small')
    ax.grid(True, alpha=0.3)
    
    # === 3. Spline cubique (degré 3) ===
    print("\n--- Spline cubique (degré 3) ---")
    res_cubic1 = SplineCubicQuant(xtab, ytab, knots, tau, monot=0, cv=0, der3=0, solver=solver)
    res_cubic2 = SplineCubicQuant(xtab, ytab, knots, tau, monot=1, cv=0, der3=0, solver=solver)
    res_cubic3 = SplineCubicQuant(xtab, ytab, knots, tau, monot=1, cv=1, der3=0, solver=solver)
    
    ax = axes[0, 2]
    ax.scatter(xtab, ytab, alpha=0.3, s=10, color='gray', label='Données')
    if res_cubic1 is not None:
        ax.plot(x_eval, res_cubic1(x_eval), 'b-', linewidth=2, label='Sans contrainte')
    if res_cubic2 is not None:
        ax.plot(x_eval, res_cubic2(x_eval), 'g-', linewidth=2, label='Croissante')
    if res_cubic3 is not None:
        ax.plot(x_eval, res_cubic3(x_eval), 'orange', linewidth=2, label='Croissante+Convexe')
    ax.plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=10, label='Nœuds')
    ax.set_title('Spline cubique (degré 3)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize='small')
    ax.grid(True, alpha=0.3)
    
    # === 4. Spline quartique (degré 4) ===
    print("\n--- Spline quartique (degré 4) ---")
    res_quart1 = SplineQuarticQuant(xtab, ytab, knots, tau, monot=0, cv=0, d3=0, solver=solver)
    res_quart2 = SplineQuarticQuant(xtab, ytab, knots, tau, monot=1, cv=0, d3=0, solver=solver)
    res_quart3 = SplineQuarticQuant(xtab, ytab, knots, tau, monot=0, cv=1, d3=0, solver=solver)
    
    ax = axes[1, 0]
    ax.scatter(xtab, ytab, alpha=0.3, s=10, color='gray', label='Données')
    if res_quart1 is not None:
        ax.plot(x_eval, res_quart1(x_eval), 'b-', linewidth=2, label='Sans contrainte')
    if res_quart2 is not None:
        ax.plot(x_eval, res_quart2(x_eval), 'g-', linewidth=2, label='Croissante')
    if res_quart3 is not None:
        ax.plot(x_eval, res_quart3(x_eval), 'orange', linewidth=2, label='Convexe')
    ax.plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=10, label='Nœuds')
    ax.set_title('Spline quartique (degré 4)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize='small')
    ax.grid(True, alpha=0.3)
    
    # === 5. Vérification des contraintes de convexité ===
    print("\n--- Vérification des contraintes ---")
    
    xtab_convex, ytab_convex = generate_convex_data(n_points=200)
    knots_convex = np.quantile(xtab_convex, np.linspace(0, 1, kn + 1))
    
    ax = axes[1, 2]
    ax.scatter(xtab_convex, ytab_convex, alpha=0.3, s=10, color='gray', label='Données convexes')
    
    res_none = SplineCubicQuant(xtab_convex, ytab_convex, knots_convex, tau, 
                                monot=0, cv=0, der3=0, solver=solver)
    res_convex = SplineCubicQuant(xtab_convex, ytab_convex, knots_convex, tau,
                                  monot=0, cv=1, der3=0, solver=solver)
    res_concave = SplineCubicQuant(xtab_convex, ytab_convex, knots_convex, tau,
                                   monot=0, cv=-1, der3=0, solver=solver)
    
    if res_none is not None:
        ax.plot(x_eval, res_none(x_eval), 'b-', linewidth=2, label='Sans contrainte')
    if res_convex is not None:
        ax.plot(x_eval, res_convex(x_eval), 'g-', linewidth=2, label='Convexe')
    if res_concave is not None:
        ax.plot(x_eval, res_concave(x_eval), 'r--', linewidth=2, label='Concave (forcée)')
    
    ax.plot(knots_convex, np.ones_like(knots_convex)*max(ytab_convex), 'r|', markersize=10, label='Nœuds')
    ax.set_title('Contraintes de convexité')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize='small')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print("\n" + "=" * 70)
    print("TESTS TERMINÉS")
    print("=" * 70)
    
    return {
        'linear': {'none': res_linear1, 'increasing': res_linear2, 'decreasing': res_linear3},
        'quadratic': {'none': res_quad1, 'increasing': res_quad2, 'convex': res_quad3},
        'cubic': {'none': res_cubic1, 'increasing': res_cubic2, 'increasing_convex': res_cubic3},
        'quartic': {'none': res_quart1, 'increasing': res_quart2, 'convex': res_quart3},
    }


def test_with_different_quantiles():
    """Test avec différents quantiles."""
    
    print("\n" + "=" * 70)
    print("TEST AVEC DIFFÉRENTS QUANTILES")
    print("=" * 70)
    
    xtab, ytab = generate_test_data(n_points=200)
    kn = 10
    knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    solver = 'CLARABEL'
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    x_eval = np.linspace(0, 1, 300)
    
    quantiles = [0.1, 0.5, 0.9]
    
    for idx, tau in enumerate(quantiles):
        ax = axes[idx]
        ax.scatter(xtab, ytab, alpha=0.2, s=10, color='gray')
        
        # Appel direct aux fonctions selon le degré
        for degree, color in zip([1, 2, 3, 4], ['blue', 'green', 'red', 'orange']):
            if degree == 1:
                res = SplineLinearQuant(xtab, ytab, knots, tau, monot=0, solver=solver)
            elif degree == 2:
                res = SplineQuadraticQuant(xtab, ytab, knots, tau, monot=0, cv=0, solver=solver)
            elif degree == 3:
                res = SplineCubicQuant(xtab, ytab, knots, tau, monot=0, cv=0, der3=0, solver=solver)
            else:  # degree == 4
                res = SplineQuarticQuant(xtab, ytab, knots, tau, monot=0, cv=0, d3=0, solver=solver)
            
            if res is not None:
                ax.plot(x_eval, res(x_eval), color=color, linewidth=2, label=f'deg {degree}')
        
        ax.plot(knots, np.ones_like(knots)*max(ytab), 'r|', markersize=8, label='Nœuds')
        ax.set_title(f'Quantile τ = {tau}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.legend(fontsize='small')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Lancer les tests
    results = test_all_degrees()
    test_with_different_quantiles()
