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

def build_bsplines_and_deriv(knots,degree):
    """Returns the derivatives of the Bspline according to the degree
    For linear and constant derivatives returns a vector of values at knots.
    When  derivatives are splines of degree 2 or 3, returns a matrix
    of normalised coefficients on the local power basis, i.e. the PP-form of the derivative of the spline
    prepared for applying Karlin-Studden constraints
    """
    kn=len(knots)-1
    N=kn+degree
    l_end = np.min(knots)
    r_end = np.max(knots)
    
    # Traitement des nœuds
    
    # Filtrage des données dans l'intervalle des nœuds
    #mask = (xtab >= np.min(knots))
    #xtab = xtab[mask]
    #ytab = ytab[mask]
    #weight = weight[mask]
    #n = len(xtab)
    
    # Extension de la séquence de nœuds
    s = np.concatenate([[l_end]*degree, knots, [r_end]*degree])
    
    List_Bsplines=[]
    #List_Bsplines_der2=[]
    
    if degree==4:
        Der1_array = np.zeros((N, degree ,kn ))   
        Der2_array = np.zeros((N,degree-1, kn))
        Der3_array = np.zeros((N,kn+1))

    if degree==3:
        Der1_array = np.zeros((N, degree ,kn ))   
        Der2_array = np.zeros((N,kn+1))
        Der3_array = np.zeros((N,kn))
        
    if degree==2:
        Der1_array = np.zeros((N,kn+1))   
        Der2_array = np.zeros((N,kn))

    if degree==1:
        Diff1=np.array([1,0])
        Der1_array = np.zeros((N,kn))   
        
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1        # j-th Bspline element of the basis 
        spline_j=BSpline(s, coefs,degree)
        List_Bsplines.append(spline_j)
        
        ppoly=PPoly.from_spline(List_Bsplines[j])
        ppoly_der1=ppoly.derivative(1)
        if degree>=2:
            ppoly_der2=ppoly.derivative(2)
            if degree>=3:
                ppoly_der3=ppoly.derivative(3)
        #coeffs_local=ppoly.c    
        for l in range(0,degree+1):
         nu=j-l+degree
         if ((nu>=degree) & (nu <N)):
            t_k, t_k1 = s[nu],s[nu+1]
            h = t_k1 - t_k
        

            #coeffs_local_der = ppoly_der1.c[:, nu]
            
            F1 = ppoly_der1.c[:, nu]
            
            if degree==4:           
                #Scale the local basis from (x-tk)^i to ((x-t_k)/(tk1_tk))^i,
                #First derivative
                B=np.array([h**3,h**2,h,1])
                Der1_array[j, :, nu-degree]= B*F1
                #Second derivative
                F2 = ppoly_der2.c[:, nu] 
                B2=[h**2,h,1]                
                Der2_array[j,:,nu-degree]=B2*F2[0:3]
            if degree==3:
                B2=np.array([h**2,h,1])
                Der1_array[j, :, nu-degree]=B2*F1 #+(t_k)*F2+(t_k)**2/2*F3)
        if degree==4:
            Der3_array[j, :]=ppoly_der3(knots)
        if degree==3:
            Der2_array[j, :]=ppoly_der2(knots)
            Der3_array[j, :]=ppoly_der3(knots[0:kn]) # coefficient de x^3       
        if degree==2:
            Der1_array[j, :]=ppoly_der1(knots)
            Der2_array[j, :]=ppoly_der2(knots[0:kn])
        if degree==1:
            Der1_array[j, :]=ppoly_der1(knots[0:kn])       
   
    if degree==1:
        return List_Bsplines, Der1_array
    if degree==2:
        return List_Bsplines, Der1_array, Der2_array
    if degree>=3:
        return List_Bsplines, Der1_array, Der2_array, Der3_array

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

def apply_karlin_constraints_quadratic(coeffs_quad, const_sign=1):
    """
    Applique les contraintes de Karlin pour un polynôme quadratique (dérivée seconde)
    selon la Proposition 5 du document PDF
    """
    # coeffs_quad = [a, b, c] pour a*u^2 + b*u + c
    # Note: selon la notation du PDF: p(x) = p0 + p1*x + p2*x^2
    # Donc: p0 = c, p1 = b, p2 = a
    p0, p1, p2 = coeffs_quad[2], coeffs_quad[1], coeffs_quad[0]
    #Karlin's matrices
    K1 = np.array([[1, 0, -1, -1],
                   [0, 1, 0, -1]])
    K2 = np.array([1, 0, 1, 1])
    
    constraints = []
    
    if const_sign != 0:
        # Variable auxiliaire z0 >= 0
        z0 = cp.Variable()
        constraints.append(z0 >= 0)
        P=const_sign*np.array([p0,p1,p2,const_sign*z0])
        x1 =    K2@P
        [x2,x3]=K1@P
            
        # Pour cv_sign > 0: p(x) >= 0 sur [0,1]
        #if cv_sign > 0:
            # Selon Proposition 5: (p0 + p2 + z0, p0 - p2 - z0, p1 - z0) ∈ Q3
            #x1 = p0 + p2 + z0
            #x2 = p0 - p2 - z0
            #x3 = p1 - z0
        #    P=cv_sign*np.array([p0,p1,p2,cv_sign*z0])
        #    x1 =    K2@P
        #    [x2,x3]=K1@P
            
            # Contrainte SOC: x1 >= ||(x2, x3)||_2
        constraints.append(cp.SOC(x1, cp.hstack([x2, x3])))  
    
    return constraints


def apply_val_constraints(const, sign=1):
    """
    Applique les contraintes pour un polynôme linéaire (dérivée troisième)
    p(u) = a*u + b >= 0 (ou <= 0) pour u ∈ [0,1]
    """
    constraints = []
    
    if sign > 0:
        # p(u) >= 0 sur [0,1] ↔ min(p(0), p(1)) >= 0
        constraints.append(const >= 0)      # p(0) >= 0
        #constraints.append(a + b >= 0)  # p(1) >= 0
    elif sign < 0:
        # p(u) <= 0 sur [0,1] ↔ max(p(0), p(1)) <= 0
        constraints.append(const <= 0)      # p(0) <= 0
        #constraints.append(a + b <= 0)  # p(1) <= 0
    
    return constraints

def SplineLinearQuant(xtab, ytab, knots, tau, monot=0, solver='GUROBI', weight=None):
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
    List_Bsplines, Der1_array = build_bsplines_and_deriv(knots,degree)
    
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
            deriv_sum = alpha @ Der1_array[:, i]
            
            # Appliquer contrainte sur la constante
            deriv_constraints = apply_val_constraints(deriv_sum, monot[i])
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

def SplineQuadraticQuant(xtab, ytab, knots, tau, monot=0, cv=0, solver='GUROBI', weight=None):
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
    List_Bsplines, Der1_array, Der2_array =  build_bsplines_and_deriv(knots,degree)
    
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
    for i in range(kn): #tous les intervallesq
        if monot[i] != 0:
            # Coefficients de la dérivée première sur cet intervalle
            coeffs_sum = alpha @ Der1_array[:,  i]  # [a, b] pour a*u + b
            
            # Appliquer contraintes linéaires
            lin_constraints = apply_val_constraints(coeffs_sum, monot[i])
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
    
    # Contraintes de convexité sur les itervalles (derivee seconde constante par morceaux
    for j in range(kn):
        if cv_knots[j] != 0:
            # La dérivée seconde au nœud j
            const_constraints = apply_val_constraints(Der2_array[:, j] @ alpha, cv_knots[j])
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

def SplineCubicQuant(xtab, ytab, knots, tau, monot=0, cv=0, der3=0, 
                                   solver='CLARABEL', weight=None):
    
    """
    Régression quantile avec splines cubiques
    Les contraintes de monotonie utilisent la même approche matricielle que les quartiques
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
    degree = 3
    N = kn + degree
    
    print(f"Degré: {degree}, Nœuds: {kn}, Fonctions de base: {N}")
    
    # Gestion des contraintes
    if isinstance(monot, int):
        monot = monot * np.ones(kn)
    else:
        monot = np.array(monot)
    print(f"Contraintes monotonie: {monot}")
    
    if np.isscalar(cv):
        cv = cv * np.ones(kn + 1)
    else:
        cv = np.array(cv)
    print(f"Contraintes convexité: {cv}")
    
    if isinstance(der3, int):
        der3 = der3 * np.ones(kn)
    else:
        der3 = np.array(der3)
    print(f"Contraintes dérivée 3: {der3}")
    
    # === CONSTRUCTION DES B-SPLINES ===
    # Construction des B-splines avec les dérivées
    l_end = np.min(knots)
    r_end = np.max(knots)
    s = np.concatenate([[l_end] * degree, knots, [r_end] * degree])
    
    List_Bsplines, Der1_array ,    Der2_array, Der3_array = build_bsplines_and_deriv(knots,degree)
    # === MATRICE DE DESIGN ===
    env_array = np.zeros((n, N))
    for j in range(N):
        env_array[:, j] = List_Bsplines[j](xtab)
    
    # === OPTIMISATION ===
    alpha = cp.Variable(N)
    z = cp.Variable(kn)  # Variables auxiliaires pour les contraintes
    
    residuals = ytab - env_array @ alpha
    objective = cp.Minimize(rhotau(residuals, tau))
    
    constraints = []
    
    # === CONTRAINTES DE MONOTONIE ===
    # Utilise la même approche matricielle que les quartiques
    for i in range(kn):
        if monot[i] != 0:
            # Coefficients de la dérivée première sur l'intervalle i
            # Der1_array[:, :, i] est (N, 3) : [a, b, c] pour a*u^2 + b*u + c
            coeffs_sum = alpha @ Der1_array[:, :, i]
            
            # Appliquer les contraintes de Karlin pour quadratique (comme pour les quartiques)
            quad_constraints = apply_karlin_constraints_quadratic(coeffs_sum, monot[i])
            constraints.extend(quad_constraints)
    
    # === CONTRAINTES DE CONVEXITÉ ===
    for j in range(kn + 1):
        if cv[j] != 0:
            coeffs_sum=Der2_array[:, j] @ alpha
            cv_constraints = apply_val_constraints(coeffs_sum, cv[j])
            constraints.extend(cv_constraints)
    
#            constraints.append(cv[j] * (Der2_array[:, j] @ alpha) >= 0)
    
    # === CONTRAINTES SUR DÉRIVÉE TROISIÈME ===
    for i in range(kn):
        if der3[i] != 0:
            coeffs_sum=Der3_array[:, i] @ alpha
            d3_constraints = apply_val_constraints(coeffs_sum, der3[i])
            constraints.extend(d3_constraints)
    
            #constraints.append(der3[i] * (Der3_array[:, i] @ alpha) >= 0)
    
    # === RÉSOLUTION ===
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False, solver=solver)
    
    if alpha.value is None:
        warnings.warn("L'optimisation n'a pas convergé")
        return None
    
    print(f"Statut: {prob.status}, Valeur objectif: {prob.value:.4f}")
    
    # Construction de la spline résultante
    polyn = BSpline(s, alpha.value, degree)
    
    return polyn
    
    # Der3_array: dérivée troisième (constante par intervalle)
    Der3_array = np.zeros((N, kn))
    
    # Création des B-splines de base
    for j in range(N):
        coefs = np.zeros(N)
        coefs[j] = 1
        spline_j = BSpline(s, coefs, degree)
        List_Bsplines.append(spline_j)
        
        # Convertir en PPoly pour obtenir les coefficients
        ppoly = PPoly.from_spline(spline_j)
        ppoly_der1 = ppoly.derivative(1)
        ppoly_der2 = ppoly.derivative(2)
        ppoly_der3 = ppoly.derivative(3)
        
        # Pour chaque intervalle
        for l in range(0, degree + 1):
            nu = j - l + degree
            if (nu >= degree) and (nu < N):
                t_k, t_k1 = s[nu], s[nu + 1]
                h = t_k1 - t_k
                
                # Coefficients de la dérivée première dans la base locale
                # ppoly_der1.c[:, nu] = [c0, c1, c2] pour c0 + c1*(x-tk) + c2*(x-tk)^2
                F1 = ppoly_der1.c[:, nu]
                
                # Normalisation en u = (x-tk)/h
                # f'(x) = c0 + c1*h*u + c2*h^2*u^2
                # Coefficients dans [u^2, u, 1] : [c2*h^2, c1*h, c0]
                Der1_array[j, :, nu - degree] = np.array([F1[2] * h**2, F1[1] * h, F1[0]])
        
        # Dérivée seconde aux nœuds (pour convexité)
        Der2_array[j, :] = ppoly_der2(knots)
        
        # Dérivée troisième aux nœuds (constante par intervalle)
        for i in range(kn):
            x_mid = (knots[i] + knots[i+1]) / 2
            Der3_array[j, i] = ppoly_der3(x_mid)
    
    # === MATRICE DE DESIGN ===
    env_array = np.zeros((n, N))
    for j in range(N):
        env_array[:, j] = List_Bsplines[j](xtab)
    
    # === OPTIMISATION ===
    alpha = cp.Variable(N)
    z = cp.Variable(kn)  # Variables auxiliaires pour les contraintes
    
    residuals = ytab - env_array @ alpha
    objective = cp.Minimize(rhotau(residuals, tau))
    
    constraints = []
    
    # === CONTRAINTES DE MONOTONIE ===
    # Utilise la même approche matricielle que les quartiques
    for i in range(kn):
        if monot[i] != 0:
            # Coefficients de la dérivée première sur l'intervalle i
            # Der1_array[:, :, i] est (N, 3) : [a, b, c] pour a*u^2 + b*u + c
            coeffs_sum = alpha @ Der1_array[:, :, i]
            
            # Appliquer les contraintes de Karlin pour quadratique (comme pour les quartiques)
            quad_constraints = apply_karlin_constraints_quadratic(coeffs_sum, monot[i])
            constraints.extend(quad_constraints)
    
    # === CONTRAINTES DE CONVEXITÉ ===
    for j in range(kn + 1):
        if cv[j] != 0:
            constraints.append(cv[j] * (Der2_array[:, j] @ alpha) >= 0)
    
    # === CONTRAINTES SUR DÉRIVÉE TROISIÈME ===
    for i in range(kn):
        if der3[i] != 0:
            constraints.append(der3[i] * (Der3_array[:, i] @ alpha) >= 0)
    
    # === RÉSOLUTION ===
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False, solver=solver)
    
    if alpha.value is None:
        warnings.warn("L'optimisation n'a pas convergé")
        return None
    
    print(f"Statut: {prob.status}, Valeur objectif: {prob.value:.4f}")
    
    # Construction de la spline résultante
    polyn = BSpline(s, alpha.value, degree)
    
    return polyn

def SplineQuarticQuant(xtab, ytab, knots, tau, monot, cv, d3=None, solver='GUROBI', weight=None):
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
        d3 = np.zeros(kn+1)
    elif np.isscalar(d3):
        d3 = d3 * np.ones(kn+1)
    elif len(d3)<(kn+1):
        d3 = np.array(d3+[0]*(kn+1-len(d3)))
    else:
        d3 = np.array(d3)
    
    print(f"Degré: {degree}, Nœuds: {kn}, Fonctions de base: {N}")
    print(f"Contraintes monotonie: {monot}")
    print(f"Contraintes convexité (intervalles): {cv_array}")
    print(f"Contraintes convexité (nœuds): {cv_knots}")
    print(f"Contraintes dérivée 3e: {d3}")
    
    # Construction des B-splines avec dérivée troisième
    List_Bsplines, Der1_array, Der2_array,Der3_array = build_bsplines_and_deriv(knots, degree)
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
            coeffs_sum = alpha @ Der1_array[:, :, i]
            
            # Appliquer contraintes de Karlin pour cubique
            karlin_constraints = apply_karlin_constraints_cubic(coeffs_sum, monot[i])
            constraints.extend(karlin_constraints)
    
    # Contraintes de convexité (sur la dérivée seconde, quadratique)
    for i in range(kn):
        if cv_array[i] != 0:
            # Coefficients de la dérivée seconde sur cet intervalle
            coeffs_sum2 = alpha @ Der2_array[:, :, i]
            
            # Appliquer contraintes de Karlin pour quadratique
            quad_constraints = apply_karlin_constraints_quadratic(coeffs_sum2, cv_array[i])
            constraints.extend(quad_constraints)
    
    # Contraintes sur la dérivée troisième (linéaire par morceau)
    for i in range(kn+1):
        if d3[i] != 0:
            # Coefficients de la dérivée troisième sur cet intervalle
            coeffs_sum = alpha @ Der3_array[:, i]
            
            # Appliquer contraintes linéaires à chaque noed
            #lin_constraints = apply_linear_constraints(coeffs_sum, d3[i])
            val3_constraints = apply_val_constraints(coeffs_sum, d3[i])
            constraints.extend(val3_constraints)
    
    # Contraintes de convexité aux nœuds
    #for j in range(kn + 1):
    #    if cv_knots[j] != 0:
    #        constraints.append(cv_knots[j] * (con_array[j, :] @ alpha) >= 0)
    
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



def compute_all_derivatives(polyn, x_eval, d):
    """Calcule toutes les dérivées jusqu'au 3ème ordre"""
    y = polyn(x_eval)
    y_der1 = polyn.derivative(1)(x_eval)
    if d>=2:
        y_der2 = polyn.derivative(2)(x_eval)
        if d>=3:
            y_der3 = polyn.derivative(3)(x_eval)
    if d==1:
        return y, y_der1
    if d==2:
        return y,y_der1, y_der2
    if d>=3:
        return y, y_der1, y_der2, y_der3

def test_spline_degree4_with_der3():
    """Test avec des splines de degré 4 et contraintes sur la dérivée troisième"""
    np.random.seed(42)
    n_points = 50
    
    # Données de test
    xtab = np.linspace(0, 1, n_points)
    ytab = 2*xtab + 0.2*np.sin(8*np.pi*xtab) + 0.05*np.random.randn(n_points)
    #ytab=xtab**3-xtab**4
    kn = 16
    knots = kn + 1
    knots = np.quantile(xtab, np.linspace(0, 1, kn + 1))
    
    # Différents tests
    test_cases = [
        {"name": "Sans contrainte", "monot": [0]*kn, "cv": [0]*(kn+1), "d3": [0]*kn},
        {"name": "Monotone croissante", "monot": [1]*kn, "cv": [0]*(kn+1), "d3": [0]*kn},
        {"name": "Convexe", "monot": [0]*kn, "cv": [1]*(kn+1), "d3": [0]*kn},
        {"name": "Dérivée 3 positive", "monot": [0]*kn, "cv": [0]*(kn+1), "d3": [1]*(kn+1)},
        {"name": "Dérivée 3 négative", "monot": [0]*kn, "cv": [0]*(kn+1), "d3": [-1]*(kn+1)},
        {"name": "Toutes contraintes", "monot": [1]*kn, "cv": [1]*(kn+1), "d3": [1]*(kn+1)},
    ]
    
    fig, axes = plt.subplots(4, len(test_cases), figsize=(20, 16))
    
    for case_idx, test_case in enumerate(test_cases):
        tau = 0.5
        print(f"\n{'='*50}")
        print(f"Test: {test_case['name']}")
        print(f"{'='*50}")
        
        spline_result1 = SplineQuarticQuant(
            xtab, ytab, knots, tau, 
            test_case["monot"], 
            test_case["cv"],
            test_case["d3"],
            solver='CLARABEL'
        )
        spline_result2 = SplineCubicQuant(
            xtab, ytab, knots, tau, 
            test_case["monot"], 
            test_case["cv"],
            test_case["d3"],
            solver='CLARABEL'
        )
        spline_result3 = SplineQuadraticQuant(
            xtab, ytab, knots, tau, 
            test_case["monot"], 
            test_case["cv"],
            solver='CLARABEL'
        )
        spline_result4 = SplineLinearQuant(
            xtab, ytab, knots, tau, 
            test_case["monot"], 
            solver='CLARABEL'
        )
        
        if spline_result2 is not None:
            x_eval = np.linspace(0, 1, 500)
            y, y1, y2, y3 = compute_all_derivatives(spline_result1, x_eval,4)
            yb, y1b, y2b, y3b = compute_all_derivatives(spline_result2, x_eval,3)
            yc, y1c,y2c = compute_all_derivatives(spline_result3, x_eval, 2)
            yd,y1d= compute_all_derivatives(spline_result4, x_eval, 1)
            # Fonction
            axes[0, case_idx].scatter(xtab, ytab, alpha=0.3, s=10)
            axes[0, case_idx].plot(x_eval, y, 'b-', linewidth=2)
            axes[0, case_idx].plot(x_eval, yb, 'b-', linewidth=2, color="green")
            axes[0, case_idx].plot(x_eval, yc, 'b-', linewidth=2, color="red")
            axes[0, case_idx].plot(x_eval, yd, 'b-', linewidth=2, color="black")
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

def test4():
    # Données d'exemple
    np.random.seed(42)
    n_points = 200
   
    #tau=0.5
    xtab = np.linspace(0, 1, n_points)
    #ytab=xtab*(1-xtab)
    #ytab=-xtab*(xtab-1)*(xtab-0.5)
    ytab=2*xtab + 0.2*np.sin(10*np.pi*xtab) + 0.05*np.random.randn(n_points)
    #ytab = 2*xtab+np.sin(10* np.pi * xtab)  +0.05 * np.random.randn(n_points)
    #ytab=np.exp((-5+10*xtab))/(1+np.exp((-5+10*xtab))) + 0.1 * np.random.randn(n_points) # logistique
    # concave

    kn=9
    knots=kn+1
    if isinstance(knots, int):
        # Si seul le nombre de nœuds est donné
        kn = knots-1
        # Calcul des quantiles pour les nœuds
        t_k_n = np.quantile(xtab, np.linspace(0, 1, kn + 1))


        knots=t_k_n.copy()
        fig, axes = plt.subplots(1, 2, figsize=(15, 4))
    # Appel de la fonction
    monot=[1]*kn
    cv=[1]*(kn+1)
    der3=0
#    for k in range(kn):
#        if k<4:
#            monot[k]=1
#        else:
#            cv[k]=0
    #monot=1      
    #cv=0
    #knots=[0.,  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ]
    solver='CLARABEL'
#    solver='Gurobi'
#    solver='Mosek'
#    solver='Cvxopt'
    #solver='Ecos'
    #solver='Scs'
    #monot=[1., 0, 0, 0., 0., 0., 0., 0., 0., 0.]
    for tau in [0.1,0.5,0.99]:
     spline_result = SplineQuantBspkn4_with_der3(xtab, ytab, knots, tau, monot, cv, der3, solver)
     # Évaluation
     x_eval = np.linspace(0, 1, 200)
     y_eval = spline_result(x_eval)
     print("Spline calculée avec succès")
     #Données et spline
     axes[0].plot(x_eval, y_eval, linewidth=1, label='tau='+str(tau))

    axes[0].bar(knots,max(ytab),width=0.02,alpha=0.5) 
    axes[0].scatter(xtab, ytab, alpha=0.5, s=20, label='Données')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_title('Régression quantile B-spline')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    N=kn+3
    for j in range(N):
            coef=np.zeros(N)
            coef[j]=1
            #polyn=BSpline(spline_result.t,coef,3)
            #axes[1].plot(x_eval,polyn(x_eval))
            #pp=PPoly.from_spline(polyn)
            #axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.5)

        
#        plt.plot(xtab,ytab,'r*')
#        plt.plot(x_eval,y_eval)
    plt.show()
    #print(spline_result.t,spline_result.c)
    return spline_result


if __name__ == "__main__":
    test_spline_degree4_with_der3()

