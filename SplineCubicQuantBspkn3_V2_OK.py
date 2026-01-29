import numpy as np
import cvxpy as cp
import cvxpy_gurobi
from scipy.interpolate import BSpline, splev, PPoly
from scipy.stats import percentileofscore
from matplotlib import pyplot as plt
import warnings
#import mosek #only for testing, 1year free academic licence. No significant imrovement at first sight.





#Ca marche a peu près.
#27/12/2025: il semble que les contraintes monotones sont appliquées avec un peu d'anticipation.
#01/01/2026, ça marche bien, mais sur les 3-4 premiers inter_knots de l'intervalle, les contraintes ne sont pas bien appliquées.
#02/01 : le problème est dû à la multiplication par t_k qui vaut 0.


def rhotau(u, tau):
    """Fonction de perte quantile équivalente à celle du MATLAB"""
    return cp.sum(cp.maximum(tau * u, (tau - 1) * u))

def build_bsplines_and_deriv(knots,degree=3):
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
    s = np.concatenate([[l_end]*3, knots, [r_end]*3])
    
    List_Bsplines=[]
    List_Bsplines_der2=[]

    DPi2 = np.zeros((N, 3,kn )) 
    
    con_array = np.zeros((kn+1, N))
    Der3_array= np.zeros(kn,N) # one value per inter-knot
    for j in range(N):  
        coefs = np.zeros(N)
        coefs[j] = 1
        
        spline_j=BSpline(s, coefs,degree)
        List_Bsplines.append(spline_j)

        spline_j_der2=spline_j.derivative(2) 
        spline_j_der3=pp
        con_array[:, j] = spline_j_der2(knots)
        
    for nu in range(3,N):
        
      t_k, t_k1 = s[nu],s[nu+1]
      print(nu,s[nu], s[nu+1])
      h = t_k1 - t_k
      
      for l in range(0,4):  
        j=nu+l-3
        print(nu,j,l)
        ppoly=PPoly.from_spline(List_Bsplines[j])
        
        ppoly_der=ppoly.derivative(1)
        
                
                # Get coefficients for the derivative on this piece
                # The derivative is already a quadratic: c[0]*(x-x0)^2 + c[1]*(x-x0) + c[2]
        coeffs_local_der = ppoly_der.c[:, nu]
               
         
        coeffs_local=ppoly.c[:,nu]
        
#        if j<4 and nu<4: #for debugging purpose
         #print('ppoly_der(j)',nu,j, ppoly_der[0:6,:])
         #print('ppoly',j,ppoly[:,0:nu])
#         print("coefs",j,nu)
#         print(coeffs_local,coeffs_local_der)
#         print('h,tk',h,t_k1)
        DPi2[j, :, nu-3] =[3*coeffs_local[0]*h**2,  # a pour 3u^2
                          2*coeffs_local[1]*h,        # b pour 2u
                          coeffs_local[2] ]            # c
        if j<4 and nu<4:
            print(nu,j)
            print(DPi2[j,:,nu-3]) # test pour verifier les  matrices de contraintes. 
        der3_array=Dpi2[:,0,:]/3
    return List_Bsplines, DPi2, con_array, der3_array


def SplineCubicQuantBspkn3(xtab, ytab, knots, tau, monot, cv, weight=None):
    #print('len',len(xtab))
    """
    Version Python de la fonction MATLAB SplineCubicQuantBspkn3
    
    Parameters:
    -----------
    xtab, ytab : array-like
        Données x et y
    knots : int or list
        Nombre de nœuds ou liste des nœuds
    tau : float
        Paramètre quantile
    monot : int or list
        Contrainte de monotonie (+1/-1/0)
    cv : int or list  
        Contrainte de convexité (+1/-1/0)
    weight : array-like, optional
        Poids des observations
        
    Returns:
    --------
    polyn : BSpline object
        Fonction spline résultante
    """
    
    # Initialisation des poids si non fournis

    if weight is None:
        weight = np.ones(len(xtab))
    
    # Tri des données
    sort_idx = np.argsort(xtab)
    xtab = np.array(xtab)[sort_idx]
    ytab = np.array(ytab)[sort_idx]
    weight = np.array(weight)[sort_idx]

    
    n = len(xtab)

    if isinstance(knots, int):
        # Si seul le nombre de nœuds est donné
        kn = knots-1
    
        # Calcul des quantiles pour les nœuds
        t_k_n = np.quantile(xtab, np.linspace(0, 1, kn + 1))


        knots=t_k_n.copy()
    kn=len(knots)-1
    
    # Gestion des contraintes de monotonie et convexité
    # Paramètres B-spline
    m = 3  # degré cubique
    N = kn + m  # nombre de fonctions de base B-spline
    print("knots:",knots)
    if isinstance(monot,int):
        monot = monot * np.ones(kn)
    else:
        monot = np.array(monot)
    print('Monotoncity contstraints',monot)    
        
    if np.isscalar(cv):
        cv = cv * np.ones(kn + 1)
    else:
        cv = np.array(cv)
    print("convexity constraints",cv)
    
    
    
    # Construction de la matrice de design B-spline
    # Matrice de design pour l'évaluation
    List_Bsplines, DPi2, con_array=build_bsplines_and_deriv(knots,degree=3)
  
    env_array = np.zeros((n, N))
    for j in range(N):
        spline_j=List_Bsplines[j]
        env_array[:, j] = spline_j(xtab)


    # Matrices de Karlin pour les contraintes de monotonie 
    

    K1 = np.array([[1, 0, -1, -1],
                   [0, 1, 0, -1]])
    K2 = np.array([1, 0, 1, 1])
    
    A_j_array = np.zeros((2, N + kn, kn))
    c_j_array = np.zeros((kn, N + kn))
    
    for j in range(kn):
        Ej = np.zeros(N+kn)
        Ej[N+j] = 1

        DPi2_j = DPi2[:, :, j].T  # (3, N) coeff des derivees Bspline normalisee deg=3,2,1
        # sur le j_ieme intervalle
        DPi2_j_flipped = np.flipud(DPi2_j)
       
        top = np.concatenate([DPi2_j_flipped, np.zeros((3, kn))], axis=1)

        bottom=Ej.reshape((1,N+kn))
        DPj2 = np.concatenate([top, bottom], axis=0)  # (4, N+kn)

        #DPj2 a 4 lignes : les trois premieres pour les coef de la base splines,
        #et la dernieres pour le coef de la variable z.
       
        A_j_array[:, :, j] = K1 @ DPj2
        c_j_array[j, :]    = K2 @ DPj2

        # Optimisation avec CVXPY
    
    alpha = cp.Variable(N)
    z = cp.Variable(kn)
    
    # Fonction objectif
    residuals = ytab - env_array @ alpha
    
    objective = cp.Minimize(rhotau(residuals, tau))
    #u_plus  = cp.Variable(n)     # Résidus positifs
    #u_minus = cp.Variable(n)    # Résidus négatifs
    # Fonction objectif
    #objective = cp.Minimize(tau * cp.sum(u_plus) + (1 - tau) * cp.sum(u_minus))
    #objective=cp.Minimize(cp.norm(residuals))
    constraints=[]
        # Contraintes de base
    #constraints0 = [
    #    u_minus - u_plus == residuals,
    #    u_plus >=0,
    #    u_minus >= 0  ]

    # Contraintes de monotonie 
    for j in range(kn):
        if monot[j] != 0:
            A_mat = A_j_array[:, :, j]  
            c_vec = c_j_array[j, :]     
            # Contrainte de cône second-order
            if np.any(A_mat != 0):  # Éviter les contraintes vides
             monotone_vec = cp.hstack([alpha*monot[j] , z])
             #constraints.append(c_vec @ monotone_vec >= cp.sqrt((cp.sum_squares(A_mat @ monotone_vec))))
             constraints.append(c_vec @ monotone_vec >= cp.norm(A_mat @ monotone_vec)) 
             constraints.append(c_vec @ monotone_vec >=0)
             constraints.append(z[j] >= 0)
    
    # Contraintes de convexité
    for j in range(kn + 1):
        if cv[j] != 0:
            constraints.append(cv[j] * (con_array[j, :] @ alpha) >= 0)
    
    # Résolution du problème
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False,solver='GUROBI')
    #prob.solve(verbose=False,solver='MOSEK')
    if alpha.value is None:
        warnings.warn("L'optimisation n'a pas convergé")
        return None
    
    # Construction de la spline résultante
    l_end = np.min(knots)
    r_end = np.max(knots)
    
    # Extension de la séquence de nœuds
    s = np.concatenate([[l_end]*3, knots, [r_end]*3])
   
    polyn = BSpline(s, alpha.value, m)
        
    return polyn

# Exemple d'utilisation
#if __name__ == "__main__":

def test():
    # Données d'exemple
    np.random.seed(42)
    n_points = 100
   
#    tau=0.5
    xtab = np.linspace(0, 1, n_points)
    #ytab=xtab*(xtab-1)*(xtab-0.5)
    ytab = 2*xtab+0.1*np.sin(10* np.pi * xtab)  +0.05 * np.random.randn(n_points)
    #ytab=np.exp((-5+10*xtab))/(1+np.exp((-5+10*xtab))) + 0.1 * np.random.randn(n_points) # logistique
    # concave

    kn=15
    knots=kn+1
    if isinstance(knots, int):
        # Si seul le nombre de nœuds est donné
        kn = knots-1
        # Calcul des quantiles pour les nœuds
        t_k_n = np.quantile(xtab, np.linspace(0, 1, kn + 1))


        knots=t_k_n.copy()
        fig, axes = plt.subplots(1, 2, figsize=(15, 4))
    # Appel de la fonction
    monot=[0]*kn
    cv=[0]*(kn+1)   
    for k in range(kn):
        if k/kn<4:
            monot[k]=0
        else:
            cv[k]=0
    #monot=1      
    #cv=0
    for tau in [0.1,0.5,0.99]:
     spline_result = SplineCubicQuantBspkn3(xtab, ytab, knots, tau, monot, cv)
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
            polyn=BSpline(spline_result.t,coef,3)
            axes[1].plot(x_eval,polyn(x_eval))
            pp=PPoly.from_spline(polyn)
            axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.5)

        
#        plt.plot(xtab,ytab,'r*')
#        plt.plot(x_eval,y_eval)
    plt.show()
    #print(spline_result.t,spline_result.c)
    return spline_result
        
