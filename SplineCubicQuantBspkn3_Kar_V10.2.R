# Chargement des librairies nécessaires
library(CVXR)
library(splines)
library(splines2)
library(fda)      # Pour les B-splines
library(pracma)

#Monotonicity and convexity constraints OK
#Version 2, 24/04/2026
#Code in R, similar to the python version, but the Bspline basis functions have
#been rewriten.


# Fonction de perte quantile. Non utilisee
rhotau <- function(u, tau) {
  return(sum(tau * pos(u) + (1 - tau) * pos(-u)))
}

#' Convertit une B-spline en représentation polynomiale 
#' par morceaux en appelant la fontion Bspline_base

bspline_to_deriv_coeffs_pp <- function(tn,degree = 3) {
  
  # Créer la base avec create.bspline.basis
  kn <- length(tn) - 1
  
  # Nombre correct de fonctions de base: kn + degree+1
  nbasis <- kn + degree
  
  # Créer la base B-spline avec fda
  norder <- degree + 1  # 4 pour cubique
  
  # Créer la base avec create.bspline.basis
  sn=c(tn[1]*ones(1,degree),tn,tn[kn+1]*ones(1,degree)) # this is the extended knots sequence
  BB <- Bspline_base(sn,degree)
  basis<-BB$base
  
  N <-BB$n_splines
  
  cat("Nombre de fonctions de base (fda):", N, "\n")
  
  # Matrice des coefficients pour la dérivée première normalisée
  
  deriv_coeffs <- array(0, dim = c(kn, N, 3))
  deriv2_val<-array(0, dim = c(kn+1, N))
  
  for (j in 1:N){
    for (nu in (degree+1):(kn+degree))
      {

      h=sn[nu+1]-sn[nu]
      a3<-3*basis[j,nu,4]*h^2
      a2<-2*basis[j,nu,3]*h
      a1<-basis[j,nu,2]
      c1<-2*basis[j,nu,3]
      
      # coeffs_poly est [a0, a1, a2, a3] a0+a1*x+a_2*x^2...
      deriv_coeffs[nu-degree,j,]=c(a3,a2,a1)
      deriv2_val[nu-degree,j]=c1
    }
    # for the last knot the second deriv is an affine function
    # c1+c2(t-t_{kn-1}) h is the last intervall space
    c2=6*basis[j,nu,4]
    deriv2_val[nu-degree+1,j]=c1+c2*h
}

  return(list(d1=deriv_coeffs, d2=deriv2_val))
}



#' Applique les contraintes de Karlin pour un polynôme quadratique
#' 
apply_karlin_constraints <- function(p0, p1, p2, z0) {
  # P2, p1, p0 sont les coefficients du polynôme quadratique: p2*u^2 + p1*u + p0
  # Dans la notation de l'article

  constraints <- list()
  constraints <- c(constraints, list(z0 >= 0))
  
  K1_vec <- vstack(p0 + p2 - z0,p1-z0)
  K2_vec <- (p0 +p2+ z0)
  
    constraints <- c(constraints, list(K2_vec >= p_norm(K1_vec, 2)))
  
  return(constraints)
}

#' Régression quantile avec splines cubiques - Version avec contraintes de Karlin
#' pour la monotonie
#' et contraintes de convexite aux noeuds.
#' 
SplineCubicQuantBspkn3_Karlin <- function(xtab, ytab, knots, tau, 
                                          monot = 0,
                                          convcons=0,
                                          solver = "CLARABEL", weight = NULL) 
  {
  
  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }
  
  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]
  
  n <- length(xtab)
  
  if (length(knots) == 1 && is.numeric(knots))
    {
    kn <- knots - 1
    knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }
  
  kn <- length(knots) - 1
  
  cat("knots:", knots, "\n")
  
  if (length(monot) == 1) {
    monot <- rep(monot, kn)
  }
  cat("Monotonicity constraints (Karlin):", monot, "\n")
  
  
  boundary_knots <- range(knots)
  degree=3
  N=length(knots)+3-1
  B <- bs(xtab, knots = knots, degree = degree)
  #        Boundary.knots = boundary_knots, intercept = TRUE)
  #kn=length(knots)-1
  #N <- kn+degree

  int_knots=knots[2:kn]
  D2B=dbs(x=knots,derivs=2,knots=int_knots, degree = 3, Boundary.knots = range(knots))
  #cat("Nombre de fonctions de base:", N, "\n")
  B=B[,1:N]
  # Calcul des coefficients normalisés des dérivées avec fda
  deriv_spline <- bspline_to_deriv_coeffs_pp(knots, degree = 3)
  deriv_coeffs <-deriv_spline$d1
  deriv_coeffs2<-deriv_spline$d2
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean
  
  alpha <- Variable(N)
  #Contraintes monotones  
  # Variables auxiliaires z
  residuals <- ytab_centered - (B %*% alpha)
  
   u_plus <- pos(residuals)
   u_minus <- pos(-residuals)
   weighted_loss <- sum(weight * (tau * u_plus + (1 - tau) * u_minus))
  
  objective <- Minimize(weighted_loss)
  
  constraints <- list()
  z_vars <- list()
  
  if (any(monot != 0)) {
    for (i in 1:(kn)) {
      if (monot[i] != 0) {
      z_vars[[i]] <- Variable(1, name = paste0("z", i))
       a_coef=sum(deriv_coeffs[i,,1]*alpha) *monot[i]
       b_coef=sum(deriv_coeffs[i,,2]*alpha) *monot[i]
       c_coef=sum(deriv_coeffs[i,,3]*alpha) *monot[i]
      
      CK<-apply_karlin_constraints(c_coef,b_coef,a_coef,z_vars[[i]])      
      constraints<-c(constraints,CK)
      }
    }
  }
 
  #"contraintes convexes
  if (length(convcons) == 1) {
    convcons <- rep(convcons, (kn+1))
  }
  # eliminate the null (unconstrained) case
  
  if (any(convcons !=0)){
    CV<-list(convcons*(deriv_coeffs2 %*% alpha)>=0) # Very simple, only use the second derivatives at the knots.
    constraints<-c(constraints,CV)
  }
  
  problem <- Problem(objective, constraints)
  
  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")
  
  for (s in unique(solvers_to_try)) {
    cat("Tentative avec solveur:", s, "\n")
    result <- tryCatch({
      solve(problem, solver = toupper(s), verbose = FALSE)
    }, error = function(e) {
      cat("Échec:", e$message, "\n")
      NULL
    })
    
    if (!is.null(result) && !is.null(result$getValue(alpha))) {
      cat("Solveur réussi:", s, "\n")
      break
    }
  }
  
  if (is.null(result) || is.null(result$getValue(alpha))) {
    warning("L'optimisation n'a pas convergé avec aucun solveur")
    return(NULL)
  }
  
  alpha_val <- result$getValue(alpha)+y_mean
  
#  cat("Statut:", result$status, "\n")
#  cat("Valeur objectif:", result$value, "\n")
#  cat("Coefficients alpha (range):", range(alpha_val), "\n")
  

  return(list(
    #spline = spline_result,
    coefficients = alpha_val,
    degree=3,
    #basis_matrix = B,
    knots = knots,
    int_knots = knots
  ))
}

#' Fonction de test simplifiée
test_karlin_simple <- function() {
  set.seed(42)
  n_points <- 50
  xtab <- linspace(0,1,n_points)

  # Données simples et clairement croissantes
  #ytab <- -3 * xtab +sin(3*2*xtab*3.14)+ 0.2 * rnorm(n_points)
  ytab <- 2* xtab + 0.5 * sin(6 * pi * xtab) + 0.05 * rnorm(n_points)
  #ytab<-xtab*(1-xtab)
  kn <- 12
#  ytab<- xtab^2
#  kn <- 6

    n=7
  monot=c(ones(1,n),zeros(1,(12-n)))
  knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  res<-SplineCubicQuantBspkn3_Karlin(xtab, ytab, knots, tau = 0.9,
                                   monot = 0, solver = "OSQP")
  cat("\n=== TEST CROISSANT PARTIEL ===\n")

#  monot=0
  res_croissant <- SplineCubicQuantBspkn3_Karlin(xtab, ytab, knots, tau = 0.5,
                                                monot = monot, convcons=0,solver = "OSQP")
  
  cat("\n=== TEST DÉCROISSANT ===\n")
  res_decroissant <- SplineCubicQuantBspkn3_Karlin(xtab, ytab, knots, tau = 0.5,
                                                   monot = -1, solver = "OSQP")
cat("\n=== TEST CONVEXE ===\n")
res_convexe <- SplineCubicQuantBspkn3_Karlin(xtab, ytab, knots,monot=0,convcons=1, tau = 0.5)
                                                 #                                                   monot = -1, solver = "OSQP")
  # Visualisation
  par(mfrow = c(2, 2))
  x_eval <- seq(0, 1, length.out = 200)
  y_sans=spline_eval(res,x_eval)
  y_croiss=spline_eval(res_croissant,x_eval)
  y_decroiss=spline_eval(res_decroissant,x_eval)
  y_convexe=spline_eval(res_convexe,x_eval)

  # Croissant
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "Contrainte croissante (sur les premiers noeuds)")
  lines(x_eval,y_croiss,lwd=2)
  
#  view_spline
    #lines(x_eval, y_croiss, col = "red", lwd = 2)
    
  abline(v = knots, col = "blue", lty = 2)
  
  # Décroissant
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "Contrainte décroissante")
  
    lines(x_eval, y_decroiss, col = "red", lwd = 2)
  
  abline(v = knots, col = "blue", lty = 2)


# Convexe
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     main = "Contrainte convexe")
if (!is.null(res_convexe)) {
  lines(x_eval, y_convexe, col = "red", lwd = 2)
}
abline(v = knots, col = "blue", lty = 2)

# sans contraintes
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     main = "sans contrainte ")
lines(x_eval,y_sans,lwd=2)
abline(v = knots, col = "blue", lty = 2)

}
# Exécuter le test simplifié
test_karlin_simple()