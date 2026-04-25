library(pracma)
library(polynom)
#include('change_base.R')
library(splines)

Omega<-function(s,j,l,k=j)#t: knots in the base t-t[j]
    {
  if(s[j]==s[j+l-1]){w<-c(0,0)}
  else
 {alpha<-1/(s[j+l-1]-s[j])
  w<-c(-alpha,s[j]-s[k]) } # This is the affine element wik of deBoor
  #wik=(t-t_i)/(t_{i+k-1}-t{i})=(t-t_j)*(q+(t_j-t_i)*q
  #p=as.polynomial(w)
  print(c(j,l,k))
  print
  return(w)
}

Omega0<-function(s,j,o)#t: knots in the base t-t[j]
{
  if(s[j]==s[j+o-1]){w<-c(0,0)}
  else
  {alpha=1/(s[j+o-1]-s[j])
  beta=(-s[j]*alpha)
  w<-c(beta,alpha) } # This is the affine element wjl of deBoor
  #wik=(t-t_i)/(t_{i+k-1}-t{i})=(t-t_j)*(q+(t_j-t_i)*q
 
  print(w)
  return(w)
}

polymul <- function(p1, p2) {
  if ((sum(abs(p1))==0) || (sum(abs(p2))==0)){ return(0)}
  else{deg1 <- length(p1) - 1
  deg2 <- length(p2) - 1
  res <- numeric(deg1 + deg2 + 1)
  for (i in 0:deg1) {
    for (j in 0:deg2) {
      res[i + j + 1] <- res[i + j + 1] + p1[i + 1] * p2[j + 1]
    }
  }
  return(res)}
}


polyadd<- function(p1, p2) {
  l1 <- length(p1)
  l2 <- length(p2)
  if (l1>l2) {p=p2
  p2=p1
  p1=p
  l=l2
  l2=l1
  l1=l}
  p1=c(p1,rep(0,(l2-l1)))
  return(p1+p2)
}



change_polynomial_base_taylor <- function(coeffs_a, a, b) {
  n <- length(coeffs_a) - 1
  
  # Calculer les valeurs du polynôme et ses dérivées au point b
  # P(t) = sum c_k (t-a)^k
  # On utilise la formule de Taylor : c'_j = P^{(j)}(b)/j!
  
  coeffs_b <- numeric(n + 1)
  
  for (j in 0:n) {
    # Calculer P^{(j)}(b) = sum_{k=j}^{n} c_k * k!/(k-j)! * (b-a)^{k-j}
    deriv_val <- 0
    for (k in j:n) {
      if (abs(coeffs_a[k + 1]) > 1e-10) {
        deriv_val <- deriv_val + coeffs_a[k + 1] * 
          factorial(k) / factorial(k - j) * 
          (b - a)^(k - j)
      }
    }
    coeffs_b[j + 1] <- deriv_val / factorial(j)
  }
  
  return(coeffs_b)
}



Bspline_base<-function(sn,degree=3,deriv=0)
{#' this computes a Bspline basis coefficients, on the local power bases. 
#'sn is the extended knot vector
#' including the ends of th interval. This means if t0..t_{kn} it the set of knots 
#' then sn should be given as a vector with  "degree" times t_0 and t_{kn} at the begining and the ends.
#'  its length is number of intervals+1+2*degree
#'This   
 
  tn=sn[(degree+1):(length(sn)-degree)] #effective knots partition
  kn=length(tn)-1 # tn is the interior knots without the extended partition.
  
  n_intervals<-kn+2*degree #Nb extended intervals
  n_splines<-kn+degree
  B<-array(0,dim=c((degree+1),n_splines,n_intervals,(degree+1))) # B is the initial B-spline basis : piecewise constant
  for (i in (degree+1):(kn+degree)){B[1,i,i,1]<-1} # initialisation of the splines
  #with degree 0 wich has dimension kn
  # 1: degree 0, kn = Nb interval; kn= Nb elements in the basis of degree 0
  #the coefficients are in increasing order, contrary to the functions in RR like ppval, etc.
  if (degree>0){
    for (o in (2:(degree+1))){
      #k : dimension of local basis=deg+1
      for (j in (1:(n_splines))) {
      #go through the elements of the basis
        for (nu in ((degree+1):(n_intervals))) #go through the  pieces of the spline of order l
{
            Bjnu<-B[(o-1),j,nu,1:(o-1)]
            
            if ((j+1)>n_splines){Bjpnu=0
            wjp1nu=0}else
            {Bjp1nu=B[(o-1),(j+1),nu,1:(o-1)]
            wjp1o=polyadd(c(1,0),-Omega0(sn,(j+1),o))}
            wjo=Omega0(sn,j,o)
              
            term1=polymul(Bjnu,wjo)
            term2=polymul(Bjp1nu,wjp1o)
            
            sumterm=polyadd(term1,term2)
            
            B[o,j,nu,(1:o)]<-sumterm

           }
        }
    }}
  
  #changement de base 
  # because the coefficients of each polynomial piece
  #are computed on the canonical basis 1,x,x^2...x^degree 
  Bn=array(0,dim=dim(B))
  for (o in 1:(degree+1)){
  for (j in 1:n_splines)
  {for (nu in 1:n_intervals)
  {
    Bn[o,j,nu,]<-change_polynomial_base_taylor(B[o,j,nu,],0,sn[nu])
  }}}
  

  B=B[(degree+1),,,]
  Bn=Bn[(degree+1),,,]
  B=round(B,10)
  Bn=round(Bn,10)
  
  return(list(base =Bn ,base0=B, knots = sn, int_knots=tn, degree = degree, n_splines = (n_splines) ) )
#  return (B)
}


polyderiv<-function(p,der=1) # this is a working polyder function, 
  #because the order higher than 1 do not work  in the R function polyder

{if (der==0)
  {return(p)}
else{
 q=rev(p) #reverse the coefficient in order
 #to match the R convention of coefficients p=c(a3,a2,a1,a0) for a0+a1x+a2x^2+a3x^3...
 l=length(p)
 dl=l-der
  for (i in 1:der)
    {
  q=polyder(q)
  }

if (length(q)==1){
if ((l-der)>0 & (q==0)){
  q<-c(zeros(1,l-der))
  }}
return(rev(q)) # reverse back to match our own convention
 #(python style) with coefficients (a0,a1,a2,a3)
}
}

Bspline_deriv<-function(bspline,der=2){
#computes the derrivative fo a Bspline basis
  Bn=bspline$base
  B0=bspline$base0
  n_splines=bspline$n_splines
  knots=bspline$knots

  degree=bspline$degree
  degree_der=max(degree-der,0)
  NS=length(knots)-1 #Nb extended intervals
  Bn_der=array(dim=c(n_splines,NS,max((degree_der+1),1) ),0)
  
  B0_der=array(dim=dim(Bn_der),0)
  
  for (j in 1:n_splines){
    for (nu in 4:NS){
      p=polyderiv(Bn[j,nu,],der)
      q=polyderiv(B0[j,nu,],der)
      if (!is.null(p)){
        print(c(j,nu))
       Bn_der[j,nu,]=p
       B0_der[j,nu,]=q
       }
      }
  }

  return(list(base =Bn_der, base0=B0_der, knots = knots, int_knots=bspline$int_knots, degree = degree_der, n_splines = (n_splines) ) )
}

Spline_der_knots<-function(Bspline,der=1)
  #compute the values of a derivatives only at the knots (simpler, only uses the coefficients)
  {
  coeff=Bspline$base
  nsplines=Bspline$n_splines
  kn=Bspline$int_knots
  m=Bspline$degree
  if (der>m){
    Der2_knots=zeros(nsplines,kn)
    }
  else{
  Der2_knots=coeff[,,(der+1)]*factorial(der)
  }
  return(Der2_knots)
}


spline_eval<-function(Bspline,xvalues)
#Bspline has a new type R container, 
#désign by its coefficients, the degree 
#and the knots
  
{
  knots=Bspline$int_knots #interior knots
  degree=Bspline$degree
  coefficients=Bspline$coefficients
  
  BB=bs(xvalues,knots=knots,degree)
  #BB2=bs_direct(xvalues,)
  N=length(knots)+degree-1
  yvalues=BB[,1:N]%*%coefficients
  return(yvalues)
}


bs_direct<-function(Basis,xvalues)
{
  #Calcule les valeurs d'une base Bspline.
  #comme la fonction bs de R, mais en utilisant la
  # base calculée sous PP-forme : coeff des polynômes sur la base locale.
  
  n=length(xvalues)
  int_knots=Basis$int_knots#interior knots
  kn=length(knots)
  degree=Basis$degree
  nsplines=Basis$n_splines
  bb=Basis$base[,(degree+1):(nsplines),(degree+1):1] #internal bspline knots/ reverse coeff because of the stupd convention in R
  yvalues=zeros(nsplines,n)
  
  for (j in 1:nsplines)
    {
    p=mkpp(t(bb[j,,]),x=c(int_knots))
    yvalues[j,]<-ppval(p,xvalues)
  }
  
  
  
  return(yvalues)
}


view_spline<-function(Bspline,xvalues)
 {#permet de tracer la spline Bspline
  #Bspline est un type ad-hoc avec les 
  #coefficients, les noeuds, le degré.
  
   yvalues=spline_eval(Bspline,xvalues)
   matplot(xvalues, yvalues)
}