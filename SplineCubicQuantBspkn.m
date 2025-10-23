function [fitt,alpha]=SplineCubicQuantMatrixBspkn(xtab,ytab,x,kn,tau,monot,cv,weight)

%cv=-1 for convexity, cv=1 for concavity
%tau is the quantile parameter
%x is the grid where the fitted funtion is evaluated
%kn is either a number of internal knots, including 0, or a list of
%internal knots including 0
%monot is either +1/-1/0 for increasing/decreasing/no constraint, or a list
%of such values corresponding to the inter-knots intervals.
% CV is either +1/-1/0 for increasing/decreasing/no constraint, or a list
%of such values corresponding to the inter-knots intervals.


if (~exist( 'weight','var'))
        weight=ones(1,length(xtab));
end
    
l_end = min(xtab);
r_end = max(xtab);

%sort the x_i's
[xtab,so]=sort(xtab);
ytab=ytab(so);
n = length(xtab);
xgrid = x;


%% interior knots

%t_kn_i = quantile(xtab(idx),(1:kn)./(kn+1));

%t_kn_i = quantile(xtab,(1:kn)./(kn+1));

%t = [0 t_kn_i 1];

%t_kn = [l_end t_kn_i r_end];

%t=t_kn;
if (floor(kn)==kn)
    t_kn = (0:kn)./(kn);
    t_kn_i=t_kn(2:kn);
    
%t_kn = [l_end t_kn_i r_end];
else
    t_kn=[kn; 1];
    kn=length(kn);
    t_kn_i=t_kn(2:kn);
end

t=t_kn;

l_end=0;
r_end=1;

%extended knot partition:
s=[l_end l_end l_end t r_end r_end r_end];

if length(monot)==1
    monot=monot*ones(1,kn);
end


    if length(cv)==1
    cv=cv*ones(1,kn+1);
    end


c_j_array =[];

%% quantile constraints %%


m=3;
N=kn+m;



%Rescaling Matrix

%Rescaling+translation to [0,1]

%Restriction matrix : projection on Bsplines
%Pri=diag([ones(1,k) zeros(1,N-k)]);
Pr=zeros(kn,N);
for k=1:kn
    Pr(k,max(1,k-3):k)=1; 
end
%coeff Bsplines

Bsp=zeros(4,4,N);
Bsp2=zeros(4,4,N);

for j=1:N
    pp=bspline(s(j:j+4));% ind corresponds to the number of pieces odf the spline.
    bsp_size=size(pp.coefs);
    %if j>=4&j<=N-3
        ind_min=1;
        ind_max=4;
        %ind=1:4;
    %end
    
    if j<4
        %ind=4-j+1:4;
        ind_min=4-j+1;
    end
    if j>N-3
%    if j>N-4
        %ind=1:N-j+1;
        ind_max=ind_min+bsp_size(1)-1;
    end
    ind=ind_min:ind_max;
    display(ind);
    display(j);
%   Bsp(:,ind,j)=flip(pp.coefs'); % les coeffs sont ceux de la base (x-t(k-1))^l
    Bsp2(ind,:,j)=pp.coefs;
end
%Bsp(:,:,j) contient les coefs de la Bsplines entre les noeuds s(j) et
%s(j+4), morceau par morceau en colonne, sur la base des puissances locales
%avec les puissance croissantes

%differentiation Matrices
%D2=diag([1 2],1);

%D3=diag([1 2 3],1);
D32=diag([3 2 1],1);

%D3_proj=D3(1:3,:); %projection on R_2[X]
D32_proj=D32(:,2:4); %projection on R_2[X]

%D2_proj=D2([1 2],:); %projection on R_1[X]
D22=D32^2;
D22_proj=D22(:,3:4); %projection on R_1[X]


%Coefficient matrix derivative projected restricted scaled

%Pi=zeros(4,N,kn+1);
%DPi=zeros(3,N,kn+1);
%DDPi=zeros(2,N,kn+1);

Pi2=zeros(N,4,kn);
DPi2=zeros(N,3,kn);
DDPi2=zeros(N,2,kn);


for k=1:kn
    %over t(k),t(k+1)
%    P0=[];
    P02=[];
    for j=1:4
    j2=j-1;
%        P0=[P0 Bsp(:,4-j+1,k+j-1)];
        P02=[P02; Bsp2(4-j2,:,k+j2)];
    end
%    Rk=(t(k+1)-t(k)).^[0:2]'; %rescaling*
    Rk2=(t(k+1)-t(k)).^[2 1 0];
%    Pi(:,k:k+3,k)=P0; %coefficients on local power basis
    Pi2(k:k+3,:,k)=P02;
    DP02=(P02*D32);
    
%    DPi(:,k:k+3,k)=Rk.*D3_proj*(P0); %derivative, takes into account scaling to [0,1]
    DPi2(k:k+3,:,k)=DP02(:,2:4).*Rk2;
%    DDP0=D3^2*P0;  %differentiation
%    DDPi(:,k:k+3,k)=DDP0(1:2,:);
    DDP02=P02*D32^2;
    DDPi2(k:k+3,:,k)=DDP02(:,3:4); %DDPI2 are the coeff of the second derivative of the Bsplines.
end


%matrix 
env_array=zeros(n,N);
env_array2=zeros(n,N);
Sxtab=[]
v=1:length(s);%=1:N+3
for i = 1:n
  k=min([v(xtab(i)<s) (N+1)]); %find the index k of the interval containing xtab(i) : [s(k-1),s(k)]
  %k=min(k); pizero=[];
  pizero2=[];
  pizero=[];
  for l=1:4
 %   pizero=[pizero polyval(flip(Pi(:,k+l-1-4,k-4)),(xtab(i)-s(k-1)))]; %bspline (k-4) to k-1, value computed at xtab(i)-xtab(k-1)
    l2=l-1;
    pizero2=[pizero2 polyval(Pi2(k-4+l2,:,k-4),(xtab(i)-s(k-1)))];
    
  end
 
  Sxtab=[Sxtab k-1];
%  env_array(i,k-4:k-1) = pizero;
  env_array2(i,k-4:k-1) = pizero2;
end;

 
con_array=zeros(kn+1,N); %this is the matrix of the values of the 2cd derivatives of the basis (coef line) at each knot (coef column)



for k=1:kn %here k is the number of the interval associated to t.
%  pitwo=[];
  pitwo2=[];
    for L=0:3
        display(k+L);
 %       pitwo=[pitwo polyval(flip(DDPi(:,k+L,k)),0)]; %bspline k+l, FOR KNOT TK slice N�L To be reviewed.
        pitwo2=[pitwo2 DDPi2(k+L,2,k)]; %value at 0 
    end;

  %pow=(xtab(i)-s(k-1)).^[0:3];DPi
  %pizero=polyval(flip(Pi(:,:,k),2),(xtab(i)-s(k-1)));
  con_array(k,k:(k+3))= pitwo2;
end
%for the last knot bspline kn to kn+3
 % pitwo=[];
  pitwo2=[];
  k=kn;
 for L=0:3
 % pitwo = [pitwo polyval(flip(DDPi(:,kn+1+L,kn+1)),t(kn+2)-t(kn+1))]; %bspline (kn+p) slice kn+1, value at t(kn+2)
  pitwo2 = [pitwo2 polyval(DDPi2(kn+L,:,kn),(t(kn+1)-t(kn)))]; %bspline number kn+1+L, slice kn+1, value at kn+2
 end;
 
con_array(k+1,k:k+3)= pitwo2;


%Matrix Karlin&Studden
K1=[[1 0 -1 -1] ; [0 1 0 -1]];
K12=[[-1 0 1 -1] ; [0 1 0 -1]];
K2=[1 0 1 1];



A_j_array = zeros(2,N+kn,kn); % 
c_j_array=[];
for j=1:kn
    if and(j>1,j<kn)
        Ej=[zeros(1,j-1) 1 zeros(1,kn-j)];
    end
        if j==1
            Ej=[1 zeros(1,kn-1)];
        end
     if j==kn+1
            Ej=[zeros(1,kn-1) 1];
     end;
  %  DPj=[[DPi(:,:,j) zeros(3,kn+1)]; [zeros(1,N) Ej]];
    DPj2=[[flip(DPi2(:,:,j)') zeros(3,kn)];[zeros(1,N) Ej]];
    A_j_array(:,:,j)=K1*DPj2;
    c_j_array=[c_j_array;K2*DPj2];
end
%Standard SOCP formulation 
    % with or without the concavity constraints
    
    
    cvx_begin
    variable alph(N);
     variable z(kn);
     minimise weight*rhotau(ytab'-env_array2*alph,tau);

     subject to 
     %if monot~=0
     for j = 1:kn%+1)
        if monot(j)~=0
            norm(A_j_array(:,:,j)*[alph*monot(j); z]) <= (c_j_array(j,:)*[alph*monot(j);z]);%*monot(j);
            z(j)>=0;
        end
       
     end
     for j=1:kn+1
     
       cv(j)*(con_array(j,:)*alph) >= 0;%convexity on interval tj tj+1
     
     end
cvx_end

alpha=alph;

[r c]=size(xgrid);
fitt=[]

for i=1:length(xgrid)
   %v=1:kn+1
   k=min([v(xgrid(i)<t)]); %k-1 is the number of the interval where xgrid(i) lives.
   %The corresponding Bsplines to be used are in Pi(:,:,k), N=kn+1+3
   
   x=xgrid(i)-t(k-1);% Here is the mistake!!!! change it with t(k) or so...
   
   %k=min(k,kn)% xgrid(i)belongs to [s(k-1),s(k)]
   pol=alpha'*Pi2(:,:,k-1);
   fitt=[fitt polyval(pol,x)]; %
   
end
end
