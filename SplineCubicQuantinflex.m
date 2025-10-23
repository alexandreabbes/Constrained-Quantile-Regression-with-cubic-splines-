function [fit,sc]=SplineCubicQuantinflex(xtab,ytab,knots,candidates,tau,monot,cv,weight)
   
   %we suppose that the cv=+1 indicates that the curve is convex then
   %concaves, and cv=-1 is the converse.
   %We assume the data are monotonic, and we impose the contrainsts with tge monot argument
   %the aim is to find the inflexion point. We impose a change of convexity
   %of the points in candidates. The best fit (measured with the asymetric
   %pseudo distance will give the result.)
   %we suppose that knots are only the interior knots
   %fit is a vector of spline object estimators
   %sc is the corresponding list of scores
   
   N=length(candidates);
   sc=zeros(N,1); 
   fit=[];
   
   for i=1:(N)
     C=candidates(i);
       knots_c=knots;


     %add the not C  
     if sum((C==knots_c))==0
         knots_c=sort([knots_c, C]);
     end    
     %remove the knot C
     %knots_c(knots_c==C)=[];
     
     
     %One adds 2 points to the set of knots: one before and one after the candidate
         C1=min(xtab(xtab>C));
         if sum((C1==knots_c))==0
         knots_c=sort([knots_c, C1]);
         end

         C0=max(xtab(xtab<C));
         
         if sum((C0==knots_c))==0
         knots_c=sort([knots_c, C0]);
         end  
     
     
     
     conv=(([0,knots_c,1]<C)*2-1)*cv;
     %set the security knots, arround the candidate C, at a 0 convexity constraint value.
     conv(knots_c==C)=0;
     conv(knots_c==C0)=0;
     conv(knots_c==C1)=0;
     display(conv);
     display(knots);
     fit0=[SplineCubicQuantBspkn3(xtab,ytab,knots_c,tau,monot,conv)];
     fit=[fit fit0];
     
     fit1=fnval(fit(i),xtab);
     %fnplt(fit0);
     %hold on;
     sc(i)=sum(rhotau(ytab-fit1,tau))/length(ytab);
     %ddfit=fnder(fit0,2);
     %sc(i,3)=min(fnzeros(ddfit))
     %sc(i,2)=candidates(i)

   end
 
%[m,argmin]=min(sc);
%fit_sol=fit(argmin);
%inflex=candidates(argmin)
%ddfit=fnder(fit_sol,2);
%display(fnzeros(ddfit));
end
