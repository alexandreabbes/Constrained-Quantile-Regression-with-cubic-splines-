%% setting 1%%

clear;


numtrial=1;
t0=clock   ;
funs=1;
%  figure(funs);
 % Please refer to the function 'Fron' with regards to the meaning of the variable 'funs'
n = 100; %number of samples

tauval=0.5;
%tauval=[0.1 0.3 0.5 0.7 0.9];
%tauval=[1:5]/5;

fit_total=zeros(length(tauval),4*2,numtrial);

%nkn=20; %bandwidth
betaa = 2;
betab=3;

rng default

for trial=1:numtrial
    
%xtab = sort(unifrnd(0,1,1,n));
xtab=[0.1,0.25,0.5,0.75,0.9]'; %attention, colonne. si les xi sont entre en ligne, ca cause un pb
ytab=1-xtab;

if funs==1
    sigma=0.1
else
    sigma=1.5
end
perturb=2;
qV=zeros(length(tauval));
%xgrid are the points where the quantile function is estimated

xgrid=xtab;


nkn=2;
q=quantile(xtab,(0:(nkn))/nkn);     %kn=4
kn=zeros(nkn,1);

%problem of chosing knots with quantiles, it can happen that no value in
%inter-knot

%The cv constraints count kn+1 values
%The monot values count kn values.

for t=1:(nkn+1)%redress the knots without interpolation
      
      kn(t)=min(xtab(xtab>=q(t)));

end




for configur=[1:1]
    %figure(configur);
   
    switch configur
    case 1
        %no constraint
        cv =0;
        monot=0;
        leg='no constaint';
         
        
        
    case 2
        cv=0;
        leg='Only monotoncity constaint';  
        monot=1;
    case 3
        cv=[ 1    -1];
        monot=1;
        leg='with monotonicity constaint, and convexity constraint';
        
    case 4
        cv=[1  -1];
        monot=0;
        leg='Only convexity constraint';
    end
    
fons=Fron(xgrid,funs);
fon=fons;%???
    
fit_t=[];%['tau'  'mse'  'made'];
fit_t2=[];
fit1=zeros(length(xgrid),length(tauval));
fit2=zeros(length(xgrid),length(tauval));
fit3=zeros(length(xgrid),length(tauval));

%knots=[0.1 0.3 0.5 0.8];
for ii=1:length(tauval)
   %parfor ii=1:length(tauval) % for paralised algorithm
   tau=tauval(ii);
   
   knmax=2
   tic
%   [spfit1 ,IC1,spMat, monotf, cvf,kmin1,kmin2 ]=SplineCubicQuantBSpkn3AutoAddGlobal(xtab,ytab,kn',knmax,tau,monot,cv);
   toc
   %fit1(:,ii)=fnval(spfit1,xgrid);
   %spfit1=SplineCubicQuantMatrix_c2(xtab,ytab,xgrid,kn,tau,monot,cv);
   %fit1(:,ii)=fnval(spfit1,xgrid);
   
   %fit1(:,ii) = SplineCubicQuantTrigo(xtab,ytab,xgrid,nkn,tau,monot,cv);
   %[fit1(:,ii),lambda] = SplineCubicQuantAutoremove(xtab,ytab,xgrid,nkn,tau,monot,cv);
   %fit2(:,ii)=SplineCubicQuant(xtab,ytab,xgrid,nkn,tau,monot,cv);
   
   SplineCubicQuantBspkn(xtab,ytab,xgrid,kn,tau,monot,cv);
   %[fit1(:,ii),alpha1]=SplineCubicQuantMatrix_c(xtab,ytab,xgrid,kn,tau,monot,cv);
   %fit1(:,ii)=SplineCubicQuantTrigo_knots(xtab,ytab,xgrid,nkn,N,tau,monot,cv);
   %spfit1=SplineCubicQuantMatrixAutoRemove2(xtab,ytab,xgrid,kn,tau,monot,cv);
   %fit1=fnval(spfit1,xgrid);
   %[fit2(:,ii),alpha2]=SplineCubicQuantMatrix_c(xtab,ytab,xgrid,kn,tau,monot,cv);
end
 


for ii=1:length(tauval)
    tau=tauval(ii);
    fit0=fit1(:,ii);
    %fit0_2=fit2(:,ii);


% %exp V is the random coeff of the y_i's

%t %true quantile function
 %true quantile function

% if trial==1
%     scatter(xtab,ytab,'.');
%     hold on;
%     plot(xgrid,fon,'-');
% 
%     
%     quantplot=plot(xgrid,fit0, 'Color','k', "Linewidth",2);
%     %quantplot=plot(xgrid,fit0,'Color','k', "Linewidth",2);
%     %quantplot2=plot(xgrid,fit0_2,'--','Color','g', "Linewidth",2);
%         
% 
% plot(xgrid,Fron(xgrid,funs),'--','color','k' );
% plot(xgrid,true,'--','color','b' );
% 
% title(leg);
% hold off;
% end
ngrid=length(xgrid);
mse=(sum((true-fit0').^2))/ngrid;
rmse=sqrt(mse);
%rnmse=sqrt((sum((true-fit0').^2)/sum((true').^2)));
%mean absolute deviation error
made=sum(abs(true-fit0'))/ngrid;
%normalised made
%nmade=sum(abs(true-fit0'))/sum(abs(true));

%mse2=(sum((true-fit0_2').^2))/length(true);
%rmse2=sqrt(mse2);
%nmse2=(sum((true-fit0_2').^2)/sum((true).^2));
%mean absolute deviation error
%made2=sum(abs(true-fit0_2'));
%normalised made
%nmade2=sum(abs(true-fit0_2'))/sum(abs(true));




% id=1:length(xgrid);
% id_plus=id(ytab-fit0'>0);
% id_minus=id(ytab-fit0'<0);
% %robustified criterion
% 
% mqe=(sum(tau*abs(ytab(id_plus)-fit0(id_plus)'))+sum((1-tau)*abs(ytab(id_minus)-fit0(id_minus)')))

%mqe=sum(rhotau(ytab-fit0',tauval(ii)))/sum(abs(ytab));
%mqe2=sum(rhotau(ytab-fit0_2',tauval(ii)))/sum(abs(ytab));
%tmqe=(sum(tau*abs(ytab(id_plus)-true(id_plus)))+sum((1-tau)*abs(ytab(id_minus)-ytab(id_minus))))/sum(abs(true-ytab));

fit_t=[fit_t ; [ rmse made]];
%fit_t2=[fit_t2 ; [nmse2 rmse2 nmade2 made2 ]];

end
%legend('data','tau=0.1,..,0.9 from bottom to top')

hold off;

display(fit_t);
%fit_total(:,1,trial)=tauval';
%fit_total2(:,1,trial)=tauval';
%array [tau,  4*4 columns] : 4 'sub arrays' of 4 columns 
%%Names={'tau','nc','monot','ccve_monot','ccve'}), each sub array [ nmse mse nmade made ]

fit_total(:,configur*2-1:configur*2,trial)=fit_t;
%fit_total2(:,configur*4-2:configur*4+1,trial)=fit_t2;
%display(fit_total);
   end
end



fit_stat_mean=[mean(fit_total,3)];
%fit_stat_mean2=[mean(fit_total2,3)];

fit_stat_std=[std(fit_total,0,3)];

display(fit_stat_mean);


t0=clock-t0;
