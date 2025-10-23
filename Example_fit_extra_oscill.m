%% setting 1%%

clear;
a=cvx_solver;
if sum(a)~=sum('SeDuMi')
    cvx_solver SeDuMi
end

    
t0=clock   ;
funs=1;
% Please refer to the function 'Fron' with regards to the meaning of the variable 'funs'

 
n = 3; %number of samples
perturb=0;

tab=[1:30]/30;
%tauval=[0.1 0.3 0.5 0.7 0.9];
%tauval=[1:5]/5;

tauval=[ 0.5 ];
numtrial=1;

fit_total=zeros(length(tauval),4*2,numtrial);
Dfit_total=zeros(length(tauval),4*2,numtrial);
DDfit_total=zeros(length(tauval),4*2,numtrial);

%nkn=20; %bandwidth

betaa = 2;
betab=3;

qV=zeros(length(tauval));

for trial=1:numtrial
%ytab=data(:,1)';
%xtab=(0:length(ytab))/length(ytab);

%xeval = min(xtab)+(max(xtab)-min(xtab))*(0:0.001:1);

%xgrid = xeval(find((min(xtab)<=xeval).*(max(xtab)>=xeval)==1)); 
n=3
xtab=[0:n]/n;
%xgrid=xtab;
%xgrid=[1:100];
%ngrid=length(xgrid);



kn=(0:9)/10;

for configur=[1]
if trial==1
    figure(configur);
end
    switch configur
    case 1
        %no constraint
        %subplot(4,1,1);
        axis([0 1 -0.2 1.2]);
        cv =0;
        monot=0;
        leg='No constaint';
        %kn=[0,0.1,0.2,0.8,0.9];
        %kn=4;     
        
    case 2
        %subplot(4,1,2);
        axis([0 1 -0.2 2]);
        cv=0;
        %monot=-1;%[-1 -1 -1 -1 -1];
        leg='Only monotoncity constaint';  
        %if funs==4
            %knots=[0 0.17 0.29 0.41 0.56 0.71 0.88 1];
            %monot=[-1 0 1 0 0 0 -1];
            %kn=[ 0.2 0.5 0.7 ];
            %kn=4;
            monot=1;
         %   kn=[0,0.1,0.2,0.8,0.9];
        %end
        % for m=1:4
        
        %     monot=[monot [-1 -1 1 1]];
        % end
    case 3
%        kn=4; 
        %subpl;ot(4,1,3);
        axis([0 1 -0.2 1.2]);
        cv=[1 0  -1 ];
        %kn=[0,0.1,0.2,0.8,0.9];
        %cv=[-1 -1 0 0 0 0 0 0 0 0 1 1 1 -1 -1];
        %cv =[1 1 1  -1  -1 -1 -1 -1 ];

        %monot=[1 1 1 1 1 1 1 1 -1 -1 -1  -1 -1 -1 -1];
        %cv=1;
        monot=1;
        leg='Monotonicity constaint, and convexity constraint';
        
    case 4
        %subplot(4,1,4);
        cv=[1  0  -1];% 1 -1 -1 -1 -1]
        %cv =[ 1 1 1  -1  -1 -1 -1 -1 ];
        %kn=[0,0.1,0.2,0.8,0.9];
        %kn=5;
        %kn=4; 
        %subpl;ot(4,1,3);
        axis([0 1 -0.2 1.2]);
        %cv=[1 1  0 -1 -1];
        %cv=1;
        monot=0;
        leg='Only convexity constraint';

    end

    
 

if funs==1
    sigma=0.1
else
    sigma=1.5
end

switch perturb
    case 0
sigma=0.1;
        V=normrnd(0,sigma,1,n); %normal(0,1)
        perturbname="Normal";
        qV=(norminv(tauval,0,1));
        qV=quantile(V,tauval);
        ytab = Fron(xtab,funs);

    case 1
        V = (betarnd(betaa,betab,1,n)); %beta()
        perturbname="Beta";
        qV=betainv(tauval,betaa,betab);
        ytab = Fron(xtab,funs).*V;

    case 2 
        sigma=0.1;
        V=normrnd(0,sigma,1,n); %normal(0,1)
        perturbname="Normal";
        qV=(norminv(tauval,0,1));
        qV=quantile(V,tauval);
        ytab = Fron(xtab,funs);
    case 3

%gamma=0.2; mu=1;
% %cauchy centered in mu, with scale parameter gamma.

        V=mu+gamma*tan((unifrnd(0,1,1,n)-0.5)*pi);
        perturbname="Cauchy";
        ytab = Fron(xtab,funs).*V;

    case 4
        funs=5;
        nf=10;
        t_alpha=tinv(0.95,nf);
        for i=1:n
         Xi=normrnd(xtab(i)*1.5,1,[nf,100]);
        
%        X=sum(normrnd(xtab(i),1,[nf,100]).^2,1);
%        ytab(i)=1-sum(X<t_alpha-xtab(i)*1.5*sqrt(nf))/100;
        for k=1:100
            X(k)=sum(Xi(:,k))*sqrt(nf-1)/std(Xi(:,k));
        end
        beta=sum(X<t_alpha)/100;
        ytab(i)=1-beta;
        end
        
        
        
end



spfit1=[];
spfit2=[];
fit_t=[];%['tau'  'mse'  'made'];
Dfit_t=[];
DDfit_t=[];
fit_t2=[];

fit1=zeros(length(xgrid),length(tauval));
Dfit1=zeros(length(xgrid),length(tauval));
DDfit1=zeros(length(xgrid),length(tauval));
fit2=zeros(length(xgrid),length(tauval));
fit3=zeros(length(xgrid),length(tauval));

%knots=[0.1 0.3 0.5 0.8];
for ii=1:length(tauval)
   %parfor ii=1:length(tauval) % for paralised algorithm
   tau=tauval(ii);
   
   [Dtrue,DDtrue]=DFron(xgrid,funs);
   if funs~=5
      true=qV(ii)+Fron(xgrid,funs);;
   end
   
   %fit1(:,ii) = SplineCubicQuantTrigo(xtab,ytab,xgrid,nkn,tau,monot,cv);
   nkn=length(kn);
   if length(kn)==1
       nkn=kn;
   end
   %spfit1 = SplineCubicQuantAutoremove(xtab,ytab,5,tau,monot,cv);
   
   %fit2(:,ii)=SplineCubicQuant(xtab,ytab,xgrid,kn,tau,monot,cv);
   spfit1=[spfit1 SplineCubicQuantBspkn3(xtab,ytab,kn,tau,monot,cv)];
   
   fit1(:,ii)=fnval(spfit1(ii),xgrid);
   
   Dspfit1=fnder(spfit1(ii));
   Dfit1(:,ii)=fnval(Dspfit1,xgrid);
   
   DDspfit1=fnder(Dspfit1);
   DDfit1(:,ii)=fnval(DDspfit1,xgrid);
   %spfit is a spline
   %[fit1x(:,ii),alpha1,fithist]=SplineCubicQuantMatrixAuto(xtab,ytab,xgrid,tau,monot,cv);
   %[fit1(:,ii),alpha1]=SplineCubicQuantMatrixAutoRemove(xtab,ytab,xgrid,kn,tau,monot,cv);
   %spfit2=[spfit2 SplineCubicQuantMatrix_c2(xtab,ytab,kn,tau,monot,cv)];
   %fit2(:,ii)=fnval(spfit2(ii),xgrid);
   %spfit2=SplineCubicQuantMatrixAutoRemove2(xtab,ytab,xgrid,kn,tau,monot,cv);   
   %fit2(:,ii)=fnval(spfit2,xgrid);
   %[fit2(:,ii),alpha2]=SplineCubicQuantMatrix_kn(xtab,ytab,xgrid,kn,tau,monot,cv);
end

for ii=1:length(tauval)
    tau=tauval(ii);
    fit0=fit1(:,ii);
    Dfit0=Dfit1(:,ii);
    DDfit0=DDfit1(:,ii);
    %fit0_2=fit2(:,ii);


% %exp V is the random coeff of the y_i's

 %true quantile function
 %true quantile function

if trial==1
    scatter(xtab,ytab,'.');
    hold on;
    %plot(xgrid,true,'-'); 
    %quantplot=plot(xgrid,fit0, 'Color','k', "Linewidth",2);
    fnplt(spfit1(ii),'k', 2,[min(xtab), max(xtab)]);
    %quantplot2=plot(xgrid,fit0_2,'--','Color','g', "Linewidth",2);
%    fnplt(spfit2(ii),'--','g',2,[min(xtab),max(xtab)]);    

%plot(xgrid,Fron(xgrid,funs),'--','color','k' );
xplot=[0:100]/100;

plot(xplot,Fron(xplot,funs),'--','color','b' );

title(leg);
%hold off;
end







mse=(sum((true-fit0').^2))/(ngrid^2);
rmse=sqrt(mse);
Dmse=(sum((Dtrue-Dfit0').^2))/ngrid^2;
Drmse=sqrt(Dmse);
DDmse=(sum((DDtrue-DDfit0').^2))/ngrid^2;
DDrmse=sqrt(DDmse);


%rnmse=sqrt((sum((true-fit0').^2)/sum((true').^2)));
%mean absolute deviation error
made=sum(abs(true-fit0'))/ngrid;
%normalised made
%nmade=sum(abs(true-fit0'))/sum(abs(true));
%mse2=(sum((true-fit0_2').^2))/ngrid;
%rmse2=sqrt(mse2);
%mse2=(sum((true-fit0_2').^2)/sum((true).^2));
%mean absolute deviation error
%made2=sum(abs(true-fit0_2'))/ngrid;
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
Dfit_t=[Dfit_t ; [ Drmse 0]];
DDfit_t=[DDfit_t ; [ DDrmse 0]];
%fit_t2=[fit_t2 ; [rmse2  made2 ]];


end
%legend('data','tau=0.1,..,0.9 from bottom to top')

hold off;

%display(fit_t);
%fit_total(:,1,trial)=tauval';
%fit_total2(:,1,trial)=tauval';
%array [tau,  4*4 columns] : 4 'sub arrays' of 4 columns 
%%Names={'tau','nc','monot','ccve_monot','ccve'}), each sub array [ nmse mse nmade made ]

fit_total(:,configur*2-1:configur*2,trial)=fit_t;
Dfit_total(:,configur*2-1:configur*2,trial)=Dfit_t;
DDfit_total(:,configur*2-1:configur*2,trial)=DDfit_t;
%fit_total2(:,configur*2-1:configur*2,trial)=fit_t2;
display(fit_total(:,:,trial));
%display(fit_total2);

end

end




fit_stat_mean=[mean(fit_total,3)]
%fit_stat_mean2=[mean(fit_total2,3)];
Dfit_stat_mean=[mean(Dfit_total,3)]
DDfit_stat_mean=[mean(DDfit_total,3)]
fit_stat_std=[std(fit_total,0,3)]
Dfit_stat_std=[std(Dfit_total,0,3)]
DDfit_stat_std=[std(DDfit_total,0,3)]


t0=clock-t0;
