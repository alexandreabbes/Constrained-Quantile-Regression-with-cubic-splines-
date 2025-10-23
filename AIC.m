function ic=AIC(vect,deg,type)
n=length(vect);
k=deg;
if (~exist( 'type','var'))
      type=1;
end

switch type
    case 1 %HS 98
    ic=log(sum(abs(vect)))+2*(k+3)/n;
    case 2 %HS95 or global-M
    ic=(sum(abs(vect)))+2*(k+3)/n;
end
