function points=MCGEn(jpri,lscale,rscale,n,imethod)
%MCGen Table based Monte Carlo Generator
%Inputs 
% jpri - 2D-Histogam whose statistics we wish to mimic:  size mxn
% lscale - bin edges of the rows size: (m+1)*1
% rscale - bin edges of the columns size: 1*(n+1)
% n - the number of points to be produced
%Last Mod 27/1/07  SD
long=sum(jpri,2)+1e-10*sum(jpri(:));%add 1e-10 to ensure monotinic increase
if nargin==4,imethod='pchip';end
lpoints=sort(interp1([0; cumsum(long)/sum(long)],lscale,rand(n,1),imethod));
rpoints=zeros(n,1);
spos=1; 
nthrowa=histc(lpoints,lscale);
for i=1:length(nthrowa)-1
         nthrow=nthrowa(i);
     if nthrow>0
         radial=jpri(i,:)+1e-8;
         rpoint=interp1([0 cumsum(radial)/sum(radial)],rscale,rand(nthrow,1),imethod);
         rpoints(spos:(spos+nthrow-1))=rpoint;
         spos=spos+nthrow;
     end
end
points=[lpoints rpoints];
