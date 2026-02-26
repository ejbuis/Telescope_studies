function tsmc=ShowerParm(r,z,Eo,type)
% Shower Parameterisation Function
% Currently Supports
% 'CORSIKA'     - average of 100 showers
% 'SVDParam'    - SVD Parameterisation
% 'Niess'       - Niess and Bertin's Parameterisation
% 'Saund'       - Saund Parameterisation
% 'Sloan'       - Terry's Parameterisation 3/1/07
% 'Cylinder'    - Cylindrical Parameterisation for simple studies sigma=3cm
%
% Inputs:       r [cm] radial bin centres
%               z [cm] longitudinal bin centres
%               Eo [GeV] Energy (between 10^5 and 10^12 GeV)
% SD/SB last mod 7/7/08
F=str2func(type);
tsmc=F(r,z,Eo);
%tsmc=feval(type,r,z,Eo);
%%
%----------------------------------------------------
function tsmc=CORSIKA(r,z,Eo)
%Loads the Mean showers from file
% r and z are dummies
s=load('Sm'); %Sm is created by GenSm
tsmc=s.Sm(:,:,round(2*log10(Eo)-9));
tsmc=tsmc';
%-----------------------------------------------
%%
function tsmc=SVDParm(r,z,Eo)
%Energies is a vector containing the energies for which the MC points should
%be generated GeV *10^energies
% Note r and z are dummy variables as a constant calling syntax is required
energy=log10(Eo);
fluct=1;% change this to zero if no fluctuations
load Herwig %This is created by GenHParam
n=sum(nv);wmc=[];breaks=[0 cumsum(nv)];
%vsregen=interp1(5:.5:12,means',energy);
if fluct==0
    vsregen=interp1(5:.5:12,means,energy);
else
    vsregen=vsregen+interp1(5:.5:12,stds,energy).*(randn(length(energy),n)*spline(5:.5:12,cmra,energy));
end
for i=1:length(nv)
    wmc1=ws(:,breaks(i)+1:breaks(i+1))*vsregen(:,breaks(i)+1:breaks(i+1))';
    wmc=[wmc wmc1(:)];
end
tsmc=max(0,wmc*v');
%%
%%----------------------------------------------------
function y=Niess(r,z,Eo)
% Niess and Bertin's radial paramaterisation function
% z distance along shower core (cm)
% r radial distance from shower core (cm)
% Eo primary energy in GeV (Used to Calculate zmax)
% Zmax - shower peak. (optional)
[rg,zg]=meshgrid(r,z);
Ec=0.05427; Xo=35.29;   b=0.56;
zpmax=0.65*log(Eo/Ec)+3.93;
zmax=Xo*zpmax;
a=b*zpmax+1;
%Radial Distribution
x=3.5./rg;n1=1.66-0.29.*zg./zmax;n2=2.7;
y=x.^n1.*(x>1)+x.^n2.*(x<=1);
y=y.*rg;y=diag(1./sqrt(sum(y.^2,2)))*y;
%longitidinal distribution
zp=zg./Xo;
long=Eo./Xo*(b*zp).^(a-1).*exp(-b*zp)/gamma(a);
y=y.*long;
%%
%-------------------------------------------------------------
function tsmc=SAUND(r,z,Eo)
% Saund Parameterisation from Sloan and Vanderbruke
r=r(:)';z=z(:);
Xo=36.1; %cm
Ec=0.0838 ;%GeV
r_m=9.04; %cm
s=1.25;
zmax=0.9*Xo*log(Eo/Ec);
lambda=130-5*log10(Eo*1e-4);
t=zmax/lambda;
k=t.^(t-1)./exp(t)./lambda./gamma(t);
a=r./r_m;
rho=1./(r_m.^2).*(a.^(s-2)).*((1+a).^(s-4.5)).*gamma(4.5-s)./(2*pi.*gamma(s).*gamma(4.5-2*s));
rho=rho.*r*2*pi;
rho=rho/sum(rho);
zdist=k.*(z./zmax).^t.*exp(t-z./lambda);
tsmc=Eo*zdist*rho;
%----------------------------------------------------------
%%
function radial = Sloan(r,z,Eo)
%Terry's Latest Parameterisation 3/1/07
r=r(:)';z=z(:);
El=log10(Eo);
% Predicts the energy deposit per 20 g cm^-2 slice in z per g cm^-2 in radius
% d2E/drdz from the parameterisation. Fit over 10^4 - 10^12 GeV but the fit
% is poor for E < 10^6.5 GeV.
% EL=log10E depth=Z radius=R
% Compute RM from fits in NKG parameterisation.
Co=[ 0.1287E-01,-0.2573E+00, 0.9636E+00;     -0.4697E-04, 0.8072E-03,  0.5404E-03; 0.7344E-07,-0.1375E-05, 0.4488E-05];
Ci=[-0.8905E-03,      0.7727E-02, 0.1969E+01;  0.1173E-04,-0.1782E-03,-0.5093E-05;  -0.1058E-07, 0.1524E-06,-0.1069E-08];
C2=[ 0.5519E-01,-0.3948E-02, 0.1490E-03;-0.2109E+03,-0.6968E-02,  0.1551E+00;-0.4150E+02, 0.1139E+03,-0.4103E+01;
    0.8012E+01, 0.1144E+02,-0.5434E+00; 0.7999E-05,-0.4843E-02, 0.2552E-03; 0.4563E-04,-0.3504E-05, 0.1315E-06];
abc=Co*(El.^[2 1 0]');                      
abc2=Ci*(El.^[2 1 0]');
rm=[ones(size(z)) z z.^2]*abc;          
s=[ones(size(z)) z z.^2]*abc2;
r_div_rm=r(ones(1,size(rm,1)),:)./rm(:,ones(1,size(r,2)));
%Analytical Integral
a1=4.5-2*s;a2=4.5-s;gt=real(rm.*gamma(a1).*gamma(s).*gamma(a2));
gt(gt<=0)=mean(gt>0)*1e6;
gt(end)=1e6;
s=s(:,ones(1,size(r,2)));
ffr=real(((r_div_rm).^(s-1.)).*((1.+(r_div_rm)).^(s-4.5)));
ffr(ffr<0)=0;
ffr=diag(1./gt)*ffr;
anorm=HILLIT(El,z,C2);
anorm=anorm*(10.0^El);  % Since function HILLIT is normalised by E.
d2Edrdz=diag(anorm)*ffr;
radial=d2Edrdz;
radial(end-1:end,:)=radial(end-3:end-2,:);% Subtle Bug at 1e10Gev???
%% c---------------------------------------------------------
function HILLITRet=HILLIT(alE,z,C2)
% Compute the shower profile for parameters P in the Hillas function.
% HILLITRet = Energy deposited at z z gm/cm^2 per 20 g cm^-2 slice.
% p1 = Energy deposited at the maximum. p(1)= p1/E  
% p2 = Shower start z.
% p3 = Maximum z of the shower. 
% p4,p5,p6 = Polynomial describing the width of the distribution.
% alE=Log of energy, z=z of the shower in water.
p=C2*(alE.^[0 1 2]');
a1=max(0,(z-p(2)));
a=((a1)/(p(3)-p(2)));
b=((p(3)-p(2)))./(p(4)+p(5)*z+p(6)*z.^2);
c=p(1)*(a.^b);
d=exp((p(3)-z)./(p(4)+p(5)*z+p(6)*z.^2));
HILLITRet=c.*d;
%% ----------------------------------------------
function tsmc=Cylinder(r,z,Eo)
tsmc=ones(size(z(:)))*exp(-(4*r/sqrt(2)).^2);
tsmc=tsmc/sum(tsmc(:));
%%-----------------------------------------------
