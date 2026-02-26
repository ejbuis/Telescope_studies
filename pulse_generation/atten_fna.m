function y=atten_fna(f,dkm,opt)
% ATTEN_FNA calculates attenuation in Sea water.
% inputs f  (Hz) and dkm (distance in km)
% outputs attenuation as a fraction of the pressure remaing at the particular frequency.
% Case 1 Learned's Technique based on Lehtinen et al. w0=2.5e10.
% Case 2 Polynomial fit to Fig 1 Lehtinen et al
% Case 3 Neiss and Bertin's Complex attenuation
% Case 4 Ainslie and McColm J. Acoust. Soc. AM. 103 (3) 1671 1998
% Case 5 No attenuation
% Case 6 Francois and Garrison
% Case 7 ACoRNe (Combined A&Mc and N&B);
% Case 8 Fischer & Simmons (Note there is a bug in this).
% Last Mod SD 26/1/07
if nargin==2;opt=1;end
switch opt
    case 1
        c =1500;omega=2.5e10;
        attendb= 1e4./log(10)/omega/c.*(2*pi*f).^2;
        y=10.^-((attendb*dkm)/20);
    case 2
        c =1500;omega=2.5e10;
        c4 =[  -0.00506349498933   0.07974277244217  -0.46228957401157   1.18153834673275  -1.13888823728799]*1e2;
        attendb=10.^polyval(c4,log10(abs(f+eps)));
        y=10.^-((attendb*dkm)/20);
        y(abs(f)<10)=1;
        attendb= 1e4./log(10)/omega/c.*(2*pi*f(abs(f)>1e5)).^2;
        yh=10.^-((attendb*dkm)/20);
        y(abs(f)>1e5)=yh;

    case 3
        c =1510;w0=0.79e12;
        w=2*pi*f;w1=2*pi*91.29e3;w2=2*pi*1.31e3;
        lambda1=157.8;lambda2=45.6e3;
        y=exp(-dkm*1e3*(w.^2/w0/c+1/lambda1*i*w./(w1+i*w)+1/lambda2*i*w./(w2+i*w)));
    case 4
        %f in kHz z in km
        T=15; S=37; pH=7.9; z=2;
        f1=0.78*sqrt(S/35)*exp(T/26);
        f2=42*exp(T/17);
        f=f*1e-3; % Convert to kHz
        attendb=0.106*f1.*f.^2./(f.^2+f1^2)*exp((pH-8)/0.56) +...
            0.52*(1+T/43)*(S/35)*f2.*f.^2./(f.^2+f2^2)*exp(-z/6)+...
            0.00049*f.^2*exp(-T/27+z/17);
        y=10.^-((attendb*dkm)/20);
    case 5
        y=ones(size(f));
    case 6
        %FG f is in kHz z im m
        % TODO  This does not seem to agree with A&Mc and N&B sttenuation
        % Suspect a Typo Somewhere Here which needs to be sortd
        % Do not rely on this at present!
        f=f*1e-3; % Convert to kHz
        T=15; S=37; pH=7.9; D=2000;   Theta=T+273;
        c=1412+3.21*T+1.19*S+0.0167*D;
        % Boric Acid Contribution
        A1=8.86/c*10^(0.78*pH-5);  P1=1; f1=2.8*sqrt(S/35)*10^(4-1245/Theta);
        % MgSO_4 Contribution
        A2=21.44*S/c*(1+0.025*T); P2=1 -1.37e-4*D + 6.2e-9*D^2;
        f2=8.17*10^(8-1990/Theta)/(1 +0.0018*(S-35));
        if T<20
            A3=polyval([-1.5e-8 9.11e-7 -2.59e-5 4.937e-4],T);
        else
            A3=polyval([-6.5e-10 1.45e-7 -1.146e-5 3.964e-4],T);
        end
        P3=polyval([4.9e-10 -3.83e-5 1],D);
        alpha_BoricA=A1*P1*f1*f.^2./(f.^2+f1^2);
        alpha_MgSO4=A2*P2*f2*f.^2./(f.^2+f2^2);
        alpha_water=A3*P3*f.^2;
        attendb=alpha_BoricA+alpha_MgSO4+alpha_water;
        y=10.^-((attendb*dkm)/20);
    case 7
        %Combined Complex with A&Mc (ACoRNE)
        T=15; S=37; pH=7.9; z=2;
        w=2*pi*f;
        s=sqrt(-1)*w;
        w_B=2*pi*0.78e3*sqrt(S/35)*exp(T/26);
        w_Mg=2*pi*42e3*exp(T/17);
        K_B=0.106e-3*exp((pH-8)/0.56)/2/pi;
        K_Mg= 0.52e-3*(1+T/43)*(S/35)*exp(-z/6)/2/pi;
        K_W=0.00049e-6*exp(-T/27+z/17)/4/pi^2;
        attendb=K_B.*w_B*s./(s+w_B) +...
            K_Mg*w_Mg.*s./(s+w_Mg)+...
            K_W*w.^2;
        y=10.^-((attendb*dkm)/20);
    case 8
        % Fischer + Simmonds
        %%TODO Need to complete this one rainy afternoon
        T=15;P=20;c=1520;S=35;
        cT=fliplr([55.9 -2.37 .0477 -0.000348]);
        cP=fliplr([1 -3.84e-4 7.57e-8]);
        a0=polyval(cT,T)*polyval(cP,P)/4/pi^2*1e-15*c;
        w1=2*pi*1320*(T+273)*exp(-1700/(T+273));
        cT1=fliplr([1.03 .0236 -.000522]);
        a1=S/35*polyval(cT1,T)*1e-8;
        l1=pi*sqrt(2)/c/a1;
        w2=2*pi*15.5e6*(T+273)*exp(-3052/(T+273));
        cP2=fliplr([1 -10.3e-4 3.7e-7]);
        cT2=fliplr([5.62 9.0752]);
        a2=S/35*polyval(cT2,T)*polyval(cP2,P)*1e-8;
        l2=pi*sqrt(2)/c/a2;

    otherwise
        error('Attenuation options must be between 1 and 7');
end
