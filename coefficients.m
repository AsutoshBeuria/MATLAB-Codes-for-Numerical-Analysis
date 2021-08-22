syms peaf rho_f cp_fT H v eta k_fT peaf rho_s
syms cp_s k_s k_fT K Req db mu_fT munot C T0
syms G ReT Req Pr hfs hfsv Afs Ac mdot Rin Din T

C=110.4 ;
T0=423;
H=0.32;
eta=0.38;
rho_s=2700;
rho_f=3.39;
cp_s=700;
Din=0.15;
Rin=Din/2;
k_s=2.5;
db=11.25/1000;
munot=24.2e-6;
mdot =24.88626192/1000;


Afs=3*pi*Din*Din*(1-eta)/(2*db);

Ac=eta *pi*Rin^2;

v=mdot/(rho_f*Ac);

peaf=rho_f*cp_fT*H*v/(k_fT*eta);
peaf=vpa(peaf,8)

peas = rho_s*cp_s*H*v/(k_s*eta);
peas=vpa(peas,8)

K=rho_f*cp_fT*eta/(rho_s*cp_s*(1-eta));
K=vpa(K,8)


Req=0.25*eta*db/(1-eta);

G=mdot/(eta*pi*Rin^2);

mu_fT=munot*(T0+C)*((T/T0)^1.5)/(T+C);

ReT=4*G*Req/mu_fT;

Pr=cp_fT*mu_fT/k_fT;

hfs=0.191*mdot*(Pr^(-2/3))*cp_fT*(ReT^(-0.278))/(eta*pi*Rin^2);

hfsv=hfs*Afs/Ac;


stfs=hfsv*H*eta/(rho_f*cp_fT*v);
stfs=vpa(stfs,8)



%k_fT=0.0138;
%T=150;
%cp_fT=1105;