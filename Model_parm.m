

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tinlet=304;                     %Tmin=150k

Tmin=150;                    %excel data

tmax_c=6*60*60;                 %charging time
tmax_s=5*60;                    %standby time
tmax_d=4.5*60*60;                 %discharge time
tmax=tmax_c+tmax_s+tmax_d;
dt=20;                          %in secs
nt=tmax/dt+1;

H=32/100;                       %height of bed
dx=0.008;                       %in m
xmax_c=H;
nx=H/dx+1;

t=linspace(0,1,nt+1);           %non dimensional vector t
x=linspace(0,1,nx+1);           %non dimensional vector X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_vessel=27.3/1000;                  %excel data in m^3...

eta=0.38;                       %eta from excel
Din=15/100;                     %Model pdf table 1
Rin=Din/2;
db=11.25/1000;                  %Model pdf table 1
C=110.4;                        %Model pdf table 1
rho_f=3.39;                     %Excel
rho_s=2688;                     %from excel data
cp_s=0.7*1000;                  %from excel data
cp_f=1.105*1000;                %from excel data
Req=0.25*eta*db/(1-eta);        %eqn 12 in model pdf
mdot=V_vessel*(eta*rho_f*cp_f+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_f);       %eqn 3 in kgp paper on exp setup                 %data from excel
mdot=1;
munot=24.2*10^(-6);             %eqn 15 model pdf
T0=423;                         %eqn 15 model pdf
muf=munot*((T0+C)/(Tinlet+C))*((Tinlet/T0)^1.5);  %eqn 15 model pdf
G=mdot/(eta*pi*Din*Din/4);      %eqn 9 model pdf
Re=4*G*Req/muf;                 %eqn 8 model pdf
Ac=eta*pi*Din*Din/4;            %eqn 3 model pdf
Afs=3*pi*Din*Din*(1-eta)/(2*db);%eqn 4 model pdf
v=mdot/(rho_f*Ac);              %eqn 10 model pdf
k_s=2.5;
pr=cp_f*muf/k_fT(Tinlet);
hfs=0.191*mdot*cp_f*(pr^(-2/3))*(Re^(-0.278))/(eta*pi*(Din*Din/4));
hfsv=hfs*Afs/Ac;
stfs=eta*hfsv*H/(rho_f*cp_fT(Tinlet)*v);
K=(rho_f*cp_f*eta)/(rho_s*cp_s*(1-eta));
peaf=rho_f*cp_f*H*v/(eta*k_fT(Tinlet));%non dimensional parameters pg 39 eqn 9.44
peas=rho_s*cp_s*H*v/(eta*k_s);          %non dimensional parameters pg 39 eqn 10.45
Dout=Din+0.02;
Rout=Dout/2;
kins=0.006;
hwall=hfsv;

parm(1)=Tinlet;
parm(2)=Tmin;
parm(3)=H;
parm(4)=V_vessel;
parm(5)=eta;
parm(6)=Din;
parm(7)=Rin;
parm(8)=db;
parm(9)=C;
parm(10)=rho_f;
parm(11)=rho_s;
parm(12)=cp_s;
parm(13)=cp_f;
parm(14)=Req;
parm(15)=mdot;
parm(16)=munot;
parm(17)=T0;
parm(18)=muf;
parm(19)=G;
parm(20)=Re;
parm(21)=Ac;
parm(22)=Afs;
parm(23)=v;
parm(24)=k_s;
parm(25)=pr;
parm(26)=hfs;
parm(27)=hfsv;
parm(28)=stfs;
parm(29)=K;
parm(30)=peaf;
parm(31)=peas;
parm(32)=Rout;
parm(33)=Dout;
parm(34)=kins;
parm(35)=hwall;


%%%
theta_f=ones(nt+1,nx+1);
theta_s=ones(nt+1,nx+1);
%i is temporal and j is positional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(1,:)=0;                 %initial condition
theta_s(1,:)=0;                 %initial condition      table 1 T_s(t=0,x)=31 degree==>Theta_s(1,:)=1


theta_f(:,1)=1;                 %boundary  condition  
theta_f(:,2)=theta_f(:,1);      %boundary  condition  d(theta_f)/dx=0 at x=0


for i=2:tmax_c/dt+1
    theta_s(i,1)=theta_s(i-1,1)+stfsT(phi(theta_f(i-1,1),parm),parm)*kT(phi(theta_s(i-1,1),parm),parm)*(dt/tmax_c)*(0-theta_s(i-1,1));   %boundary conditions from table
end

theta_s(:,2)=theta_s(:,1);      %boundary condition d(theta_s)/dx=0 at x=0;

%%%     using finite difference method       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% charging eqquations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:tmax_c/dt+1
    for j =2:nx
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax_c)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j),parm),parm)*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j)) - (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax_c)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j),parm),parm)*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j),parm),parm)*stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j),parm),parm)*(theta_s(i,j)-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% standby equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(tmax_c/dt+3:tmax_c/dt+tmax_s/dt+2,1:2)=theta_f(i,1);%cut off bc
theta_s(tmax_c/dt+3:tmax_c/dt+tmax_s/dt+2,1:2)=theta_s(i,1);%cut off bc

for i=tmax_c/dt+2:tmax_c/dt+tmax_s/dt+1
    for j =2:nx
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j),parm),parm)*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j)));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j),parm),parm)*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j),parm),parm)*stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j),parm),parm)*(theta_s(i,j)-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discharge equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_f(tmax_c/dt+tmax_s/dt+2:nt,1)=0; %inlet is at hot temp
theta_f(tmax_c/dt+tmax_s/dt+2:nt,2)=theta_f(tmax_c/dt+tmax_s/dt+2:nt,1);%d(theta_f)/dx=0 =>theta_f(t,1)=theta_f(t,2)

for i=tmax_c/dt+tmax_s/dt+2:nt
    for j =2:nx
        parm(10)=1.3;
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax_d)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j),parm),parm)*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j))- (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax_d)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j),parm),parm)*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j),parm),parm)*stfsT(phi(theta_f(i,j),parm),parm)*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j),parm),parm)*(theta_s(i,j)-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T5=26,T6=32,T7=2
clear title xlabel ylabel
hold on
%d1=xlsread('d.xlsx');
%x1=d1(:,4)/37380;
%plot(x1,d1(:,5),'blue-');
%x2=linspace(0,1,7458);
%plot(x2,d1(:,3),'red-')
plot(t,-(Tmin-Tinlet)*theta_s(:,20)+Tmin)
plot(t,-(Tmin-Tinlet)*theta_s(:,21)+Tmin)
plot(t,-(Tmin-Tinlet)*theta_s(:,22)+Tmin)
plot(t,-(Tmin-Tinlet)*theta_s(:,18)+Tmin)
plot(t,-(Tmin-Tinlet)*theta_s(:,19)+Tmin)
plot(t,-(Tmin-Tinlet)*theta_s(:,23)+Tmin)

%legend('X=0.2','X=0.5','X=0.8','location','southeast')
xlabel('Non Dimensional Time t/{\tau}')
ylabel('Dimensional Temp_f')
title('Temp_f Vs Time')
%plot(t,(Tmax-Tinlet)*theta_f(:,70)+Tinlet,'green--')
%plot(t,(Tmax-Tinlet)*theta_f(:,41)+Tinlet,'blue--')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=phi(t,parm)
    T=(parm(1)-parm(2))*t+parm(2);
end

function x=cp_fT(T)
    x=(4.475e-7)*T^3+9.584*(T^2)*(10^(-4))-0.4314*T+1061;
end

function x=k_fT(T)
    x=-(3.679e-8)*T^2 + (9.789e-5)*T - 10^(-4);
end

function x=stfsT(T,parm)
eta=parm(5);db=parm(8);Rin=parm(7);mdot=parm(15);T0=parm(17);C=parm(9);munot=parm(16);Din=parm(6);H=parm(3);rho_f=parm(10);

    Ac=eta*pi*Din*Din/4;
    Afs=3*pi*Din*Din*(1-eta)/(2*db);
    v=mdot/(rho_f*Ac);
    Req=0.25*eta*db/(1-eta);
    G=mdot/(eta*pi*Rin^2);
    mu_fT=munot*(T0+C)*((T/T0)^1.5)/(T+C);
    ReT=4*G*Req/mu_fT;
    Pr=cp_fT(T)*mu_fT/k_fT(T);
    hfs=0.191*mdot*(Pr^(-2/3))*cp_fT(T)*(ReT^(-0.278))/(eta*pi*Rin^2);
    hfsv=hfs*Afs/Ac;
    x=hfsv*H*eta/(rho_f*cp_fT(T)*v);
end

function x=peafT(T,parm)
rho_f=parm(10);H=parm(3);eta=parm(5);mdot=parm(15);Din=parm(6);
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    x=rho_f*cp_fT(T)*H*v/(k_fT(T)*eta);
end

function x=peasT(T,parm)
rho_s=parm(11); cp_s=parm(12); H=parm(3); k_s=parm(24); eta=parm(5); mdot=parm(15); rho_f=parm(10); Din=parm(6);
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    x=rho_s*cp_s*H*v/(k_s*eta);
end

function x=kT(T,parm)
rho_f=parm(10); eta=parm(5); rho_s=parm(11); cp_s=parm(12);
    x=rho_f*cp_fT(T)*eta/(rho_s*cp_s*(1-eta));
end

function x=stwT(T,parm)
 Din=parm(6); Dout=parm(33); rho_f=parm(10); kins=parm(34); hwall=parm(35); eta=parm(5); Rout=parm(32); Rin=parm(7); rho_s=parm(11); H=parm(3); cp_s=parm(12); mdot=parm(15);
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    Uloss=1/(Din*log(Dout/Din)/(2*kins)+1/hwall);
    Ulossv=Uloss*2*Rout/((1-eta)*Rin*Rin);
    x=Ulossv*H*eta/(rho_s*cp_s*v);
end
