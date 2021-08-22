
clear variables;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Tinlet;Tinlet=150;                     %Tmin=150k

global Tmax;Tmax=273+31;                    %excel data

tmax_c=6*60*60;                 %charging time
tmax_s=5*60;                    %standby time
tmax_d=4.5*60*60;                 %discharge time
tmax=tmax_c+tmax_s+tmax_d;
dt=10;                          %in secs
nt=tmax/dt+1;

global H;H=32/100;                       %height of bed
dx=0.004;                       %in m
xmax_c=H;
nx=H/dx+1;

t=linspace(0,1,nt+1);           %non dimensional vector t
x=linspace(0,1,nx+1);           %non dimensional vector X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global V_vessel;V_vessel=27.3/1000;                  %excel data in m^3...
global eta;eta=0.38;                       %eta from excel
global Din;Din=15/100;                     %Model pdf table 1
global Rin;Rin=Din/2;
global db;db=11.25/1000;                  %Model pdf table 1
global C;C=110.4;                        %Model pdf table 1
global rho_f;rho_f=3.39;                     %Excel
global rho_s;rho_s=2688;                     %from excel data
global cp_s;cp_s=0.7*1000;                  %from excel data
global cp_f;cp_f=1.105*1000;                %from excel data
global mdot;mdot=V_vessel*(eta*rho_f*cp_f+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_f);       %eqn 3 in kgp paper on exp setup                 %data from excel
global munot;munot=24.2*10^(-6);             %eqn 15 model pdf
global T0;T0=423;                         %eqn 15 model pdf
global Ac;Ac=eta*pi*Din*Din/4;            %eqn 3 model pdf
global Afs;Afs=3*pi*Din*Din*(1-eta)/(2*db);%eqn 4 model pdf
global k_s;k_s=2.5;
global Dout;Dout=Din+0.02;
global Rout;Rout=Dout/2;
global kins;kins=0.25;
global hwall;hwall=0.5;
%%%
theta_f=ones(nt+1,nx+1);
theta_s=ones(nt+1,nx+1);
%i is temporal and j is positional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(1,:)=1;                 %initial condition
theta_s(1,:)=1;                 %initial condition      table 1 T_s(t=0,x)=31 degree==>Theta_s(1,:)=1

theta_f(:,1)=0;                 %boundary  condition  
theta_f(:,2)=theta_f(:,1);      %boundary  condition  d(theta_f)/dx=0 at x=0

for i=2:tmax_c/dt+1
    theta_s(i,1)=theta_s(i-1,1)+stfsT(phi(theta_f(i-1,1)))*kT(phi(theta_s(i-1,1)))*(dt/tmax_c)*(0-theta_s(i-1,1));   %boundary conditions from table
end

theta_s(:,2)=theta_s(:,1);      %boundary condition d(theta_s)/dx=0 at x=0;

%%%     using finite difference method       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% charging eqquations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:tmax_c/dt+1
    for j =2:nx
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax_c)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j)))*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j)) - (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax_c)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j)))*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j)))*stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j)))*(theta_s(i,j)-1));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% standby equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(tmax_c/dt+3:tmax_c/dt+tmax_s/dt+2,1:2)=theta_f(i,1);%cut off bc
theta_s(tmax_c/dt+3:tmax_c/dt+tmax_s/dt+2,1:2)=theta_s(i,1);%cut off bc


for i=tmax_c/dt+2:tmax_c/dt+tmax_s/dt+1
    for j =2:nx
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j)))*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j)));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j)))*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j)))*stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j)))*(theta_s(i,j)-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discharge equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_f(tmax_c/dt+tmax_s/dt+2:nt,1)=1; %inlet is at hot temp
theta_f(tmax_c/dt+tmax_s/dt+2:nt,2)=theta_f(tmax_c/dt+tmax_s/dt+2:nt,1);%d(theta_f)/dx=0 =>theta_f(t,1)=theta_f(t,2)

for i=tmax_c/dt+tmax_s/dt+2:nt
    for j =2:nx
        rho_f=1.3;
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax_d)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j)))*(dx/H)*(dx/H))-stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j))- (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax_d)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j)))*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j)))*stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j)))*(theta_s(i,j)-1));
        a = stfsT(phi(theta_f(i,j)))*kT(phi(theta_f(i,j)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear title xlabel ylabel
hold on
plot(t,(theta_f(:,:)*(Tmax-Tinlet)+Tinlet))
%d1=xlsread('d.xlsx');
%x=d1(:,4)/37380;
%plot(x,d1(:,5));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=phi(t)
global Tinlet;global Tmax;
    T=(Tmax-Tinlet)*t+Tinlet;
end

function x=cp_fT(T)
    x=(4.475e-7)*T^3+9.584*(T^2)*(10^(-4))-0.4314*T+1061;
end

function x=k_fT(T)
    x=-(3.679e-8)*T^2 + (9.789e-5)*T - 10^(-4);
end

%function x=mdotT(T,eta,V_vessel,rho_f,rho_s,cp_s)
%    x=V_vessel*(eta*rho_f*cp_fT(T)+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_fT(T));
%end

function x=stfsT(T)
global eta;global db;global Rin;global mdot;global T0;global C;global munot;global Din;global H;global rho_f;
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

function x=kT(T)
global rho_f;global eta;global rho_s;global cp_s;
    x=rho_f*cp_fT(T)*eta/(rho_s*cp_s*(1-eta));
end

function x=peasT(T)
global rho_s;global cp_s;global H;global k_s;global eta;global mdot;global rho_f;global Din;
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    x=rho_s*cp_s*H*v/(k_s*eta);
end
function x=peafT(T)
global rho_f;global H;global eta;global mdot;global Din;
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    x=rho_f*cp_fT(T)*H*v/(k_fT(T)*eta);
end

function x=stwT(T)
global Din;global Dout;global rho_f;global kins;global hwall;global eta;global Rout;global Rin;global rho_s;global H;global cp_s;global mdot;
    Ac=eta*pi*Din*Din/4;
    v=mdot/(rho_f*Ac);
    Uloss=1/(Din*log(Dout/Din)/(2*kins)+1/hwall);
    Ulossv=Uloss*2*Rout/((1-eta)*Rin*Rin);
    x=Ulossv*H*eta/(rho_s*cp_s*v);
end
