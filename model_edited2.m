%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tinlet=150;                     %Tmin=150k

Tmax=273+31;                    %excel data

tmax_c=6*60*60;                 %charging time
tmax_s=1*60*60;                    %standby time
tmax_d=4.5*60*60;                 %discharge time
tmax=tmax_c+tmax_s+tmax_d;
dt=5;                          %in secs
nt=(tmax)/dt+1;

H=32/100;                       %height of bed
dx=0.002;                       %in m
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
munot=24.2*10^(-6);             %eqn 15 model pdf
T0=423;                         %eqn 15 model pdf
muf=munot*((T0+C)/(Tinlet+C))*((Tinlet/T0)^1.5);  %eqn 15 model pdf
G=mdot/(eta*pi*Din*Din/4);      %eqn 9 model pdf
Re=4*G*Req/muf;                 %eqn 8 model pdf
%Re=mdot*db/(pi*(din*din/4)*(1-eta)*muf);
Ac=eta*pi*Din*Din/4;            %eqn 3 model pdf
Afs=3*pi*Din*Din*(1-eta)/(2*db);%eqn 4 model pdf
v=mdot/(rho_f*Ac);              %eqn 10 model pdf
%k_f=-3.679*Tinlet*Tinlet*10^(-8)+9.789*Tinlet*10^(-5)-10^-4;
%k_s=-2.469*(Tinlet^3)*(10^(-8)) -(9.509*(Tinlet^2))-(0.124*Tinlet)+61.76;
k_s=2.5;
pr=cp_f*muf/k_fT(Tinlet);
hfs=0.191*mdot*cp_f*(pr^(-2/3))*(Re^(-0.278))/(eta*pi*(Din*Din/4));
hfsv=hfs*Afs/Ac;


stfs=eta*hfsv*H/(rho_f*cp_fT(Tinlet)*v);
K=(rho_f*cp_f*eta)/(rho_s*cp_s*(1-eta));
peaf=rho_f*cp_f*H*v/(eta*k_fT(Tinlet));%non dimensional parameters pg 39 eqn 9.44
peas=rho_s*cp_s*H*v/(eta*k_s);          %non dimensional parameters pg 39 eqn 10.45
Uloss=0;                                %considering no Uloss
stw=(Uloss*H*eta)/(rho_s*cp_s*v);

theta_f=ones(nt+1,nx+1);
theta_s=ones(nt+1,nx+1);

%i is temporal and j is positional



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(1,:)=1;                 %initial condition
theta_s(1,:)=1;                 %initial condition      table 1 T_s(t=0,x)=31 degree==>Theta_s(1,:)=1


theta_f(:,1)=0;                 %boundary  condition   
theta_f(:,2)=theta_f(:,1);      %boundary  condition  d(theta_f)/dx=0 at x=0


for i=2:tmax_c/dt+1
    theta_s(i,1)=theta_s(i-1,1)+stfs*kT(phi(theta_s(i-1,1)))*(dt/tmax_c)*(0-theta_s(i-1,1));   %boundary conditions from table
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
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax_d)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peafT(phi(theta_f(i,j)))*(dx/H)*(dx/H))-stfs*(theta_f(i,j)-theta_s(i,j))- (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax_d)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peasT(phi(theta_s(i,j)))*(dx/H)*(dx/H)) +kT(phi(theta_f(i,j)))*stfsT(phi(theta_f(i,j)))*(theta_f(i,j)-theta_s(i,j))-stwT(phi(theta_s(i,j)))*(theta_s(i,j)-1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear title xlabel ylabel
%hold on
plot(t,(theta_f(:,:)*(Tmax-Tinlet)+Tinlet))
%plot(t,(theta_f(:,round(0.5*nx))*(Tmax-Tinlet)+Tinlet))
%plot(t,(theta_f(:,142)*(Tmax-Tinlet)+Tinlet))
d1=xlsread('d.xlsx');
x=d1(:,4)/37380;
%plot(x,d1(:,5));
%legend('{\tau}','2/3{\tau}','{\tau}');
%legend('@0.8 X','(2/_3)X','X');
%title('T_f vs non-dimensional time');
%ylabel('{\theta}_s');
%ylabel('Dimensional T_f');
%xlabel('Non-dimensional Distance X')
%xlabel('Non-dimensional time (t / t_{total})');
%d1=xlsread('d.xlsx');
%x=linspace(0,1,7458);
%plot(linspace(0,1,7458),d1(:,2),'-');
%plot(linspace(0,1,7458),d1(:,3),'-');
%legend('x=16cm cal','x=16cm exp');
%legend('x=2cm (T5)cal','x=28cm (T7)cal','x=2cm (T5) exp','x=28cm (T7) exp');
%hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=phi(t)
    Tmax=273+31;Tinlet=150;
    T=(Tmax-Tinlet)*t+Tinlet;
end

function x=cp_fT(T)
    x=(4.475e-7)*T^3+9.584*(T^2)*(10^(-4))-0.4314*T+1061;
end

function x=k_fT(T)
    x=-(3.679e-8)*T^2 + (9.789e-5)*T - 10^(-4);
end

function x=mdotT(T,eta,V_vessel,rho_f,rho_s,cp_s)
   x=V_vessel*(eta*rho_f*cp_fT(T)+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_fT(T));
end

function x=stfsT(T)
    x=20.210347/(((1.9796104*T + 218.54899)/(0.0023640662*T)^(3/2))^(139/500)*((0.01290828*cp_fT(T)*(0.0023640662*T)^(3/2))/(k_fT(T)*(T + 110.4)))^(2/3));
end

function x=kT(T)
    x= 0.0000010993344*cp_fT(T);
end

function x=peasT(~)
    x=695972.61;
end
function x=stwT(~)
    x=0;
end
function x=peafT(T)
    x=(3.1208296*cp_fT(T))/k_fT(T);
end