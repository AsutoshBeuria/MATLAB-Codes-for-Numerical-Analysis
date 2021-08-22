format LONG
%%%Using the method given in Transient model pdf (Euler) Explicit method like DP%%not working for all points especially for discharge part%%
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tinlet=150;                     %Tmin=150k
Tmax=273+31;                    %excel data
tmax_c=6*60*60;                 %charging time
tmax_s=5*60;                    %standby time
tmax_d=4.5*60*60;                 %discharge time
tmax = tmax_c + tmax_s + tmax_d;
dt = 5;                          %in secs
nt=(tmax_c+tmax_s+tmax_d)/dt+1;
H=32/100;                       %height of bed
dx=0.001;                       %in m
xmax_c=H;
nx=H/dx+1;
t=linspace(0,1,nt+1);           %non dimensional vector t
x=linspace(0,1,nx+1); %non dimensional vector X
theta_f=ones(nt+1,nx+1);
theta_s=ones(nt+1,nx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_vessel=27.3/1000;                  %excel data in m^3...
eta=0.38;                       %eta from excel
din=15/100;                     %Model pdf table 1
db=11.25/1000;                  %Model pdf table 1
C=110.4;                        %Model pdf table 1
rho_f=3.39;                     %Excel
rho_s=2688;                     %from excel data
cp_s=0.7*1000;                  %from excel data
cp_f=1.105*1000;                %from excel data
Req=0.25*eta*db/(1-eta);        %eqn 12 in model pdf
Ac=eta*pi*din*din/4;            %eqn 3 model pdf
Afs=3*pi*din*din*(1-eta)/(2*db);%eqn 4 model pdf
             %eqn 10 model pdf
k_s=2.5;    %non dimensional parameters pg 39 eqn 9.44
      %non dimensional parameters pg 39 eqn 10.45


%i is temporal and j is positional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_f(1,:)=1;                 %initial condition
theta_s(1,:)=1;                 %initial condition      table 1 T_s(t=0,x)=31 degree==>Theta_s(1,:)=1

theta_f(:,1)=0;                 %boundary  condition  
theta_f(:,2)=theta_f(:,1);      %boundary  condition  d(theta_f)/dx=0 at x=0

%theta_s(:,2)=theta_s(:,1);      %boundary condition d(theta_s)/dx=0 at x=0;
%%%     using finite difference method       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% charging eqquations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:tmax_c/dt+1
    for j =1:nx
        mdot=V_vessel*(eta*rho_f*cp_fT(Tem(theta_f(i,j)))+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_fT(Tem(theta_f(i,j))));       %eqn 3 in kgp paper on exp setup                 %data from excel
        munot=24.2*10^(-6);             %eqn 15 model pdf
        T0=423;v=mdot/(rho_f*Ac);                   %eqn 15 model pdf
        muf=munot*((T0+C)/(Tem(theta_f(i,j))+C))*((Tem(theta_f(i,j))/T0)^1.5);  %eqn 15 model pdf
        G=mdot/(eta*pi*din*din/4);      %eqn 9 model pdf
        Re=4*G*Req/muf;                 %eqn 8 model pdf%Re=mdot*db/(pi*(din*din/4)*(1-eta)*muf);
       
%k_f=-3.679*Tinlet*Tinlet*10^(-8)+9.789*Tinlet*10^(-5)-10^-4;
%k_s=-2.469*(Tinlet^3)*(10^(-8)) -(9.509*(Tinlet^2))-(0.124*Tinlet)+61.76;
        pr=cp_fT(Tem(theta_f(i,j)))*muf/k_fT(Tem(theta_f(i,j)));
        hfs=0.191*mdot*cp_fT(Tem(theta_f(i,j)))*(pr^(-2/3))*(Re^(-0.278))/(eta*pi*(din*din/4));
        hfsv=hfs*Afs/Ac;
        %tmax = H*eta/v;
        stfs=hfsv*H/(rho_f*cp_fT(Tem(theta_f(i,j)))*v);
        K=(rho_f*cp_fT(Tem(theta_f(i,j)))*eta)/(rho_s*cp_s*(1-eta));
        peaf=rho_f*cp_fT(Tem(theta_f(i,j)))*H*v/(eta*k_fT(Tem(theta_f(i,j))));
        peas=rho_s*cp_s*H*v/(eta*k_s);
        Uloss=4000;                                %considering no Uloss
        stw=(Uloss*H*eta)/(rho_s*cp_s*v);
        if i ~= 1 && j == 1 && j == nx
           theta_s(i,1)=theta_s(i-1,1)+stfs*K*(dt/tmax)*(0-theta_s(i-1,1));   %boundary conditions from table
           theta_s(i,nx+1)=theta_s(i-1,nx+1)+stfs*K*(dt/tmax)*(1-theta_s(i-1,nx+1));
           %theta_f(i,nx+1) = theta_f(i,nx);
        end
        if j ~= 1
        theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peaf*(dx/H)*(dx/H))-stfs*(theta_f(i,j)-theta_s(i,j)) - (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
        theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peas*(dx/H)*(dx/H)) +K*stfs*(theta_f(i,j)-theta_s(i,j))-stw*(theta_s(i,j)-1));
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% standby equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=tmax_c/dt+2:tmax_c/dt+tmax_s/dt+1
    for j =1:nx
        mdot=V_vessel*(eta*rho_f*cp_fT(Tem(theta_f(i,j)))+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_fT(Tem(theta_f(i,j))));       %eqn 3 in kgp paper on exp setup                 %data from excel
        munot=24.2*10^(-6);             %eqn 15 model pdf
        T0=423;v=mdot/(rho_f*Ac);                   %eqn 15 model pdf
        muf=munot*((T0+C)/(Tem(theta_f(i,j))+C))*((Tem(theta_f(i,j))/T0)^1.5);  %eqn 15 model pdf
        G=mdot/(eta*pi*din*din/4);      %eqn 9 model pdf
        Re=4*G*Req/muf;                 %eqn 8 model pdf%Re=mdot*db/(pi*(din*din/4)*(1-eta)*muf);
       
%k_f=-3.679*Tinlet*Tinlet*10^(-8)+9.789*Tinlet*10^(-5)-10^-4;
%k_s=-2.469*(Tinlet^3)*(10^(-8)) -(9.509*(Tinlet^2))-(0.124*Tinlet)+61.76;
        pr=cp_fT(Tem(theta_f(i,j)))*muf/k_fT(Tem(theta_f(i,j)));
        hfs=0.191*mdot*cp_fT(Tem(theta_f(i,j)))*(pr^(-2/3))*(Re^(-0.278))/(eta*pi*(din*din/4));
        hfsv=hfs*Afs/Ac;
        %tmax = H*eta/v;
        stfs=hfsv*H/(rho_f*cp_fT(Tem(theta_f(i,j)))*v);
        K=(rho_f*cp_fT(Tem(theta_f(i,j)))*eta)/(rho_s*cp_s*(1-eta));
        peaf=rho_f*cp_fT(Tem(theta_f(i,j)))*H*v/(eta*k_fT(Tem(theta_f(i,j))));
        peas=rho_s*cp_s*H*v/(eta*k_s);
        Uloss=4000;                                %considering no Uloss
        stw=(Uloss*H*eta)/(rho_s*cp_s*v);
        if j ~= 1
           theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peaf*(dx/H)*(dx/H))-stfs*(theta_f(i,j)-theta_s(i,j)));
           theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peas*(dx/H)*(dx/H)) +K*stfs*(theta_f(i,j)-theta_s(i,j))-stw*(theta_s(i,j)-1));
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_f(tmax_c/dt+tmax_s/dt+2:nt,1)=1; %inlet is at hot temp
theta_f(tmax_c/dt+tmax_s/dt+2:nt,2)= theta_f(tmax_c/dt+tmax_s/dt+2:nt,1);

for i=tmax_c/dt+tmax_s/dt+2:nt
    for j = 1:nx
        mdot=V_vessel*(eta*rho_f*cp_fT(Tem(theta_f(i,j)))+(1-eta)*rho_s*cp_s)/(342*rho_f*cp_fT(Tem(theta_f(i,j))));       %eqn 3 in kgp paper on exp setup                 %data from excel
        munot=24.2*10^(-6);             %eqn 15 model pdf
        T0=423;v=mdot/(rho_f*Ac);                   %eqn 15 model pdf
        muf=munot*((T0+C)/(Tem(theta_f(i,j))+C))*((Tem(theta_f(i,j))/T0)^1.5);  %eqn 15 model pdf
        G=mdot/(eta*pi*din*din/4);      %eqn 9 model pdf
        Re=4*G*Req/muf;                 %eqn 8 model pdf%Re=mdot*db/(pi*(din*din/4)*(1-eta)*muf);
       
%k_f=-3.679*Tinlet*Tinlet*10^(-8)+9.789*Tinlet*10^(-5)-10^-4;
%k_s=-2.469*(Tinlet^3)*(10^(-8)) -(9.509*(Tinlet^2))-(0.124*Tinlet)+61.76;
        pr=cp_fT(Tem(theta_f(i,j)))*muf/k_fT(Tem(theta_f(i,j)));
        hfs=0.191*mdot*cp_fT(Tem(theta_f(i,j)))*(pr^(-2/3))*(Re^(-0.278))/(eta*pi*(din*din/4));
        hfsv=hfs*Afs/Ac;
        %tmax = H*eta/v;
        stfs=hfsv*H/(rho_f*cp_fT(Tem(theta_f(i,j)))*v);
        K=(rho_f*cp_fT(Tem(theta_f(i,j)))*eta)/(rho_s*cp_s*(1-eta));
        peaf=rho_f*cp_fT(Tem(theta_f(i,j)))*H*v/(eta*k_fT(Tem(theta_f(i,j))));
        peas=rho_s*cp_s*H*v/(eta*k_s);
        Uloss=4000;                                %considering no Uloss
        stw=(Uloss*H*eta)/(rho_s*cp_s*v);
        if j > 1
           theta_f(i+1,j+1)=theta_f(i,j+1)+(dt/tmax)*((theta_f(i,j+1)-2*theta_f(i,j)+theta_f(i,j-1))/(peaf*(dx/H)*(dx/H))-stfs*(theta_f(i,j)-theta_s(i,j))- (theta_f(i,j+1)-theta_f(i,j))/(dx/H));
           theta_s(i+1,j+1)=theta_s(i,j+1)+(dt/tmax)*((theta_s(i,j+1)-2*theta_s(i,j)+theta_s(i,j-1))/(peas*(dx/H)*(dx/H)) +K*stfs*(theta_f(i,j)-theta_s(i,j))-stw*(theta_s(i,j)-1));
        end
    end
            theta_s(i,1)=theta_s(i-1,1)+stfs*K*(dt/tmax)*(theta_f(i,1)-theta_s(i-1,1));   %boundary conditions from table
            theta_s(i,nx+1)=theta_s(i-1,nx+1)+stfs*K*(dt/tmax)*(theta_f(i,nx+1)-theta_s(i-1,nx+1));

    
   
     %theta_f(i,nx+1)=theta_f(i,nx);%d(theta_f)/dx=0 =>theta_f(t,1)=theta_f(t,2)
     %theta_f(i,1) = 1;
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear title xlabel ylabel
plot(t,(theta_f(:,:)*(Tmax-Tinlet)+Tinlet))
an=annotation('doublearrow',[0.13,0.903],[0.9,0.9]);
c=an.Color;
an.Color='r';
%legend('{\tau}','2/3{\tau}','{\tau}');
%legend('@0.8 X','(2/_3)X','X');
%legend('dx=0.0001 dt=0.1');
%title('Title');
%ylabel('{\theta}_s');
%ylabel('Dimensional temp_f');
%xlabel('Non-dimensional Distance X')
%xlabel('Non-dimensional time (t/{\tau})');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=cp_fT(T)
    x=4.475*(10^(-7))*T.^3+9.584*(T.^2)*(10^(-4))-0.4314*T+1061;
end
%function x=cp_sT(T)
%    x=5.467*(10^(-7))*(T^3)-1.967*(10^(-3))*T^2+2.453*T+1.913*100;
%end
%function x=k_sT(T)
%    x=-2.469*(10^(-8))*T^3-9.509*(10^(-5))*T^2  -0.124*T+61.76;
%end
function x=k_fT(T)
    x=-3.679*(10^(-8))*T.^2+9.789*(10^(-5))*T-10^(-4);
end
function x=mu_fT(T)
    x= munot*((T0+C)./(T+C))*((T/T0).^1.5);
end
function x=stT(T)
    x=0.191*H*eta*(1-eta)*3*pi*din*din*(k_fT(T).^(2/3)).*mu_fT(T).^(-2/3+0.278)*((4*mdot*db/(pi*din*din(1-eta)))^-0.278)./(2*db*(cp_fT(T).^(5/3)));
end
function x=Tem(y)
    Tmax = 273 + 31; Tinlet = 150;
    x = abs(y*(Tmax-Tinlet) + Tinlet);
end