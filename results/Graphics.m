load('variables')
Vgs=-2:0.5:2;
Vds=0:0.05:1;

 %Indicar la posición en la matriz de soluciones que se desea graficar.
 %Con estas coordenadas se establece la polarización del transistor.
 %=============== INDEXING FOR DIFFERENT GRAPHIC DISPLAY ==============
i=1; %1-21 (V_D=0V --> 1V, intervals 0.05V)
j=9; %1-9  (V_G= -2V --> 2V, intervals 0.5V)

%% Representación del gráfico Vd - I & Vg - I

figure(1)
plot(Vds,Ids(:,1))
hold on
for o=2:length(Vgs)
    plot(Vds,Ids(:,o))
end
hold off
legend('V_{G}=-2 V','V_{G}=-1.5 V','V_{G}=-1 V','V_{G}=-0.5 V','V_{G}=0 V','V_{G}=0.5 V','V_{G}=1 V','V_{G}=1.5 V','V_{G}=2 V')
xlabel('V_{D} (V)')
ylabel('I (A)')

figure(3)
plot(Vgs,Ids(1,:),'*--')
hold on
plot(Vgs,Ids(6,:),'*--')
plot(Vgs,Ids(11,:),'*--')
plot(Vgs,Ids(16,:),'*--')
plot(Vgs,Ids(21,:),'*--')
hold off
legend('V_{D}=0 V','V_{D}=0.25 V','V_{D}=0.5 V','V_{D}=0.75 V','V_{D}=1 V')
xlabel('V_{G} (V)')
ylabel('I (A)')


%% Definición de los datos del sistema
%Para ecuación de Poisson. 
tt=30e-9; %m
tg=0.6e-9; 
tb=300e-9;
L=500e-9;
LCs=50e-9;
Lexts=100e-9;
Lextd=100e-9;
LCd=50e-9;
Vs=0; %V
Vg=(j-5)*0.5;
Vd=(i-1)*0.05;
Vb=0;
e0=8.8541878176e-12; %F/m
et=5.9*e0; %TO
eb=3.9*e0; %BO
ey=3.3*e0; %GR y
ex=9*e0; %GR x
ny=70; %minimo 32
nx2=7; %número de nodos en el eje x en la lámina de grafeno.
% rho=linspace(0,0,ny);
rhoo=rho(:,i,j);

%Para ecuación diferencial física de semiconductores
h=6.626e-34/2/pi; %J/s
kB=1.38e-23; %J/K
vF=1e6; %m/s
q=1.6e-19; %C
T=300; %K

m0=0.2; %m^2/(V*s)
vsat=1e6; %m/s
W=100e-9; %m
Ids1=0; %menor (A)
Ids2=1; %mayor (A)

%% Discretización de nuestro dominio. Número de nodos en cada dirección. 
nx=ceil(0.01*nx2*(tt+tg+tb)/tg); %garantizamos nx2>7. 
%nx=30;


figure(2)

% plot(0:0.025:1,Ids(:,1))
% hold on
% plot(0:0.025:1,Ids(:,2))
% plot(0:0.025:1,Ids(:,3))
% hold off
% legend('Vbg=0V','Vbg=1V','Vbg=2V')

%Resolvemos las ecuaciones. Sacamos las variables para hacer los gráficos.


[phi,x,y,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rhoo,nx2);
[V,n,p,m,rhoo1,Ids]=semic(phiG,kB,T,h,vF,Vs,Vd,dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rhoo);

phigraf=zeros(ny,9);
Vgraf=zeros(ny,9);
for i=1:9
    [~,~,~,phigraf(:,i),dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,1,(i-5)*0.5,Vb,et,eb,ex,ey,rho(:,21,i),nx2);
    [Vgraf(:,i),~,~,~,~,~]=semic(phigraf(:,i),kB,T,h,vF,Vs,1,dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rho(:,21,i));
end


%Creamos los gráficos. 
figure(2)
tiledlayout(3,1)
nexttile
contourf(y,x,phi)
colorbar
title('Potencial electrostático')
xlabel('Eje y (m)')
ylabel('Eje x (m)')
zlabel('Potencial electrostático (V)')
nexttile
plot(V)
hold on
plot(phiG)
hold off
title('Potencial electroquímico y potencial electrostático en el grafeno')
xlabel('Eje y (m)')
ylabel('Potencial (V)')
legend('V','\phi')
nexttile
plot(rhoo)
title('Densidad de carga')
xlabel('Eje y (m)')
ylabel('Densidad de carga (C/m^3)')

figure(4)
tiledlayout(2,1)
nexttile
for i=1:9
    plot(y,Vgraf(:,i))
    hold on
end
yl = ylim;
xl = xlim;
xBox = [xl(1), xl(1), xl(1)+LCs, xl(1)+LCs];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.05);
xBox = [xl(1)+LCs, xl(1)+LCs, xl(1)+Lexts+LCs, xl(1)+Lexts+LCs];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'blue', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs, xl(1)+Lexts+LCs, xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'yellow', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L+Lextd, xl(1)+Lexts+LCs+L+Lextd];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs+L+Lextd, xl(1)+Lexts+LCs+L+Lextd, xl(2), xl(2)];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.05);
hold off
title('Potencial electroquímico en el grafeno (V_{D}=1 V)')
xlabel('Eje y (m)')
ylabel('V (V)')
legend('V_{G}=-2 V','V_{G}=-1.5 V','V_{G}=-1 V','V_{G}=-0.5 V','V_{G}=0 V','V_{G}=0.5 V','V_{G}=1 V','V_{G}=1.5 V','V_{G}=2 V')


nexttile
for i=1:9
    plot(y,phigraf(:,i))
    hold on
end
yl = ylim;
xl = xlim;
xBox = [xl(1), xl(1), xl(1)+LCs, xl(1)+LCs];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.05);
xBox = [xl(1)+LCs, xl(1)+LCs, xl(1)+Lexts+LCs, xl(1)+Lexts+LCs];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'blue', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs, xl(1)+Lexts+LCs, xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'yellow', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L, xl(1)+Lexts+LCs+L+Lextd, xl(1)+Lexts+LCs+L+Lextd];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.05);
xBox = [xl(1)+Lexts+LCs+L+Lextd, xl(1)+Lexts+LCs+L+Lextd, xl(2), xl(2)];
yBox = [yl(1), yl(2), yl(2), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.05);
hold off
title('Potencial electrostático en el grafeno (V_{D}=1 V)')
xlabel('Eje y (m)')
ylabel('\phi (V)')
legend('V_{G}=-2 V','V_{G}=-1.5 V','V_{G}=-1 V','V_{G}=-0.5 V','V_{G}=0 V','V_{G}=0.5 V','V_{G}=1 V','V_{G}=1.5 V','V_{G}=2 V')

%% Funciones
function Rho=RHO(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho0,kB,T,h,vF,vsat,q,m0,Ids1,Ids2,W,nx2)

[~,~,~,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho0,nx2);
[~,~,~,~,Rho]=semic(phiG,kB,T,h,vF,Vs,Vd,dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rho0);

end


function [phi,x,y,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho,nx2)
%Definición del mallado
% x=linspace(0,tt+tg+tb,nx);
% dx=x(2)-x(1);
% nx1=floor(nx*tt/(tt+tg+tb)); %Importante definir número de nodos por cada región
% nx2=floor(nx*tg/(tt+tg+tb));
% nx7=nx-nx1-nx2;
% dx1=dx;
% dx2=dx;
% dx7=dx;
tt=tt*1e9;
tg=tg*1e9;
tb=tb*1e9;
L=L*1e9;
Lexts=Lexts*1e9;
LCs=LCs*1e9;
Lextd=Lextd*1e9;
LCd=LCd*1e9;
rho=rho/1e9^2;



nx1=floor((nx-1)*tt/(tt+tg+tb))+1; %Importante definir número de nodos por cada región
x1=linspace(0,tt,nx1);
x2=linspace(tt,tt+tg,nx2+1);
nx7=nx-nx2-nx1;
x7=linspace(tt+tg,tt+tg+tb,nx7+1);
dx1=x1(2)-x1(1);
dx2=x2(2)-x2(1);
dx7=x7(2)-x7(1);
x=[x1,x2(2:end),x7(2:end)];

y=linspace(0,L+LCs+LCd+Lexts+Lextd,ny);
dy=y(2)-y(1);
ny7=ny;
ny3=floor((ny-1)*Lexts/(L+LCs+LCd+Lexts+Lextd))+1;
ny2=floor((ny-1)*LCs/(L+LCs+LCd+Lexts+Lextd))+1;
ny5=floor((ny-1)*Lextd/(L+LCs+LCd+Lexts+Lextd))+1;
ny6=floor((ny-1)*LCd/(L+LCs+LCd+Lexts+Lextd))+1;
ny4=ny7-ny2-ny3-ny5-ny6;
ny1=ny4;
dy1=dy;
dy2=dy;
dy7=dy;

%Inicialización de matriz
M=zeros(nx1*ny1+nx2*ny7+nx7*ny7);


%R1
for i=1:ny1*nx1%Valores de la matriz M que siguen apartado 4.2.1 del documento. Recorremos
%la matriz completa. Luego redefiniremos aquellos puntos frontera que no
%cumplen esta configuración
    M(i,i)=-2/dy1^2-2/dx1^2;
    M(i,i+1)=1/dy1^2;
    if i>1 %Omitimos valores que exceden índices del array. En este caso es frontera Neuman. 
        M(i,i-1)=1/dy1^2;
    end
    if i>ny1
        M(i,i-ny1)=1/dx1^2;
    end
    if i<ny1*(nx1-1)+1 %Tenemos en cuenta la diferencia de número de nodos en el salto entre la región 1 y 2
        M(i,i+ny1)=1/dx1^2;
    else
        M(i,i+ny1+ny2+ny3)=1/dx1^2;
    end
end

for i=1:ny1 %Redefinimos matriz para fronteras Dirichlet 4.2.2.1
    M(i,:)=0;
    M(i,i)=1;
end


%R2+...+R6
for i=ny1*nx1+1:ny1*nx1+ny7*nx2 %Valores de la matriz M que siguen apartado 4.2.1 del documento.
    M(i,i)=-2*ey/dy2^2-2*ex/dx2^2;
    M(i,i+1)=ey/dy2^2;
    if i>ny1*nx1+1
        M(i,i-1)=ey/dy2^2;
    end
    if i>ny1*nx1+ny7
        M(i,i-ny7)=ex/dx2^2;
    end
    M(i,i+ny7)=ex/dx2^2;
end

for i=ny1*nx1+ny2+ny3+1:ny1*nx1+ny2+ny3+ny4 %Cambio de material TO-->GR. 4.2.3.
    M(i,i-1)=0;
    M(i,i)=et/dx1+ex/dx2; 
    M(i,i+1)=0;
    M(i,i-ny1-ny2-ny3)=-et/dx1;
    M(i,i+ny7)=-ex/dx2;
end
for i=nx1*ny1+1:nx1*ny1+ny2 %Dirichlet
    M(i,:)=0;
    M(i,i)=1;
end
for i=nx1*ny1+ny2+ny3+ny4+ny5:nx1*ny1+ny7 %Dirichlet
    M(i,:)=0;
    M(i,i)=1;
end
for i=ny1*nx1+ny7*nx2-ny7+2:ny1*nx1+ny7*nx2-1 %Cambio de material GR-->BO. 4.2.3.
    M(i,i-1)=0;
    M(i,i)=eb/dx7+ex/dx2;
    M(i,i+1)=0;
    M(i,i+ny7)=-eb/dx7;
    M(i,i-ny7)=-ex/dx2;
end

%R7
for i=ny1*nx1+ny7*nx2+1:ny1*nx1+ny7*nx2+nx7*ny7 %Valores de la matriz M que siguen apartado 4.2.1 del documento.
    M(i,i)=-2/dy7^2-2/dx7^2;
    M(i,i-1)=1/dy7^2;
    if i<ny1*nx1+ny7*nx2+nx7*ny7
        M(i,i+1)=1/dy7^2;
    end
    if i<ny1*nx1+ny7*nx2+nx7*ny7-ny7+1
        M(i,i+ny7)=1/dx7^2;
    end
    M(i,i-ny7)=1/dx7^2;
end

for i=ny1*nx1+ny7*nx2+ny7*nx7-ny7+1:ny1*nx1+ny7*nx2+ny7*nx7 %Dirichlet
    M(i,:)=0;
    M(i,i)=1;
end
%condiciones Neumann (4.2.2.2)
%R1
for i=ny1+1:ny1:ny1*nx1
    %Lado izquierdo
    M(i,:)=0;
    M(i,i)=1;
    M(i,i+1)=-1;
    %Lado derecho
    M(i+ny1-1,:)=0;
    M(i+ny1-1,i+ny1-1)=1;
    M(i+ny1-1,i+ny1-2)=-1;
end
%R2+R3+R4+R5+R6
for i=ny1*nx1+1:ny7:ny1*nx1+ny7*nx2
    %Lado izquierdo
    M(i,:)=0;
    M(i,i)=1;
    M(i,i+1)=-1;
    %Lado derecho
    M(i+ny7-1,:)=0;
    M(i+ny7-1,i+ny7-1)=1;
    M(i+ny7-1,i+ny7-2)=-1;
    
end
%R7
for i=ny1*nx1+ny7*nx2+1:ny7:ny1*nx1+ny7*nx2+nx7*ny7
    %Lado izquierdo
    M(i,:)=0;
    M(i,i)=1;
    M(i,i+1)=-1;
    %Lado derecho
    M(i+ny7-1,:)=0;
    M(i+ny7-1,i+ny7-1)=1;
    M(i+ny7-1,i+ny7-2)=-1;
    
end
%Cuidado, en el bucle for anterior hemos quitado condiciones Dirichlet
%definidas anteriormente. Aquí las reincorporamos. 
M(end-ny7+1,:)=0;
M(end-ny7+1,end-ny7+1)=1;
M(end,:)=0;
M(end,end)=1;
%R3 por arriba
for i=ny1*nx1+ny2+1:ny1*nx1+ny2+ny3
    M(i,:)=0;
    M(i,i)=1;
    M(i,i+ny7)=-1;
end
%R5 por arriba
for i=ny1*nx1+ny2+ny3+ny4+1:ny1*nx1+ny2+ny3+ny4+ny5
    M(i,:)=0;
    M(i,i)=1;
    M(i,i+ny7)=-1;
end


%Inicializo término independiente del sistema
v=zeros(ny1*nx1+ny7*nx2+nx7*ny7,1);
%Dirichlet 4.3.2.1
for i=1:ny1
    v(i)=Vg;
end
%Dirichlet
for i=ny1*nx1+1:ny1*nx1+ny2
    v(i)=Vs;
end
%Dirichlet
for i=ny1*nx1+ny7-ny6+1:ny1*nx1+ny7
    v(i)=Vd;
end
%Dirichlet
for i=ny1*nx1+ny7*nx2+ny7*nx7-ny7+1:ny1*nx1+ny7*nx2+ny7*nx7
    v(i)=Vb;
end
%Interior material 4.3.1
for i=ny1*nx1+1:ny1*nx1+ny7*nx2
    k=floor((i-ny1*nx1-1)/ny7);
    v(i)=v(i)+rho(i-ny7*k-ny1*nx1);
end
%Cuidado, en bucle for anterior definimos valores que siguen apartado 4.3.1
%cuando deben seguir otros apartados. Aquí volvemos a poner valor nulo. 
for i=ny1*nx1+1:ny7:ny1*nx1+ny7*nx2 %Neuman 4.3.2.2 
    v(i)=0;%izquierda
    v(i+ny7-1)=0;%derecha
end
for i=ny1*nx1+ny2+1:ny1*nx1+ny2+ny3%cambio de material 4.3.2.2
    v(i)=0;
end

for i=ny1*nx1+ny2+ny3+ny4+1:ny1*nx1+ny2+ny3+ny4+ny5%Neuman 
    v(i)=0; %arriba R5
end
for i=ny1*nx1+ny2+1:ny1*nx1+ny2+ny3%Neuman 
    v(i)=0; %arriba R3
end
%Resolución ecuación
a=M\v;
%Inicializo matriz solución
phi=zeros(nx,ny);
%Reorganizo valores de a para tener el valor según x e y
for i=1:nx1
    for j=ny2+ny3+1:ny2+ny3+ny4
        phi(i,j)=a(j+(i-1)*ny1-ny2-ny3);
    end
end

for i=nx1+1:nx2+nx1
    for j=1:ny7
        phi(i,j)=a(j+(i-1)*ny7-nx1*(ny7-ny1));
    end
end

for i=nx1+nx2+1:nx2+nx1+nx7
    for j=1:ny7
        phi(i,j)=a(j+(i-1)*ny7-nx1*(ny7-ny1));
    end
end
%Dejamos como valor nulo aquellas posiciones de la matriz que equivalen a posiciones espaciales que
%no pertenecen al dominio. 
phiG=phi(nx1+ceil(nx2/2),:);
x=x*1e-9;
y=y*1e-9;
dy=dy*1e-9;
end


function [V,n,p,m,rho,Ids]=semic(phiG,kB,T,h,vF,Vs,Vd,dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rho)
%ecuación (10)
Ng=2/pi*(kB*T/h/vF)^2;

n=@(V,y) Ng*fermi(1,((-q*V+q*phiG(y))/kB/T));
p=@(V,y) Ng*fermi(1,((q*V-q*phiG(y))/kB/T));

%ecuación (9)
m=zeros(size(phiG));
m(1)=m0/sqrt(1+(abs(m0/vsat*(phiG(2)-phiG(1))/dy))^2); %en los extremos no podemos usar la ecuación (12)
m(end)=m0/sqrt(1+(abs(m0/vsat*(phiG(end)-phiG(end-1))/dy))^2);
for i=2:length(phiG)-1
    m(i)=m0/sqrt(1+(abs(m0/vsat*(phiG(i+1)-phiG(i-1))/2/dy))^2);
end

%% Método disparo
for i=1:10 %con este bucle corregimos los datos iniciales del método de la bisección, en caso de que no cumplan con las condiciones necesarias
    V=Heun(Ids1,W,q,n,p,m,ny2,ny6,ny,dy,Vs);
    
    V1=V(length(phiG)-ny6+1);
    V=Heun(Ids2,W,q,n,p,m,ny2,ny6,ny,dy,Vs);
    
    V2=V(length(phiG)-ny6+1);
    if ((V2-Vd)*(V1-Vd))>0
        fprintf("Datos iniciales incorrectos: %d, %d\n",Ids1, Ids2)
        Ids1=Ids1/2;
        Ids2=Ids2*2;
    else
        %fprintf("Datos bien elegidos\n")
        %i-1
        break
    end
end




for k=1:1000 %con este bucle ejecutamos el método del disparo mediante secante
    Ids=(Ids1/(V1-Vd)-Ids2/(V2-Vd))/(1/(V1-Vd)-1/(V2-Vd));
    if or(Ids<Ids1,Ids>Ids2)
        fprintf("Diverge\n")
    end
    if V2==Vd
        Ids=Ids2;
        V=Heun(Ids,W,q,n,p,m,ny2,ny6,ny,dy,Vs);
        break
    elseif V1==Vd
        Ids=Ids1;
        V=Heun(Ids,W,q,n,p,m,ny2,ny6,ny,dy,Vs);
        break
    end
    V=Heun(Ids,W,q,n,p,m,ny2,ny6,ny,dy,Vs);
    if abs(V(length(phiG)-ny6+1)-Vd)<1e-2 %Error absoluto del método comparado con el valor real conocido Vd
        break
    end
    if V(length(phiG)-ny6+1)<Vd
        Ids1=Ids; 
        V1=V(length(phiG)-ny6+1); 
    elseif V(length(phiG)-ny6+1)>Vd
        Ids2=Ids; 
        V2=V(length(phiG)-ny6+1); 
    end
    if ((V2-Vd)*(V1-Vd))>0
        fprintf("Método de la secante no da resultado: V1=%d, V2=%d\n, Ids1=%d, Ids2=%d, Ids=%d\n", V1,V2,Ids1,Ids2,Ids)  
    end
end

V(ny-ny6+1:ny)=Vd;


for i=1:length(V)
    rho(i)=q*(-n(V(i),i)+p(V(i),i));
end
end

function V=Heun(Ids,W,q,n,p,m,ny2,ny6,ny,dy,Vs)
V=zeros(ny,1); %definimos vector
V(1:ny2)=Vs; %definimos el valor a lo largo de SOURCE, que es conocido. 
for i=ny2:ny-ny6+1 %todos los puntos siguientes a SOURCE, hasta el primer punto de DRAIN
    k1=Ids/W/q/(n(V(i),i)+p(V(i),i))/m(i); %f(y,V)
    k2=Ids/W/q/(n(V(i)+k1*dy,i+1)+p(V(i)+k1*dy,i+1))/m(i+1); %f(y+h,V+k1*h)
    V(i+1)=V(i)+1/2*dy*(k1+k2); 
end
end
%% Fermi
function xx=fermi(aj,eta)
if length(eta)~=1
    error("La dimensión del vector es incorrecta")
end
%Program begins

format long e;
%==============================================================
% Evaluation of Trapezoidal sum begins
range=8.;
if eta > 0.
   range=sqrt(eta+64.);end;
h=0.5;
nmax=range/h;
sum=0.;
if aj== (-0.5)
   sum=1./(1.+exp(-eta));end;
for i=1:nmax
   u=i*h;
   ff=2.*(u^(2.*aj+1))/(1.+exp(u*u-eta));
   sum=sum+ff;end;

%Trapezoidal Summation ends
%==============================================================

% Pole correction for trapezoidal sum begins
pol=0.;
npole=0;
% Fix the starting value of  BK1 to start while loop
bk1=0;
while bk1 <= 14*pi
   npole=npole+1;
   bk=(2*npole -1)*pi;
   %fprintf("%d\n%d\n\n",eta,bk)
   rho=sqrt(eta*eta+bk*bk);
   t1=1;
   t2=0;
   if eta < 0;
      tk=- aj*(atan(-bk/eta)+pi);
   elseif eta ==0;
      tk=0.5*pi*aj;
   else;
      eta > 0;
      tk=aj*atan(bk/eta);
   end;
   rk=- (rho^aj);
   tk=tk+0.5*atan(t2/t1);
   if eta < 0
      rk= -rk;
   end;
   ak=(2.*pi/h)*sqrt(0.5*(rho+eta));
   bk1=(2.*pi/h)*sqrt(0.5*(rho-eta));
   if bk1 <= (14.*pi)
   gama=exp(bk1);
   t1=gama*sin(ak+tk)-sin(tk);
   t2=1.-2.*gama*cos(ak)+gama*gama;
   pol=pol+4.*pi*rk*t1/t2;
   end; %ends if loop above
  end; % Top while loop ends
  npole=npole-1;
  fdp=sum*h+pol;
% Program ends with the following output
xx=fdp/gamma(1+aj);
%   disp('Fermi-Dirac Integral Value');
%   disp(fdp/gamma(1+aj));
%   disp('Number of trapezoidal points & number of poles');
%   disp([round(nmax),npole]);
end
