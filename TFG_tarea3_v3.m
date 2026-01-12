clear
%% Definición de los datos del sistema
%Para ecuación de Poisson
tt=30e-9; %en metros
tg=0.6e-9; 
tb=300e-9;
L=500e-9;
LCs=50e-9;
Lexts=100e-9;
Lextd=100e-9;
LCd=50e-9;
Vs=0; %Voltios
Vg=1;
Vd=0;
Vb=0;
e0=8.8541878176e-12; %F/m
et=5.9*e0; %TO
eb=3.9*e0; %BO
ey=3.3*e0; %GR y
ex=9*e0; %GR x
ny=70; 
nx2=7; %número de nodos en el eje x en la lámina de grafeno.


%Para ecuación diferencial física de semiconductores
h=6.626e-34/2/pi; %J/s
kB=1.38e-23; %J/K
vF=1e6; %m/s
q=1.6e-19; %C
T=300; %K

m0=0.2; %m^2/(V*s)
vsat=1e6; %m/s
W=100e-9; %m

%Intensidad de corriente, límites para la resolución del método de la
%secante para resolver la ecuación de física de semiconductores. 
Ids1=0; %menor (A)
Ids2=1; %mayor (A)


K=1e-1;
%% Discretización de nuestro dominio. Número de nodos en cada dirección. 
nx=ceil((nx2-1)*0.01*(tt+tg+tb)/tg)+1; %garantizamos nx2>7. 
%En este caso, el tamaño en la dirección de los nodos en el grafeno son 100
%veces más pequeños que fuera de él. Para ello el factor 0,01.

%ny ya está definido en los datos del problema. 

%% fmincon
% f=@(rho)fobjt(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho*1e-8,kB,T,h,vF,vsat,q,m0,Ids1,Ids2,W)*1e13;

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% 
% tic
% rho=fmincon(f,rho*1e-5,[],[],[],[],[],[],[],options)*1e-5;
% toc
%% fminsearch 
% options = optimset('Display','iter');
% tic
% rho=fminsearch(f,rho*1e8,options)*1e-8; %no utiliza gradiente
% toc

%% fsolve 
%Se definen todas las posibles situaciones de carga que se van a plantear en la simulación
% Vds=0:0.025:1;
% Vgs=[-1,0,1,2];

Vgs=-2:0.5:2;
Vds=0:0.05:1;


%Inicializamos algunas variables

%rho es la matriz con los resultados de las optimizaciones resueltas con 
%fsolve para las distintas polaridades (densidad de carga)

%rho1 es la matriz con las distribuciones de carga calculadas con rho para 
%las distintas polaridades (densidad de carga, ha de ser igual que rho, con
%cierto error admisible)

%ef es el "exit flag" de fsolve. Indica la bondad del resultado de fsolve

%Ids es la matriz con las intensidades que atraviesan el transistor para
%cada polarización



rho=zeros(ny,length(Vds),length(Vgs));
rho1=zeros(ny,length(Vds),length(Vgs));
Ids=zeros(length(Vds),length(Vgs));



%Establecemos las opciones de optimización. Ciertos parámetros, como el
%TypicalX, son vitales para la convergencia del modelo. 

tic %Comenzamos a contar el tiempo de ejecución
for j=1:length(Vgs)
    
    for i=1:length(Vds) %Doble bucle para barrer todas las distribuciones de carga
        fprintf('%d %d\n',Vgs(j),Vds(i))
        error=1e10;
        %Semilla de fsolve. Valor inicial de la optimización.
        if i==1 %Aquí se evitan problemas indiciales (índices menores que 1)
            rh=rho(:,i,j);
        else
            rh=rho(:,i-1,j);
        end
        
        %Método de punto fijo
        for l=1:10
            rho(:,i,j)=RHO(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rh,kB,T,h,vF,vsat,q,m0,Ids1,Ids2,W,nx2);
            if sum(abs(rho(:,i,j)-rh))*1e5<3.23e-1
                break
            else
                if sum(abs(rho(:,i,j)-rh))<error
                    error=sum(abs(rho(:,i,j)-rh));
                    fprintf("%d\n\n",error)
                    rh=rho(:,i,j);
                else
                    fprintf("Punto fijo diverge\n")
                    break
                end
                
                
                
                
            end
        end
        [phi,x,y,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vds(i),Vgs(j),Vb,et,eb,ex,ey,rho(:,i,j),nx2);
        [~,~,~,~,rho1(:,i,j),Ids(i,j)]=semic(phiG,kB,T,h,vF,Vs,Vds(i),dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rho(:,i,j));
       

    end
end
toc %Muestra por pantalla el tiempo total de ejecución

%%


%Representación solución
figure(1)
plot(Vds,Ids(:,1))
hold on
for i=2:length(Vgs)
    plot(Vds,Ids(:,i))
end
hold off
legend('V_{G}=-2 V','V_{G}=-1.5 V','V_{G}=-1 V','V_{G}=-0.5 V','V_{G}=0 V','V_{G}=0.5 V','V_{G}=1 V','V_{G}=1.5 V','V_{G}=2 V')
xlabel('V_{D} (V)')
ylabel('I (A)')

%% Guardar variables

save('variables','rho','rho1','Ids')
%% F

function Rho=RHO(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho0,kB,T,h,vF,vsat,q,m0,Ids1,Ids2,W,nx2)

[~,~,~,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho0,nx2);
[~,~,~,~,Rho]=semic(phiG,kB,T,h,vF,Vs,Vd,dy,ny2,ny6,vsat,q,m0,Ids1,Ids2,W,ny,rho0);

end


function [phi,x,y,phiG,dy,ny2,ny6]=poisson(nx,ny,tt,tg,tb,L,Lexts,LCs,Lextd,LCd,Vs,Vd,Vg,Vb,et,eb,ex,ey,rho,nx2)

%Para evitar problemas de redondeo, se cambian las unidades de longitud m-->nm

tt=tt*1e9;
tg=tg*1e9;
tb=tb*1e9;
L=L*1e9;
Lexts=Lexts*1e9;
LCs=LCs*1e9;
Lextd=Lextd*1e9;
LCd=LCd*1e9;
rho=rho/1e9^2;



%Definición del mallado

nx1=floor((nx-1)*tt/(tt+tg+tb))+1; %Importante definir número de nodos por cada región. Determino el número de intervalos (nx-1) totales y los reparto proporcionalmente a la longitud que ocupan. Sumo uno para transformarlo al número de nodos.
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
%Inicializo matriz solución para tener los valores de forma bidimensional
phi=zeros(nx,ny);
%Reorganizo valores de 'a' para tener el valor según x e y
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


%Método Runge-Kutta de Heun aplicado al problema. 
function V=Heun(Ids,W,q,n,p,m,ny2,ny6,ny,dy,Vs)
V=zeros(ny,1); %definimos vector
V(1:ny2)=Vs; %definimos el valor a lo largo de SOURCE, que es conocido. 
for i=ny2:ny-ny6+1 %todos los puntos siguientes a SOURCE, hasta el primer punto de DRAIN
    k1=Ids/W/q/(n(V(i),i)+p(V(i),i))/m(i); %f(y,V)
    k2=Ids/W/q/(n(V(i)+k1*dy,i+1)+p(V(i)+k1*dy,i+1))/m(i+1); %f(y+h,V+k1*h)
    V(i+1)=V(i)+1/2*dy*(k1+k2); 
end
end








%% Fermi, usado en la función semic(..)
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