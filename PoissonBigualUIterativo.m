clear all


%EJEMPLO DE FUNCIONES CONOCIDAS PARA USARLAS SOBRE LA FRONTERA Y PARA
%COMPARAR CON LAS CALCULADAS CON EL PROGRAMA
Ufront= @(x,y)  sin(x);   %Función U
Qfront= @(x,y)  cos(x)*x; %Derivada normal de U sobre la frontera. Este ejemplo es suponiendo que la fontera es una circunferencia unitaria. Si la frontera cambia esta formula cambia pues el vector normal cambia

%FUNCIÓN DE ENTRADA
%funb= @(x,y) x+y;

%FAMILIA DE SOLUCIONES PARTICULARES A USAR (U sombrero).
funf= @(x,y) 1+norm([x,y]);
hatU= @(x,y) (norm([x,y])^2)/4+(norm([x,y])^3)/9; %Si usamos f(r)=1+r
GradU= @(x,y) (1/2+norm([x,y])/3)*[x,y];

N=input('Numero de nodos frontera:');%NÚMERO DE PUNTOS SOBRE LA FRONTERA QUE SE VAN A CONSIDERAR 50

%DEFINICIÓN DE LOS NODOS SOBRE LA FRONTERA: En este caso son dados por medio de una función pero podriamos incluso meterlos desde una taba de valores o incluso desde excel (hasta donde sé eso se puede hacer)
%En este caso el dominio es una circunferencia.
for i=1:N
    xfront(i)=cos(2*pi*(i-1)/N);
    yfront(i)=sin(2*pi*(i-1)/N);
    Nodo(i,:)=[xfront(i), yfront(i)];
    Nodoz(i,:)=[Nodo(i,:),0]; %Para verlos como vectores de 3 componentes. La tercera componente se define como cero. 
end
Nodo(N+1,:)=Nodo(1,:); %Identificamos el nodo 1 como el nodo N+1, esto para poder trabajar sobre el último elemento
Nodoz(N+1,:)=Nodoz(1,:);

for i=1:N
    L(i,:)=Nodo(i+1,:)-Nodo(i,:); %para simplificar un poco la notación más abajo
end

%Calculo de los vectores normales unitarios. Por escribimos los nodos como con 3 componentes y hacemos producto cruz con el vector [0,0,1]. 
for i=1:N
    normalz(i,:)=cross(Nodoz(i+1,:)-Nodoz(i,:),[0,0,1]); 
    norma(i)=norm(normalz(i,:));
    normaluniz(i,:)=normalz(i,:)/norma(i);
    end
normaluniz(:,3)=[]; %Elimina la tercera componente. Sabemos que todas las terceras componentes son cero
normaluni=normaluniz; %Vuelvo a pasar todo a dos componentes.

%CALCULO DEL PUNTO SOBRE EL CUAL VA A ESTAR CENTRADA LA SOLUCIÓN
%FUNDAMENTAL Y SU DERIVADA NORMAL. 
%SE ESTÁ TOMANDO EN EL MEDIO DEL SEGMENTO LINEAL
for i=1:N
    P(i,:)= Nodo(i,:);
end

%DEFINICIÓN DE NODOS INTERIORES
NI=43; %Número de nodos interiores a usar
%Elección de los puntos interiores
for i=1:10
    xint(i)=0.5*cos(2*pi*(i-1)/10);
    yint(i)=0.5*sin(2*pi*(i-1)/10);
    Pint(i,:)= [xint(i),yint(i)];
end

for i=1:10
    xint(i)=0.1*cos(2*pi*(i-1)/10);
    yint(i)=0.1*sin(2*pi*(i-1)/10);
    Pint(10+i,:)= [xint(i),yint(i)];
end

for i=1:10
    xint(i)=0.7*cos(2*pi*(i-1)/10);
    yint(i)=0.7*sin(2*pi*(i-1)/10);
    Pint(20+i,:)= [xint(i),yint(i)];
end

for i=1:10
    xint(i)=0.3*cos(2*pi*(i-1)/10);
    yint(i)=0.3*sin(2*pi*(i-1)/10);
    Pint(30+i,:)= [xint(i),yint(i)];
end


Pint(43,:)=[-0.05,0];
Pint(42,:)=[0,0];
Pint(41,:)=[0.05,0];

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%DEFINO PUNTOS AUXILIARES PARA USAR CON LAS FUNCIONES hatU y hatQ
for i=1:N+NI
    if 1<=i & i<=N
        hatP(i,:)=P(i,:);
    else
        hatP(i,:)=Pint(i-N,:);
    end
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%DEFINICIÓN DE LAS MATRICES H-barra Y G
a=-1;
b=1;
for i=1:N
         for j=1:N
             
         if i==1 && j==N
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
         else    
             if j==i-1
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
                                
            else
                if j==i
                    barh1(i,j)=0;
                    barh2(i,j)=0;
                    g1(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-3);
                    g2(i,j)=(-1/(8*pi))*norm(L(j,:))*(2*log(norm(L(j,:)))-1);
                   
                else 
                    
                    gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
                    uij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-P(i,:)));
                    qij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
                    uij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-P(i,:)));
                    qij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
            
                    %INTEGRACIÓN POR CUADRATURA DE GAUSS
                    [t,w]=Gausslp(4);
                    x=((b-a)*t+a+b)/2;
                    M=length(x);
                   for k=1:M
                        h1aux(k,:)=qij1(x(k));
                        g1aux(k,:)=uij1(x(k));
                        h2aux(k,:)=qij2(x(k));
                        g2aux(k,:)=uij2(x(k));
                   end
                   barh1(i,j)=w*h1aux*(b-a)/2;
                   barh2(i,j)=w*h2aux*(b-a)/2; %Resultado final del proceso de integración
                   g1(i,j)=w*g1aux*(b-a)/2; 
                   g2(i,j)=w*g2aux*(b-a)/2;    %Resultado final del proceso de integración
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
         end                
        end
end

for i=1:N
    for j=1:N
        if j==1
            barH(i,1)=barh1(i,1)+barh2(i,N);
            G(i,1)=g1(i,1)+g2(i,N);
        else
            barH(i,j)=barh1(i,j)+barh2(i,j-1);
            G(i,j)=g1(i,j)+g2(i,j-1); 
        end
    end
end

%CÁLCULO INDIRECTO DE LOS COEFICIENTES CORRESPONDIENTES A LAS ESQUINAS DEL
%DOMINIO
for i=1:N
    c(i)=-sum(barH(i,:));
end

idn=diag(c);

H=idn+barH;

%CÁLCULO DE U EN LOS PUNTOS INTERIORES
%Calculo De H y G en los puntos interiores
a=-1;
b=1;
for i=1:NI
    for j=1:N
        
             gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
             uij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-Pint(i,:)));
             qij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
             uij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-Pint(i,:)));
             qij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
                                   
            %INTEGRACIÓN POR CUADRATURA DE GAUSS
             [t,w]=Gausslp(4);
              x=((b-a)*t+a+b)/2;
              M=length(x);
            for k=1:M
                h1auxint(k,:)=qij1int(x(k));
                g1auxint(k,:)=uij1int(x(k));
                h2auxint(k,:)=qij2int(x(k));
                g2auxint(k,:)=uij2int(x(k));
                            
            end
            barh1int(i,j)=w*h1auxint*(b-a)/2;
            barh2int(i,j)=w*h2auxint*(b-a)/2; %Resultado final del proceso de integración
            g1int(i,j)=w*g1auxint*(b-a)/2; 
            g2int(i,j)=w*g2auxint*(b-a)/2;    %Resultado final del proceso de integración
            
       
    end
end

for i=1:NI
    for j=1:N
        if j==1
            barHint(i,j)=barh1int(i,j)+barh2int(i,N);
            Gint(i,j)=g1int(i,j)+g2int(i,N);
        else
            
            barHint(i,j)=barh1int(i,j)+barh2int(i,j-1);
            Gint(i,j)=g1int(i,j)+g2int(i,j-1); 
             
        end
    end
end

%%%%%%%%%%%%
%%%%%%%%%%%%
%DEFINICIÓN DEFINICIÓN DE LOS VECTORES hatU, hatQ, Y DE LAS MATRICES hatMatU 
%y hatMatQ
for k=1:N+NI
    for j=1:N
        hatUprom(j,k)=hatU(P(j,1)-hatP(k,1),P(j,2)-hatP(k,2));
        hatQprom(j,k)=dot(GradU(P(j,1)-hatP(k,1),P(j,2)-hatP(k,2)),normaluni(j,:));
    end
end

for k=1:N+NI
    for j=1:NI
        hatUpromint(j,k)=hatU(Pint(j,1)-hatP(k,1),Pint(j,2)-hatP(k,2));
    end
end

%Calculo de la matrix F
for i=1:N+NI
    for j=1:N+NI
        MatF(i,j)=funf(hatP(i,1)-hatP(j,1),hatP(i,2)-hatP(j,2));
    end
end


k=input('Digite 0 ó 2 para las condiciones de frontera Tipo Dirichlet y Mixta :');


%PROBLEMA TIPO DIRICHLET-> K=0 SI CONOCEMOS TODOS LOS VALORES DE U SOBRE LA FRONTERA.

if k==0
    A=G;
    for i=1:N
        Ufrontaux(i)=Ufront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
    end
    Ufrontaux=Ufrontaux'; 
end

% %PROBLEMA TIPO NEUMANN-> K=1 SI CONOCEMOS TODOS LOS VALORES DE Q (DERIVADA NORMAL DE U) SOBRE LA FRONTERA. 
% if k==1
%     A=H;
%     for i=1:N
%         Qfrontaux(i)=Qfront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
%     end
%     Qfrontaux=Qfrontaux';
% end

%PROBLEMA MIXTO-> K=2 SI CONOCEMOS UNA PARTE DE U Y Q SOBRE LA FRONTERA. 
if k==2
    r=input('Cuantos nodos Ufront conoce:');
    
    for i=1:r
        T(i)=Ufront(P(i,1),P(i,2)); 
    end

    for i=r+1:N
        T(i)=Qfront(P(i,1),P(i,2));
    end
    T=T';
end

vecbold=zeros(1,N+NI);
vecbold=vecbold';%Para volverlo columna
tol=2;
Iteraciones=0;
while tol>0.0001
       
%Cálculo de \alpha
vecalpha=MatF\vecbold;
%Cálculo vector d
vecd=(H*hatUprom-G*hatQprom)*vecalpha;
%INSERTAR LOS DATOS NOCOCIDOS EN LA FRONTERA Y DESPEJAR PARA PONER LA
%ECUACION DE LA FORMA AX=Y.

if k==0
    Y0=H*Ufrontaux-vecd;
    X0=A\Y0;
    Ufrontera=Ufrontaux;
    Qfrontera=X0;
end
if k==1
    Y1=vecd+G*Qfrontaux;
    X1=pinv(A)*Y1; %X=A\Y; 
    Ufrontera=X1;
    Qfrontera=Qfrontaux;
end
if k==2
    A1=[H(:,1:r) -G(:,r+1:N)];
    A2=[G(:,1:r) -H(:,r+1:N)];
    Y2=A1*T;
    X2=A2\Y2;
    for i=1:r
        Ufrontera(i)=Ufront(P(i,1),P(i,2));
        Qfrontera(i)=X2(i);
        Ufrontera=Ufrontera';
        Qfrontera=Qfrontera';
    end
    
    for i=r+1:N  
        Qfrontera(i)=Qfront(P(i,1),P(i,2));
        Ufrontera(i)=X2(i);
        Ufrontera=Ufrontera';
        Qfrontera=Qfrontera';
    end
end

vecdint=(barHint*hatUprom-Gint*hatQprom+hatUpromint)*vecalpha; %%%%%% Revisar

Uint=Gint*Qfrontera-barHint*Ufrontera+vecdint;
bigU=[Ufrontera;Uint];
tol=norm(bigU+vecbold)/norm(bigU);
vecbold=-bigU;

Iteraciones=Iteraciones+1;
end

Iteraciones


for i=1:N
   Ufrontexac(i)=Ufront(P(i,1),P(i,2));
end
for i=1:NI
   UintExac(i)=Ufront(Pint(i,1),Pint(i,2));
end
bigUexact=[Ufrontexac,UintExac]';

bigU-bigUexact;
ErrorbigU=norm(bigU-bigUexact)/norm(bigUexact)



%GRÁFICA
%Esto es para gráficar. Para que se vea bonito como una superficie primero
%se debe crear una malla uniformemente espaciada. Los valores desconocidos
%se calculan interpolando los valores conocidos. 

xlin = linspace(min(hatP(:,1)),max(hatP(:,1)),33); 
ylin = linspace(min(hatP(:,2)),max(hatP(:,2)),33);
[X,Y] = meshgrid(xlin,ylin); %creación de la malla uniforme
Z = griddata(hatP(:,1),hatP(:,2),bigU,X,Y,'linear');
figure
mesh(X,Y,Z) %Grafica de la malla uniforme con los valores interpolados
axis tight; hold on
plot3(hatP(:,1),hatP(:,2),bigU,'.','MarkerSize',15) %esto grafica los puntos originales (SOLO LOS PUNTOS)