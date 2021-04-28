clear all


%EJEMPLO DE FUNCIONES CONOCIDAS PARA USARLAS SOBRE LA FRONTERA Y PARA
%COMPARAR CON LAS CALCULADAS CON EL PROGRAMA
Ufront= @(x,y)  exp(-x);   %Funci�n U
Qfront= @(x,y)  -exp(-x)*x; %Derivada normal de U sobre la frontera. Este ejemplo es suponiendo que la fontera es una circunferencia unitaria. Si la frontera cambia esta formula cambia pues el vector normal cambia



%FUNCI�N DE ENTRADA
%funb= -du/dx;

%FAMILIA DE SOLUCIONES PARTICULARES A USAR (U sombrero).
funf= @(x,y) 1+(norm([x,y]))^2 +(norm([x,y]))^3; % %EN ESTE PROBLEMA ES NECESARIO CAMBIAR LA FUNCI�N f PARA EVITAR SINGULARIDADES AL MOMENTO DE DERIVARLA
fundfdx= @(x,y) 2*x+3*x*norm([x,y]); 
hatU= @(x,y) ((norm([x,y]))^2)/4+((norm([x,y]))^4)/16 +((norm([x,y]))^5)/25; % %RECUERDE QUE DEPENDE DE LA f ESCOGIDA
GradU= @(x,y) (1/2+((norm([x,y]))^2)/4+((norm([x,y]))^3)/5)*[x,y];

N=input('Numero de nodos frontera:');%N�MERO DE PUNTOS SOBRE LA FRONTERA QUE SE VAN A CONSIDERAR 40

%DEFINICI�N DE LOS NODOS SOBRE LA FRONTERA: En este caso son dados por medio de una funci�n pero podriamos incluso meterlos desde una taba de valores o incluso desde excel (hasta donde s� eso se puede hacer)
%En este caso el dominio es una circunferencia.
for i=1:N
    xfront(i)=cos(2*pi*(i-1)/N);
    yfront(i)=sin(2*pi*(i-1)/N);
    Nodo(i,:)=[xfront(i), yfront(i)];
    Nodoz(i,:)=[Nodo(i,:),0]; %Para verlos como vectores de 3 componentes. La tercera componente se define como cero. 
end
Nodo(N+1,:)=Nodo(1,:); %Identificamos el nodo 1 como el nodo N+1, esto para poder trabajar sobre el �ltimo elemento
Nodoz(N+1,:)=Nodoz(1,:);

for i=1:N
    L(i,:)=Nodo(i+1,:)-Nodo(i,:); %para simplificar un poco la notaci�n m�s abajo
end


%Calculo de los vectores normales unitarios. Por escribimos los nodos como con 3 componentes y hacemos producto cruz con el vector [0,0,1].
for i=1:N
    normalz(i,:)=cross(Nodoz(i+1,:)-Nodoz(i,:),[0,0,1]); 
    norma(i)=norm(normalz(i,:));
    normaluniz(i,:)=normalz(i,:)/norma(i);
    end
normaluniz(:,3)=[]; %Elimina la tercera componente. Sabemos que todas las terceras componentes son cero
normaluni=normaluniz; %Vuelvo a pasar todo a dos componentes.

%CALCULO DEL PUNTO SOBRE EL CUAL VA A ESTAR CENTRADA LA SOLUCI�N FUNDAMENTAL Y SU DERIVADA NORMAL. SE EST� TOMANDO EN LOS EXTREMOS DE LOS SEGMENTOS
for i=1:N
    P(i,:)= Nodo(i,:);
end

%DEFINICI�N DE NODOS INTERIORES
NI=43; %N�mero de nodos interiores a usar
%Elecci�n de los puntos interiores
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
%DEFINO PUNTOS AUXILIARES PARA USAR CON LAS FUNCIONES hatU y hatQ. Los
%llamo hatP, son los de la frontera, seguidos de los interiores
for i=1:N+NI
    if 1<=i & i<=N
        hatP(i,:)=P(i,:);
    else
        hatP(i,:)=Pint(i-N,:);
    end
end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%DEFINICI�N DE LAS MATRICES H-barra Y G
a=-1;
b=1;
for i=1:N
         for j=1:N
             
         if i==1 & j==N
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
                    
                    gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrizaci�n. Es una recta que va de nodo a nodo
                    uij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-P(i,:)));
                    qij1 = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
                    uij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-P(i,:)));
                    qij2 = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
            
                    %INTEGRACI�N POR CUADRATURA DE GAUSS
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
                   barh2(i,j)=w*h2aux*(b-a)/2; %Resultado final del proceso de integraci�n
                   g1(i,j)=w*g1aux*(b-a)/2; 
                   g2(i,j)=w*g2aux*(b-a)/2;    %Resultado final del proceso de integraci�n
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
         end                
        end
end

for i=1:N
    for j=1:N
        if j==1
            barHfron(i,1)=barh1(i,1)+barh2(i,N);
            Gfron(i,1)=g1(i,1)+g2(i,N);
        else
            barHfron(i,j)=barh1(i,j)+barh2(i,j-1);
            Gfron(i,j)=g1(i,j)+g2(i,j-1); 
        end
    end
end

%C�LCULO INDIRECTO DE LOS COEFICIENTES CORRESPONDIENTES A LAS ESQUINAS DEL
%DOMINIO
for i=1:N
    c(i)=-sum(barHfron(i,:));
end

idn=diag(c);
Hfron=idn+barHfron;

%C�lculo De H y G en los puntos interiores
a=-1;
b=1;
for i=1:NI
    for j=1:N
        
             gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrizaci�n. Es una recta que va de nodo a nodo
             uij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*log(norm(gammaj(t)-Pint(i,:)));
             qij1int = @(t) (-1/(8*pi))*norm(L(j,:))*(1-t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
             uij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*log(norm(gammaj(t)-Pint(i,:)));
             qij2int = @(t) (-1/(8*pi))*norm(L(j,:))*(1+t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
                                   
            %INTEGRACI�N POR CUADRATURA DE GAUSS
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
            barh2int(i,j)=w*h2auxint*(b-a)/2; %Resultado final del proceso de integraci�n
            g1int(i,j)=w*g1auxint*(b-a)/2; 
            g2int(i,j)=w*g2auxint*(b-a)/2;    %Resultado final del proceso de integraci�n
            
       
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

%DEFINICI�N DE LAS NUEVAS H Y G AUMENTADAS

bigH=[Hfron,zeros(N,NI);barHint,eye(NI)];
bigG=[Gfron,zeros(N,NI);Gint,zeros(NI,NI)];

%DEFINICI�N DEFINICI�N DE LOS VECTORES hatU, hatQ, Y DE LAS MATRICES hatMatU 
%y hatMatQ

for k=1:N+NI
    for j=1:N
        hatUfron(j,k)=hatU(P(j,1)-hatP(k,1),P(j,2)-hatP(k,2));
        hatQfron(j,k)=dot(GradU(P(j,1)-hatP(k,1),P(j,2)-hatP(k,2)),normaluni(j,:));
    end
end

for k=1:N+NI
    for j=1:NI
        hatUint(j,k)=hatU(Pint(j,1)-hatP(k,1),Pint(j,2)-hatP(k,2));
        hatQint(j,k)=0;
    end
end

%DEFINCI�N DE LOS VECTORES hatU Y hatQ AUMENTADOS

bighatU=[hatUfron;hatUint];
bighatQ=[hatQfron;hatQint];

%C�lculo de la matrix F
for i=1:N+NI
    for j=1:N+NI
        MatF(i,j)=funf(hatP(i,1)-hatP(j,1),hatP(i,2)-hatP(j,2));
    end
end

%C�lculo de la matrix \partial{F}\partial{x}
for i=1:N+NI
    for j=1:N+NI
        MatdFdx(i,j)=fundfdx(hatP(i,1)-hatP(j,1),hatP(i,2)-hatP(j,2));
    end
end


%CALCULO DE LA MATRIZ bigS

bigS=[bigH*bighatU-bigG*bighatQ]*inv(MatF);
bigR=bigS*MatdFdx*inv(MatF);
Aaux=bigH+bigR;

%INGRESO DE VALORES CONOCIDOS EN LA FRONTERA

k=input('Digite 0 � 2 para las condiciones de frontera Tipo Dirichlet y Mixta :');
%SOLO CONSIDERAMOS DOS CASOS: 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k=0 SI CONOCEMOS TODOS LOS VALORES DE U 
if k==0
    for i=1:N
        Ufrontaux(i)=Ufront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
    end
    Ufrontaux=[Ufrontaux,zeros(1,NI)];
    Ufrontaux=Ufrontaux';
    Y=-Aaux*Ufrontaux;
    
    A1=Aaux(:,N+1:N+NI);
    A2=-bigG(:,1:N);
    A=[A2 A1];
    %X=pinv(A)*Y;
    X=A\Y;
            
    bigU=[Ufrontaux(1:N);X(N+1:N+NI)];
    bigQ=[X(1:N);zeros(NI,1)];
end

% PROBLEMA TIPO NEUMANN
% if k==1
%      A=Aaux;
%     for i=1:N
%         Qfrontaux(i)=Qfront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
%     end
%     Qfrontaux=[Qfrontaux,zeros(1,NI)];
%     Qfrontaux=Qfrontaux';
%     Y=bigG*Qfrontaux;
%     %X=pinv(A)*Y; 
%     X=A\Y; 
%     bigU=X;
%     bigQ=Qfrontaux;
% end

    %PROBLEMA MIXTO-> K=2 SI CONOCEMOS UNA PARTE DE U Y Q SOBRE LA FRONTERA.
    %Una buena aproximaci�n es N=50 Y 25 nodos Ufron conocidos
    
    if k==2
        
        B=Aaux;
        T=input('Cuantos nodos Ufront conoce:');
    
        for i=1:T
            V(i)=Ufront(P(i,1),P(i,2)); 
        end

        for i=T+1:N
            V(i)=Qfront(P(i,1),P(i,2));
        end
    
        V=[V,zeros(1,NI)];
        bigV=V';
        A2=[-bigG(:,1:N-T) B(:,N-T+1:N+NI)];        
        Y2=[-B(:,1:T) bigG(:,T+1:N+NI)]*bigV;       
        X2=A2\Y2;
        Qfrontera1=X2(1:T); %VALORES ENCONTRADOS
        Qfrontera=[Qfrontera1', V(T+1:N)];
        bigU1=X2(T+1:N+NI); %VALORES ENCONTRADOS 
        bigU=[V(1:T),bigU1'];
        bigU=bigU';
    end
    
for i=1:N
   Ufrontexac(i)=Ufront(P(i,1),P(i,2));
end
for i=1:NI
   UintExac(i)=Ufront(Pint(i,1),Pint(i,2));
end
bigUexact=[Ufrontexac,UintExac]';




ErrorbigU=norm(bigU-bigUexact)/norm(bigUexact)

%GR�FICA
%Esto es para gr�ficar. Para que se vea bonito como una superficie primero
%se debe crear una malla uniformemente espaciada. Los valores desconocidos
%se calculan interpolando los valores conocidos.

xlin = linspace(min(hatP(:,1)),max(hatP(:,1)),33); 
ylin = linspace(min(hatP(:,2)),max(hatP(:,2)),33);
[X,Y] = meshgrid(xlin,ylin); %creaci�n de la malla uniforme
Z = griddata(hatP(:,1),hatP(:,2),bigU,X,Y,'linear');
figure
mesh(X,Y,Z) %Grafica de la malla uniforme con los valores interpolados
axis tight; hold on
plot3(hatP(:,1),hatP(:,2),bigU,'.','MarkerSize',15) %esto grafica los puntos originales (SOLO LOS PUNTOS)

