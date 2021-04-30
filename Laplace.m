
clear all
%diary Laplace
%Suponiendo conocidas las funciones U o Q sobre la frontera.
Ufront= @(x,y)     exp(x)*sin(y); %x*x-y*y;
Qfront= @(x,y)     exp(x)*sin(y)*x+ exp(x)*cos(y)*y; %2*x*x-2*y*y;
%Qfront= @(x,y) x*y;

N=input('Numero de nodos frontera:'); %Número de puntos sobre la frontera que vamos a considerar 12

%DEFINICIÓN DE LOS NODOS SOBRE LA FRONTERA. En este caso son dados por medio de una función pero 
%podriamos incluso meterlos desde una taba deva lores o incluso desde excel.

%En este caso el dominio es una circunferencia
for i=1:N
    cont1(i)=i;
    xfront(i)=cos(2*pi*(i-1)/(N));
    yfront(i)=sin(2*pi*(i-1)/(N));
    Nodo(i,:)=[xfront(i), yfront(i)];
    Nodoz(i,:)=[Nodo(i,:),0]; %para verlos como vectores de 3 componentes
end

Nodo(N+1,:)=Nodo(1,:); %identificamos el nodo 1 como el nodo N+1, esto para poder trabajar sobre el último elemento
Nodoz(N+1,:)=Nodoz(1,:);

for i=1:N
    V(i,:)=Nodo(i+1,:)-Nodo(i,:); %para simplificar un poco la notación más abajo
end

%Calculo de los vectores normales unitarios. Por eso escribo los nodos como con 3 componentes, creo que se puede hacer directo. 
for i=1:N
    normalz(i,:)=cross(Nodoz(i+1,:)-Nodoz(i,:),[0,0,1]); 
    norma(i)=norm(normalz(i,:));
    normaluniz(i,:)=normalz(i,:)/norma(i);
end
normaluniz(:,3)=[]; %elimina la tercera componente, ya sé que todas las terceras componentes son cero
normaluni=normaluniz; %Vuelvo a pasar todo a dos componentes.

%CALCULO DEL PUNTO SOBRE EL CUAL VA A ESTAR CENTRADA LA SOLUCIÓN
%FUNDAMENTAL Y SU DERIVADA. EN ESTE CASO SE ESTÁ TOMANDO EN LOS EXTREMOS DE
%LOS SEGMENTOS LINEALES
for i=1:N
    P(i,:)= Nodo(i,:);
end

%DEFINICIÓN DE LAS MATRICES H-barra Y G
a=-1;
b=1;
for i=1:N
         for j=1:N
             
         if i==1 && j==N
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-3);
         else    
             
            if j==i-1
                barh1(i,j)=0;
                barh2(i,j)=0;
                g1(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-1);
                g2(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-3);
                                
            else
                 
                if j==i
                    barh1(i,j)=0;
                    barh2(i,j)=0;
                    g1(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-3);
                    g2(i,j)=(-1/(8*pi))*norm(V(j,:))*(2*log(norm(V(j,:)))-1);
                   
                else 
                    
                    gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
                    uij1 = @(t) (-1/(8*pi))*norm(V(j,:))*(1-t)*log(norm(gammaj(t)-P(i,:)));
                    qij1 = @(t) (-1/(8*pi))*norm(V(j,:))*(1-t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
                    uij2 = @(t) (-1/(8*pi))*norm(V(j,:))*(1+t)*log(norm(gammaj(t)-P(i,:)));
                    qij2 = @(t) (-1/(8*pi))*norm(V(j,:))*(1+t)*dot(gammaj(t)-P(i,:),normaluni(j,:))/((norm(gammaj(t)-P(i,:)))^2);
            
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


%barh1;
%barh2;

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
barH;
G;

for i=1:N
    c(i)=-sum(barH(i,:));
end

idn=diag(c);
H=idn+barH; %CREACION DE LA MATRIZ H


%Insertar los datos conocidos en la frontera y despejar para poner la ecuaciónn de la forma AX=Y


k=input('Digite 0 ó 2 para las condiciones de frontera Tipo Dirichlet y Mixta :');

 

%PROBLEMA TIPO DIRICHLET-> K=0 SI CONOCEMOS TODOS LOS VALORES DE U SOBRE LA FRONTERA.

if k==0
    A=G;
    for i=1:N
        Ufrontaux(i)=Ufront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
    end
    Ufrontaux=Ufrontaux';
    Y0=H*Ufrontaux;
    X=A\Y0;
    Ufrontera=Ufrontaux;
    Qfrontera=X;    
end

% %PROBLEMA TIPO NEUMANN-> K=1 SI CONOCEMOS TODOS LOS VALORES DE Q (DERIVADA NORMAL DE U) SOBRE LA FRONTERA. 
% if k==1
%     A=H;
%     for i=1:N
%         Qfrontaux(i)=Qfront(P(i,1),P(i,2)); %Creacion de los vectores auxiliares con los datos conocidos.
%     end
%     Qfrontaux=Qfrontaux';
%     Y1=G*Qfrontaux;
%     X1=pinv(A)*Y1; %X=A\Y; 
%     Ufrontera=X1;
%     Qfrontera=Qfrontaux;
% end


%PROBLEMA MIXTO-> K=2 SI CONOCEMOS UNA PARTE DE U Y Q SOBRE LA FRONTERA. 

if k==2
    
    r=input('Cuantos nodos Ufront conoce:');
    
    for i=1:r
        T(i)=Ufront(P(i,1),P(i,2)); 
        cont3(i)=i;
    end

    for i=r+1:N
        T(i)=Qfront(P(i,1),P(i,2));
        cont4(i-r)=i;
    end
    
    T=T';
    
    A1=[H(:,1:r) -G(:,r+1:N)];
    A2=[G(:,1:r) -H(:,r+1:N)];
    Y2=A1*T;
    X2=A2\Y2;
    
    for i=1:r
        
        Ufrontera(i)=Ufront(P(i,1),P(i,2));
        Qfrontera(i)=X2(i);
    end
    
    for i=r+1:N
        
        Qfrontera(i)=Qfront(P(i,1),P(i,2));
        Ufrontera(i)=X2(i);
    end
    
   Ufrontera=Ufrontera';
   Qfrontera=Qfrontera';
end


%Se calcula el valor exacto

for i=1:N
   UfrontExac(i)=Ufront(P(i,1),P(i,2));
end
for i=1:N
   QfrontExac(i)=Qfront(P(i,1),P(i,2));
end


%Se calcula el error numerico
ErrorUfront=norm(Ufrontera'-UfrontExac)/norm(UfrontExac) 
ErrorQfront=norm(Qfrontera'-QfrontExac)/norm(QfrontExac) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CÁLCULO DE U EN LOS PUNTOS INTERIORES
NI=23; %Número de nodos interiores a usar
%Elección de los puntos interiores 23
for i=1:10
    xint(i)=0.5*cos(2*pi*(i-1)/10);
    yint(i)=0.5*sin(2*pi*(i-1)/10);
    Pint(i,:)= [xint(i),yint(i)];
end

for i=1:10
    xint(i)=0.2*cos(2*pi*(i-1)/10);
    yint(i)=0.2*sin(2*pi*(i-1)/10);
    Pint(10+i,:)= [xint(i),yint(i)];
end
Pint(23,:)=[-0.1,0];
Pint(22,:)=[0,0];
Pint(21,:)=[0.1,0];



            %CÁLCULO DE H Y G EN LOS PUNTOS INTERIORES

for i=1:NI
    for j=1:N
        
             gammaj = @(t) ((Nodo(j+1,:)+Nodo(j,:))+t*(Nodo(j+1,:)-Nodo(j,:)))/2; %Esta es la parametrización. Es una recta que va de nodo a nodo
             uij1int = @(t) (-1/(8*pi))*norm(V(j,:))*(1-t)*log(norm(gammaj(t)-Pint(i,:)));
             qij1int = @(t) (-1/(8*pi))*norm(V(j,:))*(1-t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
             uij2int = @(t) (-1/(8*pi))*norm(V(j,:))*(1+t)*log(norm(gammaj(t)-Pint(i,:)));
             qij2int = @(t) (-1/(8*pi))*norm(V(j,:))*(1+t)*dot(gammaj(t)-Pint(i,:),normaluni(j,:))/((norm(gammaj(t)-Pint(i,:)))^2);
                                   
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

Uint=Gint*Qfrontera-barHint*Ufrontera;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:NI
    cont2(i)=N+i;
    UintExac(i)=Ufront(Pint(i,1),Pint(i,2));
end

ErrorUint=norm(Uint'-UintExac)/norm(UintExac) %%%%%%Para imprimir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hatP=[P;Pint];
bigU=[Ufrontera;Uint];


%IMPRESION DE DATOS
tabla1=[cont1',Ufrontera,UfrontExac',Qfrontera,QfrontExac']';
tabla2=[cont2',Uint,UintExac']';

fprintf('%10s %12s %12s %12s %12s\n','Nodo','Ufrontera','UfrontExac','Qfrontera','QfrontExac');
fprintf('%8.0f %12.4f %12.4f %12.4f %12.4f\n',tabla1);
disp('--------------------------------------------------')
fprintf('%6s %11s %12s\n','Nodo',' Uint','UintExac');
fprintf('%6.0f %12.4f %12.4f\n',tabla2);

%GRÁFICA

xlin = linspace(min(hatP(:,1)),max(hatP(:,1)),33); 
ylin = linspace(min(hatP(:,2)),max(hatP(:,2)),33);
[X,Y] = meshgrid(xlin,ylin); %creación de la malla uniforme
Z = griddata(hatP(:,1),hatP(:,2),bigU,X,Y,'linear');
figure
mesh(X,Y,Z) %Grafica de la malla uniforme con los valores interpolados
axis tight; hold on
plot3(hatP(:,1),hatP(:,2),bigU,'.','MarkerSize',15) %esto grafica los puntos originales (SOLO LOS PUNTOS)
title('Solución numérica de U utilizando 50 nodos frontera y 43 nodos interiores')










