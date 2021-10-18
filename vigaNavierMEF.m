##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO ELEMENTOS FINITOS (VIGA PORTICO DE NAVIER)
## /_____/_/  |_/___/___/___/     |
##                                |    metodoM: script principal
##---------CICLO LECTIVO 2020----------------------------------------------------

## CONFIGURACION

clear

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

syms x

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

t=1.27e-3;

A1=t*2*10e-3;

A2=t*6e-3;

E=70.36e9;

L=100e-3;

Fe=17.27/4;

Fs=.9*85.8375/4;

Fg=.1*85.8375/2;

Ipc=2.20108e-10; # Momento de inercia desde centroide

xc=0.00296549; # Centroide

f=1;

#Ip=Ipc+A1*(f*10e-3-xc)^2;

Ip=Ipc;

Is=2.28600e-11;

##---- CARACTERISTICAS FISICO-GEOMETRICAS

NODO=[0 0;0 L;L L;L 0;L/2 L/2]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E A1 Ip;2 3 E A1 Ipc;4 3 E A1 Ip;4 1 E A1 Ipc;1 5 E A2 Is;2 5 E A2 Is;3 5 E A2 Is;4 5 E A2 Is]; % [NodoInicial NodoFinal E A I] en fila i define ubicacion de la viga "i-esima" y sus propiedades

#%--------------CARGAS ---> L O C A L E S

cargasLocales=[1 -Fs/L 0;3 -Fs/L 0]; # Nelemento Qx Qy

#########%  CONDICIONES DE CONTORNO Y CARGAS ---> G L O B A L E S

%-------- COND./CARGAS NULAS SE DEJAN V A C I A S = [] -------------#


CCx=[1 0;4 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 0;4 0]; % [Nodo Uy] define condicion de contorno en nodo i 

CCw=[]; % [Nodo W] define condicion de contorno en nodo i 

CARGAx=[]; % [Nodo Px] define carga en nodo i

CARGAy=[2 -Fe;3 -Fe;5 -Fg]; % [Nodo Py] define carga en nodo i

MOMENTO=[]; % [Nodo M] define MOMENTO en nodo i ( P O S I T I V O sentido anti-horario)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SCRIPT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cargaLocal]=proyeccionCargas(ELEMENTO,NODO,cargasLocales);

[KG,fq]=rigidezGlobal(ELEMENTO,NODO,cargaLocal);

GL=size(KG,1); % Cant. de grados de libertad globales

[P,U]=vectorCargas(GL,CARGAx,CARGAy,MOMENTO,CCx,CCy,CCw);

				% Guyan

INDEX=linspace(1,GL,GL);

i=1;guyan=[1 1];

CC=sum(isnan(U));

				#while (isempty(guyan)<1 && i<CC+1)
while (i<CC+1)
  
  Null=INDEX(isnan(U));

  Zeros=find(U==0);

  if (Null(i)>Zeros(1))

    guyan=[Null(i) Zeros(1)];

    try

      registroGuyan=[registroGuyan;guyan];

    catch

      registroGuyan=[guyan];

    end_try_catch

    [KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,guyan,1,2);
        
  endif

  i++;
  
endwhile


				% RESOLUCION



K11=KG(1:CC,1:CC);

K12=KG(1:CC,CC+1:GL);

K21=KG(CC+1:GL,1:CC);

K22=KG(CC+1:GL,CC+1:GL);

PII=P(1:CC,:); % Cargas conocidas

Fl=fq(1:CC,:); % Cargas locales asociadas a P conocidas

Fp=fq(CC+1:end,:); % Cargas locales NO asociadas a P conocidas

UI=U(CC+1:GL,:); % Desplazamientos conocidos

determinante=det(K11);

if determinante==0

  UII=zeros(size(PII));

else
  
  UII=inv(K11)*(PII-Fl-K12*UI);  % Desplazamientos desconocidos

endif


PI=K21*(UII)+K22*UI+Fp; % Cargas desconocidas


				% REARMADO DE RESULTADOS

j=1;h=j;i=j;

while (j<=size(P,1))

  
  
  if isnan(P(j))==1

    P(j)=PI(i);

    i++;
    
  endif

  if isnan(U(j))==1

    U(j)=UII(h);

    h++;
    
  endif

  j++;
  
endwhile


[KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,registroGuyan,2,1);


##---------- VERIFICACION y POST-PROCESADO ---------------

GLY=zeros(1,GL);  % Grados de libertad globales discriminados

GLY(linspace(2,GL-1,GL/3))=1;

GLW=mod(linspace(1,GL,GL),3)==0;   

GLX=mod(linspace(1,GL,GL),1)==0-GLW-GLY;


postProcesado(ELEMENTO,NODO,U,10000)


				% VERIFICACION

FX=sum(P(GLX==1))-sum(fq(GLX==1)) % sumatoria en X

FY=sum(P(GLY==1))-sum(fq(GLY==1)) % sumatoria en Y

