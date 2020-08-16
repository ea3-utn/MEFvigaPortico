##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO ELEMENTOS FINITOS (VIGA PORTICO DE NAVIER)
## /_____/_/  |_/___/___/___/     |
##                                |    rigidezGlobal: Calculo de Matriz de Rig. Global y local
##                                |                     Globalizacion de Cargas locales
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [KG,fq]=rigidezGlobal(ELEMENTO,NODO,cargaLocal)

  
  k= @(E,A,L,I) [(A*E)/L 0 0 -(A*E)/L 0 0;0 (12*E*I)/L^3 (6*E*I)/L^2 0 -(12*E*I)/L^3 (6*E*I)/L^2;0 (6*E*I)/L^2 (4*E*I)/L 0 -(6*E*I)/L^2 (2*E*I)/L;-(A*E)/L 0 0 (A*E)/L 0 0;0 -(12*E*I)/L^3 -(6*E*I)/L^2 0 (12*E*I)/L^3 -(6*E*I)/L^2;0 (6*E*I)/L^2 (2*E*I)/L 0 -(6*E*I)/L^2 (4*E*I)/L];

  T= @(cos,sin) [cos sin 0 0 0 0;-sin cos 0 0 0 0;0 0 1 0 0 0;0 0 0 cos sin 0;0 0 0 -sin cos 0;0 0 0 0 0 1];
  
  for u=1:size(ELEMENTO,1)

				% CALCULOS GEOMETRICOS POR ELEMENTO

    nodoI=ELEMENTO(u,1);

    nodoF=ELEMENTO(u,2);

    Ufx=NODO(nodoF,1);
    
    Ufy=NODO(nodoF,2);

    Uix=NODO(nodoI,1);
    
    Uiy=NODO(nodoI,2);

    ladoY=Ufy-Uiy;

    ladoX=Ufx-Uix;
    
    Largo=sqrt(ladoY^2+ladoX^2);

    coseno=ladoX/Largo;

    seno=ladoY/Largo;
    
    
    

				% MATRIZ LOCAL A GLOBAL
   
    
    conectividad=[3*nodoI-2 3*nodoI-1 3*nodoI 3*nodoF-2 3*nodoF-1 3*nodoF];
    
    K=transpose(T(coseno,seno))*k(ELEMENTO(u,3),ELEMENTO(u,4),Largo,ELEMENTO(u,5))*T(coseno,seno); % Matriz local en cuatro grados de libertad
    
    fqGlobal=transpose(T(coseno,seno))*transpose(cargaLocal(find(cargaLocal(:,1)==u),2:7));
    

    ################# ENSAMBLADO #####################################
    
    for i=1:6

      if (isempty(fqGlobal)==1)

	try

	  fq(conectividad(i))=fq(conectividad(i))+0;
	  
	catch

	  fq(conectividad(i))=0;
	  
	end_try_catch
	
      else
	
	try

	  fq(conectividad(i))=fq(conectividad(i))+fqGlobal(i);
	  
	catch

	  fq(conectividad(i))=fqGlobal(i);
	  
	end_try_catch

     endif


      
      for j=1:6

	try

	  KG(conectividad(i),conectividad(j))=KG(conectividad(i),conectividad(j))+K(i,j);
	  
	catch

	  KG(conectividad(i),conectividad(j))=K(i,j);
	  
	end_try_catch

      endfor
    endfor

  endfor

  fq=transpose(fq);
  
endfunction

  
