##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO ELEMENTOS FINITOS (VIGA PORTICO DE NAVIER)
## /_____/_/  |_/___/___/___/     |
##                                |    proyeccionCargas: Proyeccion de las cargas distribuidas locales
##                                |                     en los grados de libertad 
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [cargaLocal]=proyeccionCargas(ELEMENTO,NODO,cargasLocales)

  if (isempty(cargasLocales)==0)
        
    syms x
    
    Nv= @(x,L) [0;(2*x^3)/L^3-(3*x^2)/L^2+1;x^3/L^2-(2*x^2)/L+x;0;(3*x^2)/L^2-(2*x^3)/L^3;x^3/L^2-x^2/L];

    Nb= @(x,L) [1-x/L;0;0;x/L;0;0];

    
    for u=1:size(cargasLocales,1)

				% CALCULOS GEOMETRICOS POR ELEMENTO

      elem=double(cargasLocales(u,1));
      
      nodoI=ELEMENTO(elem,1);

      nodoF=ELEMENTO(elem,2);

      Ufx=NODO(nodoF,1);
      
      Ufy=NODO(nodoF,2);

      Uix=NODO(nodoI,1);
      
      Uiy=NODO(nodoI,2);

      ladoY=Ufy-Uiy;

      ladoX=Ufx-Uix;
      
      Largo=sqrt(ladoY^2+ladoX^2);

				# ---------- Proyecciones --------#

      barra=-integrador(Nb(x,Largo)*cargasLocales(u,2),0,Largo);
      
      viga=-integrador(Nv(x,Largo)*cargasLocales(u,3),0,Largo);

    				# ---------- Ensamblado --------#

      try

	cargaLocal=[cargaLocal;elem,(viga+barra)'];

      catch

	cargaLocal=[elem,(viga+barra)'];

      end_try_catch
      
      
    endfor

  else

    cargaLocal=zeros(1,7);

  endif
  
endfunction
