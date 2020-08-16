##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO VIGA-PORTICO DE NAVIER)
## /_____/_/  |_/___/___/___/     |
##                                |    vectorCargas: Conformado del Vector de cargas y CC
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [P,U]=vectorCargas(GL,CARGAx,CARGAy,MOMENTO,CCx,CCy,CCw)
  
  
  P=NaN(GL,1);

  for i=1:size(CARGAx,1)

    if (isempty(CARGAx)==0)
      
      GLx=3*CARGAx(i,1)-2;

      P(GLx)=CARGAx(i,2);

    endif
    
  endfor

  for i=1:size(CARGAy,1)

    if (isempty(CARGAy)==0)
      
      GLy=3*CARGAy(i,1)-1;

      P(GLy)=CARGAy(i,2);

    endif

  endfor

  for i=1:size(MOMENTO,1)

    if (isempty(MOMENTO)==0)

      GLw=3*MOMENTO(i,1);

      P(GLw)=MOMENTO(i,2);

    endif
    
  endfor



  
				% VECTOR DE DESPLAZAMIENTOS
  U=NaN(GL,1);

  
  for i=1:size(CCx,1)

    if (isempty(CCx)==0)

      GLx=3*CCx(i,1)-2;

      U(GLx)=CCx(i,2);

    endif
    

  endfor

  for i=1:size(CCy,1)

    if (isempty(CCy)==0)
      
      GLy=3*CCy(i,1)-1;

      U(GLy)=CCy(i,2);

    endif
    
  endfor

  for i=1:size(CCw,1)

    if (isempty(CCw)==0)
    
      GLw=3*CCy(i,1);

      U(GLw)=CCw(i,2);

    endif
    
  endfor

endfunction

  
