##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO ELEMENTOS FINITOS (VIGA PORTICO DE NAVIER)
## /_____/_/  |_/___/___/___/     |
##                                |    integrador: Integrador Numerico, metodo trapezoidal
##                                |                     
##---------CICLO LECTIVO 2020----------------------------------------------------------------


function d=integrador(z,xi,xf)
  
  PRESICION=10;
  
  dom=linspace(xi,xf,PRESICION);
  
  d=zeros(size(z));
  
  for i=1:size(z,1)

    for j=1:size(z,2)	     

      test=sym(z(i,j));
      
      if (isempty(symvar(test))==1)

	y=ones(size(dom))*double(z(i,j)); # Si el integrando es constante

      else
      
	fp=function_handle(z(i,j));
	
	y=fp(dom);
	
      endif  
	
	
      try
	
	d(i,j)=trapz(dom,y);

      catch  # Si es nulo

	d(i,j)=0;

      end_try_catch

      
    endfor
    
  endfor
  

endfunction

