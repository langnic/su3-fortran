
function [Ekin, Epot] = getE(X,Nt,m,o)
  Ektmp = 0; 
  Eptmp = 0;
  for j = 1: Nt
    jm = mod(j-2, Nt)+1;    
    jp = mod(j, Nt)+1;   
    Ektmp += m/2 * (X(jp) - X(j))^2;
    Eptmp += m*o^2/2 * X(j)^2;
  end
  Ekin = Ektmp/Nt;
  Epot = Eptmp/Nt;
end  
    
    