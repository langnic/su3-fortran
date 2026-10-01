# Metropolis:
#(i)   Suche rand site i € {1,...,N_xN_y}
#(ii)  Berechne E_old
#(iii) Schlage update s_i,old -> s_i,new vor (Ising: -s_i,old)
#(iv)  Berechne E_new  
#(v)   Akzeptiere s_i,new falls E_new < E_old,
#      akzeptiere s_i,new mit p = exp(-Delta E), falls E_new > E_old, wo Delta E = E_new - E_old
#

# See example under www.durr.itp.unibe.ch/BUW/2017_lattice/durr_ising.m
# See example under www.durr.itp.unibe.ch/BUW/2017_lattice/wrap_ising.m

function [Acc,M,E] = ising_main(Nx, Ny, b, h, do_random_init, do_thermalize)
  
  U(:,:) = ising_init(Nx, Ny, do_random_init)
  acc = 0
  
  # Thermalisierung
  if do_thermalize
    [U,acc]= ising_sweep(U,Nx, Ny, b, h, 10, false);
    acc = 0
  end  
  
  for i = 1 : 100
    [U,acc]= ising_sweep(U,Nx, Ny, b, h, 10, false);
    
    # measure
    Acc(i) = acc/(10*Nx*Ny)
    acc = 0
    
    M(i) = sum(U(:))/(Nx*Ny);
    E(i) = ising_getE(U,Nx,Ny,b,h)/(Nx*Ny)
    
  end
  
  image(U.*10)
  drawnow()
  
end

