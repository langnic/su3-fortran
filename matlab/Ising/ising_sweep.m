# Metropolis:
#(i)   Suche rand site i € {1,...,N_xN_y}
#(ii)  Berechne E_old
#(iii) Schlage update s_i,old -> s_i,new vor (Ising: -s_i,old)
#(iv)  Berechne E_new  
#(v)   Akzeptiere s_i,new falls E_new < E_old,
#      akzeptiere s_i,new mit p = exp(-Delta E), falls E_new > E_old, wo Delta E = E_new - E_old
#

# See example under www.durr.itp.unibe.ch/BUW/2017_lattice/durr_ising.m

function [U,acc] = ising_sweep(U, Nx, Ny, b, h, k, do_live_graphics)
  
  acc = 0
  
  for j = 1: Nx*Ny*k
    
    ix = randi(Nx)
    iy = randi(Ny)
    
    imx = mod(ix-2, Nx)+1
    imy = mod(iy-2, Ny)+1
    
    ipx = mod(ix, Nx)+1
    ipy = mod(iy, Ny)+1
    
    Eold = b*U(ix,iy) * ( U(imx,iy) + U(ipx,iy) + U(ix,ipy) + U(ix,imy) ) + U(ix,iy)*h
    Enew = Eold * (-1) #b*U(ix,iy)*(-1) * ( U(imx,iy) + U(ipx,iy) + U(ix,ipy) + U(ix,imy) )
    deltaE = Enew - Eold
    
    #if Enew < Eold      # könnte man weglassen, nur unteres if benutzen
    #  U(ix,iy) = U(ix,iy) * (-1)
    #else
    r = rand(1)
    if r < exp(-deltaE)
      U(ix,iy) = U(ix,iy) * (-1)
      acc = acc + 1
    end
    #end
    
    if do_live_graphics
      image(U.*10)
      drawnow()
    end  
    
  end
  
end

