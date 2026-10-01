

# Measure Magnetisation

function E = ising_getE(U,Nx,Ny,b,h)
  
  E = 0
  
  for ix = 1:Nx
    for iy = 1:Ny
    
      imx = mod(ix-2, Nx)+1
      imy = mod(iy-2, Ny)+1
      
      ipx = mod(ix, Nx)+1
      ipy = mod(iy, Ny)+1
      E += b*U(ix,iy) * ( U(imx,iy) + U(ipx,iy) + U(ix,ipy) + U(ix,imy) ) + U(ix,iy)*h
      
    end
  end
  E = E/2
  
end