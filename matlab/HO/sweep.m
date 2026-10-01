# Metropolis:
#(i)   Suche rand site i € {1,...,N_xN_y}
#(ii)  Berechne E_old
#(iii) Schlage update s_i,old -> s_i,new vor (Ising: -s_i,old)
#(iv)  Berechne E_new  
#(v)   Akzeptiere s_i,new falls E_new < E_old,
#      akzeptiere s_i,new mit p = exp(-Delta E), falls E_new > E_old, wo Delta E = E_new - E_old
#

function [X,acc,h] = sweep(X, Nt, o, m, h, k, do_live_graphics)
  
  acc = 0
  
  for l = 1: k
    for j = 1: Nt
      
      # calculate indices
      it = randi(Nt)
      imt = mod(it-2, Nt)+1    
      ipt = mod(it, Nt)+1
      
      # calculate shift
      xshift = h * (rand - 0.5)
      xprop = X(it) + xshift
      
      # calculate old and new energy
      Eold = m/2 * ( (X(ipt) - X(it))^2 + ( ( X(it) - X(imt) )^2 ) ) + m*o^2/2 * X(it)^2 
      Enew = m/2 * ( (X(ipt) - xprop)^2 + ( ( xprop - X(imt) )^2 ) ) + m*o^2/2 * xprop^2 
      deltaE = Enew - Eold
      
      # update if accepted
      r = rand(1)
      if r < exp(-deltaE)
        X(it) = xprop    
        acc = acc + 1
      end 
      
      # plot X if life_graphics
      if do_live_graphics
        image(X.*10)
        drawnow()
      end  
      
    end
    
    # adjust h after one sweep
    h = max(h*acc/(l*Nt*0.8), 1e6);
  end

  acc = acc/(Nt*k)  
  
end

