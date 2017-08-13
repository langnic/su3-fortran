# Metropolis:
#(i)   Suche rand site i € {1,...,N_xN_y}
#(ii)  Berechne E_old
#(iii) Schlage update s_i,old -> s_i,new vor (Ising: -s_i,old)
#(iv)  Berechne E_new  
#(v)   Akzeptiere s_i,new falls E_new < E_old,
#      akzeptiere s_i,new mit p = exp(-Delta E), falls E_new > E_old, wo Delta E = E_new - E_old
#

# www.durr.itp.unibe.ch/BUW/2017_lattice/durr_qmosc.m

function [Acc,E,x1,x2,x3,x4,h] = main(Nt, o, m, a, do_random_init, do_thermalize)
  
  X = init(Nt, do_random_init)
  acc = 0
  h = 1
  
  # Thermalisierung
  if do_thermalize
    [X,acc,h]= sweep(X, Nt, o, m, h, 100, false);
  end  
  
  for i = 1 : 1000
    [X,acc,h]= sweep(X, Nt, o, m, h, 10, false);
    
    # measure:
    
    # acceptance  
    Acc(i) = acc/(10*Nt)

    # energy
    [Ekin, Epot] = getE(X,Nt,m,o)
    E(i) = Ekin + Epot
    
    # moments
    x1(i) = mean(X)
    x2(i) = mean(X.^2)
    x3(i) = mean(X.^3)
    x4(i) = mean(X.^4)
    

    
  end
  
  image(X.*10)
  drawnow()
  
end

