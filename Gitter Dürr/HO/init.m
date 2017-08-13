function X = init(Nt, do_random_init)

  if do_random_init
    X = rand(1,Nt).-0.5
  else
    X = zeros(1,Nt)
  end  

end