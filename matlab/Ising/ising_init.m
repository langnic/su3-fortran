function U = ising_init(Nx, Ny, do_random_init)

  if do_random_init
    U(:,:) = randi(2,Nx,Ny).*2.-3
  else
    U(:,:) = ones(Nx, Ny);
  end  

end