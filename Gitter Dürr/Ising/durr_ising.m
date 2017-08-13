function [cen,err,hist,susc,U]=durr_ising(Nx,Ny,b,h,n_ther,n_meas,n_sepa,U)

   b_crit=log(1+sqrt(2))/2; %%% b_crit=0.44069
   if nargin<1, Nx=       64, end
   if nargin<2, Ny=       64, end
   if nargin<3, b=    b_crit, end
   if nargin<4, h=         0, end
   if nargin<5, n_ther=  100, end
   if nargin<6, n_meas= 1000, end
   if nargin<7, n_sepa=  100, end
   if nargin<8
      do_randinit=real(b<b_crit), U=ising_initialize(Nx,Ny,do_randinit);
   else
      [Nx,Ny]=size(U), do_randinit=nan; %%% note: overwrite [Nx,Ny]
   end
   figure(01); spy(U+ones(Nx,Ny),'k')

   for loop_ther=1:n_ther
      for loop_sepa=1:n_sepa
         [U,acc]=ising_metropolis(U,b,h);
      end; clear loop_sepa
      fprintf('%6i: acc=%6.4f\n',loop_ther-n_ther,acc)
   end; clear loop_ther acc
   figure(02); spy(U+ones(Nx,Ny),'c')

   hist_accept=nan(n_meas,1); %%% note: initialize as column vector of nan
   hist_energy=nan(n_meas,1); %%% note: initialize as column vector of nan
   hist_magnet=nan(n_meas,1); %%% note: initialize as column vector of nan
   for loop_meas=1:n_meas
     %U=ising_flipflop(U,b,h);
      for loop_sepa=1:n_sepa
         [U,acc]=ising_metropolis(U,b,h);
      end; clear loop_sepa ; hist_accept(loop_meas)=acc;
      e=ising_energy(U,b,h); hist_energy(loop_meas)=e;
      m=mean(mean(U))      ; hist_magnet(loop_meas)=m;
      fprintf('%6i: acc=%6.4f, energy=%8.4f, magnet=%8.4f\n',loop_meas,acc,e,m)
   end; clear loop_meas acc e m
  %orig_energy=Nx*Ny*(hist_energy.^2-mean(hist_energy)^2); %%% large variance
  %orig_magnet=01*01*(hist_magnet.^2-mean(hist_magnet)^2); %%% large variance
   susc_energy=Nx*Ny*(hist_energy-mean(hist_energy)).^2;
   susc_magnet=01*01*(hist_magnet-mean(hist_magnet)).^2;
   susc_magimp=01*01*(hist_magnet-real(h==0)/real(h==0)+1).^2;
   cumu_binder=1-mean(hist_magnet.^4)/3/mean(hist_magnet.^2)^2; %%% ~0 or ~2/3

   hist=[hist_accept,hist_energy,hist_magnet];
   susc=[susc_energy,susc_magnet,susc_magimp];
   cen=mean([hist,susc]);             %%% note: works along rows --> [1,6]
   err=std([hist,susc])/sqrt(n_meas); %%% note: works along rows --> [1,6]
   figure(03); spy(U+ones(Nx,Ny),'m')
   figure(04); plot([1:n_meas],hist_energy,'b.'); legend('hist\_energy')
   figure(05); plot([1:n_meas],hist_magnet,'r.'); legend('hist\_magnet')
   figure(06); plot([1:n_meas],susc_energy,'b.'); legend('susc\_energy')
   figure(07); plot([1:n_meas],susc_magnet,'r.'); legend('susc\_magnet')
   figure(08); plot([1:n_meas],susc_magimp,'r.'); legend('susc\_magimp')
   figure(09); histogram(hist_energy); legend(sprintf('hist\\_energy, b=%g',b))
   figure(10); histogram(hist_magnet); legend(sprintf('hist\\_magnet, b=%g',b))
   figure(11); histogram(susc_energy); legend(sprintf('susc\\_energy, b=%g',b))
   figure(12); histogram(susc_magnet); legend(sprintf('susc\\_magnet, b=%g',b))
   figure(13); histogram(susc_magimp); legend(sprintf('susc\\_magimp, b=%g',b))

   fprintf('Nx=%i, Ny=%i, b=%g, h=%g, do_randinit=%i\n',Nx,Ny,b,h,do_randinit)
   fprintf('average of hist_accept = %9.6f +/- %8.6f\n',cen(1),err(1))
   fprintf('average of hist_energy = %9.6f +/- %8.6f\n',cen(2),err(2))
   fprintf('average of hist_magnet = %9.6f +/- %8.6f\n',cen(3),err(3))
   fprintf('average of susc_energy = %9.6f +/- %8.6f\n',cen(4),err(4))
   fprintf('average of susc_magnet = %9.6f +/- %8.6f\n',cen(5),err(5))
   fprintf('average of susc_magimp = %9.6f +/- %8.6f\n',cen(6),err(6))
   exact=(1-1/sinh(2*b)^4)^0.125; %%% note: Onsager solution for magnetization
   fprintf('exact_magnet^2= %8.6f, cumu_binder= %8.6f\n',exact^2,cumu_binder)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U]=ising_initialize(Nx,Ny,do_randinit)
   if do_randinit
      U=2*randi(2,Nx,Ny)-3; %%% spin configuration with random +1/-1
   else
      U=ones(Nx,Ny);        %%% spin configuration with all spins +1
   end
end

function [U,acc]=ising_metropolis(U,b,h)
   [Nx,Ny]=size(U);
   table_i=randi(Nx,[Nx*Ny,1]); %%% lookup table for i values (x-coordinates)
   table_j=randi(Ny,[Nx*Ny,1]); %%% lookup table for j values (y-coordinates)
   table_r=rand([Nx*Ny,1]);     %%% lookup table will be used in acc/rej-test
   if 0 %%% note: switch to overwrite with random permutation of site numbers
      obj=randperm(Nx*Ny);         %%% random permutation of [1:Nx*Ny]
      table_i=mod(obj-1,Nx)+1;     %%% compute pertinent x-coordinates
      table_j=floor((obj-1)/Nx)+1; %%% compute pertinent y-coordinates
      rec=(table_j-1)*Nx+table_i;  %%% checking for correct (j-1)*Nx+i
      if (nnz(rec-obj)>0), error('ising_metropolis: randperm was misused'), end
   end
   acc=0;
   for touch=1:Nx*Ny
      i=table_i(touch); i_plu=mod(i,Nx)+1; i_min=mod(i-2,Nx)+1;
      j=table_j(touch); j_plu=mod(j,Ny)+1; j_min=mod(j-2,Ny)+1;
      neighborsum=U(i_plu,j)+U(i_min,j)+U(i,j_plu)+U(i,j_min);
      e_old=-b*U(i,j)*neighborsum-h*U(i,j);
     %e_new=+b*U(i,j)*neighborsum+h*U(i,j); %%% note: generic would-be energy
      e_new=-e_old;                         %%% note: shortcut in ising model
      boltz=exp(-(e_new-e_old));
      if table_r(touch)<boltz
         U(i,j)=-U(i,j); %%% note: perform accepted spin flip
         acc=acc+1;      %%% note: step-up acceptance counter
      end
   end; clear touch
   acc=acc/(Nx*Ny); %%% acceptance counter --> acceptance ratio
end

function [U,acc]=ising_flipflop(U,b,h)
   acc=0;
   e=ising_energy(U,b,h);         %%% current energy
   U_new=-U;                      %%% after potential flip
   e_new=ising_energy(U_new,b,h); %%% after potential flip
   if rand<exp(-(e_new-e))
      U=U_new;
      e=e_new;
      acc=acc+1;
   end
end

function [e]=ising_energy(U,b,h)
   if (1)
      [Nx,Ny]=size(U);
      e=0;
      for j=1:Ny; j_plu=mod(j,Ny)+1;
      for i=1:Nx; i_plu=mod(i,Nx)+1;
         e=e-b*U(i,j)*(U(i_plu,j)+U(i,j_plu))-h*U(i,j);
      end; clear i
      end; clear j
      e=e/(Nx*Ny); %%% normalize per site
   else
      V=U.*circshift(U,-1,1); %%% note: multiply with spin in position x+1
      W=U.*circshift(U,-1,2); %%% note: multiply with spin in position y+1
      e=-b*mean(V(:))-b*mean(W(:))-h*mean(U(:)); %%% mean(U(:))=mean(mean(U))
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng('default'); [cen,err]=durr_ising(32,32,0.4) [OHNE flipflop, 15.2 sec]
%      0.32345     -0.44245    -0.015674      0.95351     0.061515     0.061761
%   0.00096636   0.00096545    0.0078471     0.044873    0.0024611    0.0024676
%
% rng('default'); [cen,err]=durr_ising(32,32,0.4) [ MIT flipflop, 15.3 sec]
%      0.32303     -0.44304    0.0056113      0.89813     0.061735     0.061767
%   0.00093818   0.00093699    0.0078611     0.039864    0.0024077    0.0024069
%
% rng('default'); [cen,err]=durr_ising(32,32,0.5) [OHNE flipflop, 15.2 sec]
%     0.079635     -0.87183      0.91051      0.73224   0.00075537      0.82978
%   0.00053193   0.00084604   0.00086955     0.035276   5.4131e-05    0.0015555
%
% rng('default'); [cen,err]=durr_ising(32,32,0.5) [ MIT flipflop, 15.3 sec]
%     0.080172     -0.87206  -0.00029883      0.80984      0.82963      0.82963
%    0.0005443   0.00088975     0.028818     0.035969    0.0016684    0.0016684
