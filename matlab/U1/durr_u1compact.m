function []=durr_u1compact(Nx,Ny,Nz,Nt,be,n_ther,n_meas,n_sepa,n_hits,n_over,U)
% usage: []=durr_u1compact(Nx,Ny,Nz,Nt,be,n_ther,n_meas,n_sepa,n_hits,n_over,U)

   if nargin< 1, Nx=       8, end
   if nargin< 2, Ny=       8, end
   if nargin< 3, Nz=       8, end
   if nargin< 4, Nt=       8, end
   if nargin< 5, be=     1.1, end
   if nargin< 6, n_ther= 100/10, end %%% number of configs to be skipped
   if nargin< 7, n_meas=1000/10, end %%% number of configs for measurements
   if nargin< 8, n_sepa=  10/10, end %%% number of packages between configs
   if nargin< 9, n_hits=   8, end %%% number of hits in MH-multihit update
   if nargin<10, n_over=   4, end %%% OR-to-MH sweep ratio in one package
   if nargin<11, U=complex(ones(4,Nx,Ny,Nz,Nt));
      elseif  5~=ndim(U)  , error('durr_u1compact:  5~=ndim(U)'),
      elseif  4~=size(U,1), error('durr_u1compact:  4~=size(U,1)'),
      elseif Nx~=size(U,2), error('durr_u1compact: Nx~=size(U,2)'),
      elseif Ny~=size(U,3), error('durr_u1compact: Ny~=size(U,3)'),
      elseif Nz~=size(U,4), error('durr_u1compact: Nz~=size(U,4)'),
      elseif Nt~=size(U,5), error('durr_u1compact: Nt~=size(U,5)'),
   end
   figure(1); [~,~,~,s_loc]=action_wils(U);
   if max(s_loc(:)-mean(s_loc(:)))>0
      contourf(squeeze(max(max(s_loc)))); axis ij; axis equal; clear s_loc
   else
      contourf([1:Nz]'*[1:Nt],[-1.0,-1.0]);  axis ij; axis equal;
   end
   title('max(max(s\_loc)) after initialization'); xlabel('t/a'); ylabel('z/a')

   hitsize=1;
   for loop_ther=1:n_ther
      for loop_sepa=1:n_sepa
         [U,acc_mh,hitsize]=u1compact_metropolis(U,be,hitsize,n_hits);
      end; clear loop_sepa
      U=U./abs(U); %%% note: backproject to U(1) occasionally
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f\n',...
       loop_ther-n_ther,hitsize,acc_mh)
   end; clear loop_ther acc_mh
   figure(2); [~,~,~,s_loc]=action_wils(U);
    contourf(squeeze(max(max(s_loc)))); axis ij; axis equal; clear s_loc
    title('max(max(s\_loc)) after thermalization'); xlabel('t/a');ylabel('z/a')

   s_avg=nan(n_meas,1); s_spa=nan(n_meas,1); s_tem=nan(n_meas,1);
   corr_spaspa=nan(n_meas,Nt); corr_temtem=nan(n_meas,Nt);
   for loop_meas=1:n_meas
      %%% new configuration means n_sepa packages of n_over*OR+1*MH sweeps
      for loop_sepa=1:n_sepa; acc_or=nan(1,n_over);
         for loop_over=1:n_over
            [U,acc_or(loop_over)]=u1compact_overrelax(U,be);
         end; clear loop_over
         [U,acc_mh,hitsize]=u1compact_metropolis(U,be,hitsize,n_hits);
      end; clear loop_sepa
      U=U./abs(U); %%% note: backproject to U(1) occasionally
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f, acc_or=%6.4f\n',...
       loop_meas,hitsize,acc_mh,mean(acc_or))
      %%% start measurements on configuration number loop_meas
      [obj_avg,obj_spa,obj_tem]=action_wils(U);
      s_avg(loop_meas)=obj_avg;
      s_spa(loop_meas)=mean(obj_spa);
      s_tem(loop_meas)=mean(obj_tem);
      [obj_spa,~,~,obj_tem]=calc_corr(obj_spa,obj_tem);
      corr_spaspa(loop_meas,:)=obj_spa;
      corr_temtem(loop_meas,:)=obj_tem;
   end; clear loop_meas acc_mh acc_or obj_avg obj_spa obj_tem
   figure(3); [~,~,~,s_loc]=action_wils(U);
    contourf(squeeze(max(max(s_loc)))); axis ij; axis equal; clear s_loc
    title('max(max(s\_loc)) after finishing'); xlabel('t/a'); ylabel('z/a')

   fprintf('s_spa: %9.6f +/- %8.6f\n',mean(s_spa),std(s_spa)/sqrt(n_meas))
   fprintf('s_tem: %9.6f +/- %8.6f\n',mean(s_tem),std(s_tem)/sqrt(n_meas))
   fprintf('s_avg: %9.6f +/- %8.6f\n',mean(s_avg),std(s_avg)/sqrt(n_meas))

   figure(4); plot(s_avg); title('history of << s\_avg >>')

   corr_cen=mean(corr_spaspa);              %%% note: avg along 1st dim --> row
   corr_err=std (corr_spaspa)/sqrt(n_meas); %%% note: std along 1st dim --> row
   corr_cen(find(corr_cen<0))=nan;          %%% note: overwrite neg.entries nan
   corr_err( isnan(corr_cen))=nan;          %%% note: overwrite correspond. nan
   corr_cen(Nt/2:Nt) % corr_cen([1:5,Nt-5:Nt])
   corr_err(Nt/2:Nt) % corr_err([1:5,Nt-5:Nt])
   figure(5); errorbar([1:Nt],corr_cen(1:Nt),corr_err(1:Nt),'x:');
    title('correlator <s\_spa(t)*s\_spa(0)>')

   conn_cen=mean(corr_spaspa)-mean(s_spa)^2;%%% note: avg along 1st dim --> row
   conn_err=std (corr_spaspa)/sqrt(n_meas); %%% note: std along 1st dim --> row
   conn_cen(find(conn_cen<0))=nan;          %%% note: overwrite neg.entries nan
   conn_err( isnan(conn_cen))=nan;          %%% note: overwrite correspond. nan
   conn_cen(Nt/2:Nt) % conn_cen([1:5,Nt-5:Nt])
   conn_err(Nt/2:Nt) % conn_err([1:5,Nt-5:Nt])
   figure(6); errorbar([1:Nt],conn_cen(1:Nt),conn_err(1:Nt),'x:');
    title('connected <s\_spa(t)*s\_spa(0)>-<s\_spa>^2')
%   figure(7); semilogy(nan,nan); hold on; errorbar([1:Nt],corr_cen(1:Nt),...
%    corr_err(1:Nt).*(corr_cen(1:Nt)>corr_err(1:Nt)),corr_err(1:Nt),'+:')
%   figure(8); plot([1:Nt-1],log(corr_cen([Nt,1:Nt-2])./corr_cen(2:Nt))/2,'.');
%   axis([0,Nt,-4*sqrt(msq),4*sqrt(msq)])

return

   system('rm -f durr_u1compact.ps; umask 022')
   for n=1:6
      figure(n); refresh; pause(1)
      eval(sprintf('print -dpsc -append durr_u1compact.ps'))
   end; clear n

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,acc,hitsize]=u1compact_metropolis(U,be,hitsize,n_hits)
   [~,Nx,Ny,Nz,Nt]=size(U);
   table_m=randi( 4,[1,4*Nx*Ny*Nz*Nt]);
   table_x=randi(Nx,[1,4*Nx*Ny*Nz*Nt]);
   table_y=randi(Ny,[1,4*Nx*Ny*Nz*Nt]);
   table_z=randi(Nz,[1,4*Nx*Ny*Nz*Nt]);
   table_t=randi(Nt,[1,4*Nx*Ny*Nz*Nt]);
   acc=0;
   for touch=1:4*Nx*Ny*Nz*Nt %%% in total 1 MH sweep with randomized access
      m=table_m(touch); %%% direction
      x=table_x(touch); xplu=mod(x-0,Nx)+1; xmin=mod(x-2,Nx)+1;
      y=table_y(touch); yplu=mod(y-0,Ny)+1; ymin=mod(y-2,Ny)+1;
      z=table_z(touch); zplu=mod(z-0,Nz)+1; zmin=mod(z-2,Nz)+1;
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      switch m
       case 1
         staple=U(2,x,y   ,z,t) *U(1,x,yplu,z,t)*U(2,xplu,y   ,z,t)' ...
               +U(2,x,ymin,z,t)'*U(1,x,ymin,z,t)*U(2,xplu,ymin,z,t)  ...
               +U(3,x,y,z   ,t) *U(1,x,y,zplu,t)*U(3,xplu,y,z   ,t)' ...
               +U(3,x,y,zmin,t)'*U(1,x,y,zmin,t)*U(3,xplu,y,zmin,t)  ...
               +U(4,x,y,z,t   ) *U(1,x,y,z,tplu)*U(4,xplu,y,z,t   )' ...
               +U(4,x,y,z,tmin)'*U(1,x,y,z,tmin)*U(4,xplu,y,z,tmin);
       case 2
         staple=U(1,x   ,y,z,t) *U(2,xplu,y,z,t)*U(1,x   ,yplu,z,t)' ...
               +U(1,xmin,y,z,t)'*U(2,xmin,y,z,t)*U(1,xmin,yplu,z,t)  ...
               +U(3,x,y,   z,t) *U(2,x,y,zplu,t)*U(3,x,yplu,z   ,t)' ...
               +U(3,x,y,zmin,t)'*U(2,x,y,zmin,t)*U(3,x,yplu,zmin,t)  ...
               +U(4,x,y,z,   t) *U(2,x,y,z,tplu)*U(4,x,yplu,z,   t)' ...
               +U(4,x,y,z,tmin)'*U(2,x,y,z,tmin)*U(4,x,yplu,z,tmin);
       case 3
         staple=U(1,   x,y,z,t) *U(3,xplu,y,z,t)*U(1,x   ,y,zplu,t)' ...
               +U(1,xmin,y,z,t)'*U(3,xmin,y,z,t)*U(1,xmin,y,zplu,t)  ...
               +U(2,x,   y,z,t) *U(3,x,yplu,z,t)*U(2,x,y   ,zplu,t)' ...
               +U(2,x,ymin,z,t)'*U(3,x,ymin,z,t)*U(2,x,ymin,zplu,t)  ...
               +U(4,x,y,z,   t) *U(3,x,y,z,tplu)*U(4,x,y,zplu,t   )' ...
               +U(4,x,y,z,tmin)'*U(3,x,y,z,tmin)*U(4,x,y,zplu,tmin);
       case 4
         staple=U(1,x   ,y,z,t) *U(4,xplu,y,z,t)*U(1,x   ,y,z,tplu)' ...
               +U(1,xmin,y,z,t)'*U(4,xmin,y,z,t)*U(1,xmin,y,z,tplu)  ...
               +U(2,x,y   ,z,t) *U(4,x,yplu,z,t)*U(2,x,y   ,z,tplu)' ...
               +U(2,x,ymin,z,t)'*U(4,x,ymin,z,t)*U(2,x,ymin,z,tplu)  ...
               +U(3,x,y,z   ,t) *U(4,x,y,zplu,t)*U(3,x,y,z   ,tplu)' ...
               +U(3,x,y,zmin,t)'*U(4,x,y,zmin,t)*U(3,x,y,zmin,tplu);
      end
      r=rand(2,n_hits); %%% note: sufficiently large for table at each touch
      link_old=U(m,x,y,z,t);
      s_old=1-real(link_old*staple');
      for loop_hits=1:n_hits
         link_new=exp(complex(0,hitsize*(r(1,loop_hits)-0.5)))*link_old;
         s_new=1-real(link_new*staple');
         if (r(2,loop_hits)<exp(-be*(s_new-s_old)))
            link_old=link_new;
            s_old=s_new;
            acc=acc+1;
         end
      end; clear loop_hits
      U(m,x,y,z,t)=link_old;
   end; clear touch
   acc=acc/(4*Nx*Ny*Nz*Nt*n_hits); %%% note: turn into acceptancerate
   hitsize=acc/0.8*hitsize; %%% note: adjustment for next call
end

function [U,acc]=u1compact_overrelax(U,be)
   [~,Nx,Ny,Nz,Nt]=size(U);
   table_m=randi( 4,[1,4*Nx*Ny*Nz*Nt]);
   table_x=randi(Nx,[1,4*Nx*Ny*Nz*Nt]);
   table_y=randi(Ny,[1,4*Nx*Ny*Nz*Nt]);
   table_z=randi(Nz,[1,4*Nx*Ny*Nz*Nt]);
   table_t=randi(Nt,[1,4*Nx*Ny*Nz*Nt]);
   table_r=rand([1,4*Nx*Ny*Nz*Nt]);
   acc=0;
   for touch=1:4*Nx*Ny*Nz*Nt %%% on average 1 OR update with randomized access
      m=table_m(touch); %%% direction
      x=table_x(touch); xplu=mod(x-0,Nx)+1; xmin=mod(x-2,Nx)+1;
      y=table_y(touch); yplu=mod(y-0,Ny)+1; ymin=mod(y-2,Ny)+1;
      z=table_z(touch); zplu=mod(z-0,Nz)+1; zmin=mod(z-2,Nz)+1;
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      switch m
       case 1
         staple=U(2,x,y   ,z,t) *U(1,x,yplu,z,t)*U(2,xplu,y   ,z,t)' ...
               +U(2,x,ymin,z,t)'*U(1,x,ymin,z,t)*U(2,xplu,ymin,z,t)  ...
               +U(3,x,y,z   ,t) *U(1,x,y,zplu,t)*U(3,xplu,y,z   ,t)' ...
               +U(3,x,y,zmin,t)'*U(1,x,y,zmin,t)*U(3,xplu,y,zmin,t)  ...
               +U(4,x,y,z,t   ) *U(1,x,y,z,tplu)*U(4,xplu,y,z,t   )' ...
               +U(4,x,y,z,tmin)'*U(1,x,y,z,tmin)*U(4,xplu,y,z,tmin);
       case 2
         staple=U(1,x   ,y,z,t) *U(2,xplu,y,z,t)*U(1,x   ,yplu,z,t)' ...
               +U(1,xmin,y,z,t)'*U(2,xmin,y,z,t)*U(1,xmin,yplu,z,t)  ...
               +U(3,x,y,   z,t) *U(2,x,y,zplu,t)*U(3,x,yplu,z   ,t)' ...
               +U(3,x,y,zmin,t)'*U(2,x,y,zmin,t)*U(3,x,yplu,zmin,t)  ...
               +U(4,x,y,z,   t) *U(2,x,y,z,tplu)*U(4,x,yplu,z,   t)' ...
               +U(4,x,y,z,tmin)'*U(2,x,y,z,tmin)*U(4,x,yplu,z,tmin);
       case 3
         staple=U(1,   x,y,z,t) *U(3,xplu,y,z,t)*U(1,x   ,y,zplu,t)' ...
               +U(1,xmin,y,z,t)'*U(3,xmin,y,z,t)*U(1,xmin,y,zplu,t)  ...
               +U(2,x,   y,z,t) *U(3,x,yplu,z,t)*U(2,x,y   ,zplu,t)' ...
               +U(2,x,ymin,z,t)'*U(3,x,ymin,z,t)*U(2,x,ymin,zplu,t)  ...
               +U(4,x,y,z,   t) *U(3,x,y,z,tplu)*U(4,x,y,zplu,t   )' ...
               +U(4,x,y,z,tmin)'*U(3,x,y,z,tmin)*U(4,x,y,zplu,tmin);
       case 4
         staple=U(1,x   ,y,z,t) *U(4,xplu,y,z,t)*U(1,x   ,y,z,tplu)' ...
               +U(1,xmin,y,z,t)'*U(4,xmin,y,z,t)*U(1,xmin,y,z,tplu)  ...
               +U(2,x,y   ,z,t) *U(4,x,yplu,z,t)*U(2,x,y   ,z,tplu)' ...
               +U(2,x,ymin,z,t)'*U(4,x,ymin,z,t)*U(2,x,ymin,z,tplu)  ...
               +U(3,x,y,z   ,t) *U(4,x,y,zplu,t)*U(3,x,y,z   ,tplu)' ...
               +U(3,x,y,zmin,t)'*U(4,x,y,zmin,t)*U(3,x,y,zmin,tplu);
      end

      if (0)
         link_old=U(m,x,y,z,t);
         s_old=1-real(link_old*staple');
         staple_proj=staple/abs(staple);
         link_new=staple_proj*link_old'*staple_proj;
         s_new=1-real(link_new*staple');
         if (table_r(touch)<exp(-be*(s_new-s_old)))
            U(m,x,y,z,t)=link_new;
            acc=acc+1;
         end
      else
         staple_proj=staple/abs(staple);
         U(m,x,y,z,t)=staple_proj*U(m,x,y,z,t)'*staple_proj;
         acc=acc+1;
      end

   end; clear touch
   acc=acc/(4*Nx*Ny*Nz*Nt);
end

function [s_avg,s_spa,s_tem,s_loc]=action_wils(U)
   [~,Nx,Ny,Nz,Nt]=size(U);
   s_avg=nan; s_spa(Nt)=0; s_tem(Nt)=0; s_loc=nan(Nx,Ny,Nz,Nt); obj(1:6)=nan;
   for t=1:Nt; tplu=mod(t-0,Nt)+1;
   for z=1:Nz; zplu=mod(z-0,Nz)+1;
   for y=1:Ny; yplu=mod(y-0,Ny)+1;
   for x=1:Nx; xplu=mod(x-0,Nx)+1;
      obj(1)=real(U(1,x,y,z,t)*U(2,xplu,y,z,t)*U(1,x,yplu,z,t)'*U(2,x,y,z,t)');
      obj(2)=real(U(1,x,y,z,t)*U(3,xplu,y,z,t)*U(1,x,y,zplu,t)'*U(3,x,y,z,t)');
      obj(3)=real(U(1,x,y,z,t)*U(4,xplu,y,z,t)*U(1,x,y,z,tplu)'*U(4,x,y,z,t)');
      obj(4)=real(U(2,x,y,z,t)*U(3,x,yplu,z,t)*U(2,x,y,zplu,t)'*U(3,x,y,z,t)');
      obj(5)=real(U(2,x,y,z,t)*U(4,x,yplu,z,t)*U(2,x,y,z,tplu)'*U(4,x,y,z,t)');
      obj(6)=real(U(3,x,y,z,t)*U(4,x,y,zplu,t)*U(3,x,y,z,tplu)'*U(4,x,y,z,t)');
      spa=mean(ones(1,3)-obj([1,2,4])); %%% note: average spatial  plaquette
      tem=mean(ones(1,3)-obj([3,5,6])); %%% note: average temporal plaquette
      s_spa(t)=s_spa(t)+spa;
      s_tem(t)=s_tem(t)+tem;
      s_loc(x,y,z,t)=(spa+tem)/2; %%% note: average plaquette in (x,y,z,t)
   end
   end
   end
   end
   s_avg=mean(s_loc(:));   %%% note: average plaquette in present lattice
   s_spa=s_spa/(Nx*Ny*Nz); %%% note: average spa-plaquette per time-slice
   s_tem=s_tem/(Nx*Ny*Nz); %%% note: average tem-plaquette per time-slice
end

function [corr_spaspa,corr_spatem,corr_temspa,corr_temtem]=calc_corr(s_spa,s_tem)
   Nt=size([s_spa;s_tem],2);
   corr_spaspa=nan(1,Nt);
   corr_spatem=nan(1,Nt);
   corr_temspa=nan(1,Nt);
   corr_temtem=nan(1,Nt);
   shift_spa=s_spa;
   shift_tem=s_tem;
   for delta=1:Nt
      shift_spa=circshift(shift_spa,-1); %%% note: one to the left
      shift_tem=circshift(shift_tem,-1); %%% note: one to the left
      corr_spaspa(delta)=mean(s_spa.*shift_spa);
      corr_spatem(delta)=mean(s_spa.*shift_tem);
      corr_temspa(delta)=mean(s_tem.*shift_spa);
      corr_temtem(delta)=mean(s_tem.*shift_tem);
   end
end

