function durr_phiwils(Nx,Ny,Nz,Nt,msq,lam,n_ther,n_meas,n_sepa,n_over,phi)
% usage: durr_phiwils(Nx,Ny,Nz,Nt,msq,lam,n_ther,n_meas,n_sepa,n_over,phi);

   %%% solu:=solve({ka=(1-2*la)/(4+msq/2),la=lam/4*ka^2},{msq,lam});
   %%% eval(solu,{la=0.275376,ka=0.288998}); %%% DebashishBanerjee
   %%% {lam = 13.18853038, msq = -4.890995786}

   if nargin< 1, Nx=      24-8, end
   if nargin< 2, Ny=      24-8, end
   if nargin< 3, Nz=      24-8, end
   if nargin< 4, Nt=      48-16, end
   if nargin< 5, msq=    0.1, end
   if nargin< 6, lam=    0.2, end
   if nargin< 7, n_ther= 100/10, end %%% number of configs to be skipped
   if nargin< 8, n_meas=1000/100, end %%% number of configs for measurements
   if nargin< 9, n_sepa=  10, end %%% number of packages between configs
   if nargin<10, n_over=   4, end %%% OR-to-MH sweep ratio in one package
   if nargin<11, phi=ones(Nx,Ny,Nz,Nt);
      elseif Nx~=size(phi,1), error('durr_phiwils: Nx~=size(phi,1)'),
      elseif Ny~=size(phi,2), error('durr_phiwils: Ny~=size(phi,2)'),
      elseif Nz~=size(phi,3), error('durr_phiwils: Nz~=size(phi,3)'),
      elseif Nt~=size(phi,4), error('durr_phiwils: Nt~=size(phi,4)'),
   end
   for n=1:16
      figure(n); clf
   end; clear n
   figure(01)
   if max(abs(phi(:))-mean(abs(phi(:))))>0
      contourf(squeeze(max(max(abs(phi))))); axis ij; axis equal;
   else
      contourf([1:Nz]'*[1:Nt],[-1.0,-1.0]);  axis ij; axis equal;
   end
   title('max(max(|phi|)) after initialization'); xlabel('t/a'); ylabel('z/a')

   hitsize=1.0;
   for loop_ther=1:n_ther
      for loop_sepa=1:n_sepa
         [phi,acc_mh,hitsize]=phiwils_metropolis(phi,msq,lam,hitsize);
      end; clear loop_sepa
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f\n',...
       loop_ther-n_ther,hitsize,acc_mh)
   end; clear loop_ther acc_mh acc_or
   figure(02); contourf(squeeze(max(max(abs(phi))))); axis ij; axis equal;
    title('max(max(|phi|)) afterthermalization'); xlabel('t/a'); ylabel('z/a')

   phi0=nan(n_meas,1); phi1=phi0; phi2=phi0; phi3=phi0; phi4=phi0;
   skin=nan(n_meas,1); spot=nan(n_meas,1); corr=nan(n_meas,Nt);
   for loop_meas=1:n_meas
      %%% new configuration means n_sepa packages of n_over*OR+1*MH sweeps
      for loop_sepa=1:n_sepa; acc_or=nan(1,n_over);
         for loop_over=1:n_over
            [phi,acc_or(loop_over)]=phiwils_overrelax(phi,msq,lam);
         end; clear loop_over
         [phi,acc_mh,hitsize]=phiwils_metropolis(phi,msq,lam,hitsize);
      end; clear loop_sepa
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f, acc_or=%6.4f\n',...
       loop_meas,hitsize,acc_mh,mean(acc_or))
      %%% start measurements on configuration number loop_meas
      phi0(loop_meas)=mean(abs(phi(:)));
      phi1(loop_meas)=mean(phi(:).^1);
      phi2(loop_meas)=mean(phi(:).^2);
      phi3(loop_meas)=mean(phi(:).^3);
      phi4(loop_meas)=mean(phi(:).^4);
      [~,skin(loop_meas),spot(loop_meas)]=action_wils(phi,msq,lam);
      corr(loop_meas,:)=calc_corr(phi);
   end; clear loop_meas acc_mh acc_or
   skin=skin/(Nx*Ny*Nz*Nt); %%% per site
   spot=spot/(Nx*Ny*Nz*Nt); %%% per site
   figure(03); contourf(squeeze(max(max(abs(phi))))); axis ij; axis equal;
    title('max(max(abs(phi))) after finishing'); xlabel('t/a'); ylabel('z/a')

   fprintf('|phi| : %9.6f +/- %8.6f\n',mean(phi0),std(phi0)/sqrt(n_meas))
   fprintf('phi^1 : %9.6f +/- %8.6f\n',mean(phi1),std(phi1)/sqrt(n_meas))
   fprintf('phi^2 : %9.6f +/- %8.6f\n',mean(phi2),std(phi2)/sqrt(n_meas))
   fprintf('phi^3 : %9.6f +/- %8.6f\n',mean(phi3),std(phi3)/sqrt(n_meas))
   fprintf('phi^4 : %9.6f +/- %8.6f\n',mean(phi4),std(phi4)/sqrt(n_meas))
   fprintf('susc  : %9.6f +/- %8.6f\n',mean(phi2)-mean(phi0)^2,nan)
   fprintf('skin  : %9.6f +/- %8.6f\n',mean(skin),std(skin)/sqrt(n_meas))
   fprintf('spot  : %9.6f +/- %8.6f\n',mean(spot),std(spot)/sqrt(n_meas))

   figure(04); plot(phi0); title('history of << |phi| >>')
   figure(05); plot(phi1); title('history of << phi^1 >>')
   figure(06); plot(phi2); title('history of << phi^2 >>')
   figure(07); plot(phi3); title('history of << phi^3 >>')
   figure(08); plot(phi4); title('history of << phi^4 >>')
   figure(09); plot(skin); title('history of << skin >>')
   figure(10); plot(spot); title('history of << spot >>')

   corr_cen=mean(corr);              %%% note: avg along 1st dim --> row-vector
   corr_err=std (corr)/sqrt(n_meas); %%% note: std along 1st dim --> row-vector
   corr_cen(find(corr_cen<0))=nan;   %%% note: overwrite neg. entries with nan
   corr_err( isnan(corr_cen))=nan;   %%% note: overwrite corre. error with nan
   corr_cen(1:Nt/4) % corr_cen([1:5,Nt-5:Nt])
   corr_err(1:Nt/4) % corr_err([1:5,Nt-5:Nt])
   figure(11); errorbar([1:Nt],corr_cen(1:Nt),corr_err(1:Nt),'x:');
    title('correlator <phi(t)*phi(0)>')
   figure(12); semilogy(nan,nan); hold on; errorbar([1:Nt],corr_cen(1:Nt),...
    corr_err(1:Nt).*(corr_cen(1:Nt)>corr_err(1:Nt)),corr_err(1:Nt),'+:')
   figure(13); plot([1:Nt-1],log(corr_cen([Nt,1:Nt-2])./corr_cen(2:Nt))/2,'.');
   axis([0,Nt,-4*sqrt(msq),4*sqrt(msq)])

   conn_cen=mean(corr)-mean(phi1)^2; %%% note: avg along 1st dim --> row-vector
   conn_err=std (corr)/sqrt(n_meas); %%% note: std along 1st dim --> row-vector
   conn_cen(find(conn_cen<0))=nan;   %%% note: overwrite neg. entries with nan
   conn_err( isnan(conn_cen))=nan;   %%% note: overwrite corre. error with nan
   conn_cen(1:Nt/4) % conn_cen([1:5,Nt-5:Nt])
   conn_err(1:Nt/4) % conn_err([1:5,Nt-5:Nt])
   figure(14); errorbar([1:Nt],conn_cen(1:Nt),conn_err(1:Nt),'x:');
    title('connected <phi(t)*phi(0)>-<phi>^2')
   figure(15); semilogy(nan,nan); hold on; errorbar([1:Nt],conn_cen(1:Nt),...
    conn_err(1:Nt).*(conn_cen(1:Nt)>conn_err(1:Nt)),conn_err(1:Nt),'+:')
   figure(16); plot([1:Nt-1],log(conn_cen([Nt,1:Nt-2])./conn_cen(2:Nt))/2,'.');
    axis([0,Nt,-4*sqrt(msq),4*sqrt(msq)])

return

   system('rm -f durr_phiwils.ps; umask 022')
   for n=1:16
      figure(n); refresh; pause(1)
      eval(sprintf('print -dpsc -append durr_phiwils.ps'))
   end; clear n

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stot,skin,spot]=action_wils(phi,msq,lam)
   tmp=phi.*laplace_wils(phi);
   skin=-0.5*sum(tmp(:));
   spot=0.5*msq*sum(phi(:).^2)+0.25*lam*sum(phi(:).^4);
   stot=skin+spot;
end

function [lapphi]=laplace_wils(phi)
   [Nx,Ny,Nz,Nt]=size(phi);
   lapphi=nan(Nx,Ny,Nz,Nt);
   for t=1:Nt; tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
   for z=1:Nz; zplu=mod(z-0,Nz)+1; zmin=mod(z-2,Nz)+1;
   for y=1:Ny; yplu=mod(y-0,Ny)+1; ymin=mod(y-2,Ny)+1;
   for x=1:Nx; xplu=mod(x-0,Nx)+1; xmin=mod(x-2,Nx)+1;
      lapphi(x,y,z,t)=-8*phi(x,y,z,t)+phi(xmin,y,z,t)+phi(xplu,y,z,t) ...
                                     +phi(x,ymin,z,t)+phi(x,yplu,z,t) ...
                                     +phi(x,y,zmin,t)+phi(x,y,zplu,t) ...
                                     +phi(x,y,z,tmin)+phi(x,y,z,tplu);
   end
   end
   end
   end
end

function [phi,acc,hitsize]=phiwils_metropolis(phi,msq,lam,hitsize)
   [Nx,Ny,Nz,Nt]=size(phi);
   table_x=randi(Nx,[1,Nx*Ny*Nz*Nt]);
   table_y=randi(Ny,[1,Nx*Ny*Nz*Nt]);
   table_z=randi(Nz,[1,Nx*Ny*Nz*Nt]);
   table_t=randi(Nt,[1,Nx*Ny*Nz*Nt]);
   table_r=rand([2,Nx*Ny*Nz*Nt]);
   acc=0;
   for touch=1:Nx*Ny*Nz*Nt %%% in total 1 MH sweep with randomized access
      x=table_x(touch); xplu=mod(x-0,Nx)+1; xmin=mod(x-2,Nx)+1;
      y=table_y(touch); yplu=mod(y-0,Ny)+1; ymin=mod(y-2,Ny)+1;
      z=table_z(touch); zplu=mod(z-0,Nz)+1; zmin=mod(z-2,Nz)+1;
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      %%% symmetric proposal for new phi(x,y,z,t)
      phiold=phi(x,y,z,t);
      phinew=phiold+hitsize*(table_r(1,touch)-0.5);

      switch(3)
       case(1)
         phi(x,y,z,t)=phinew;
         snew=action_wils(phi,msq,lam);
         phi(x,y,z,t)=phiold;
         sold=action_wils(phi,msq,lam);
       case(2)
         sold=(4+0.5*msq)*phiold^2+0.25*lam*phiold^4 ...
             +phiold*(-0.5*phi(xmin,y,z,t)-0.5*phi(xplu,y,z,t) ...
                      -0.5*phi(x,ymin,z,t)-0.5*phi(x,yplu,z,t) ...
                      -0.5*phi(x,y,zmin,t)-0.5*phi(x,y,zplu,t) ...
                      -0.5*phi(x,y,z,tmin)-0.5*phi(x,y,z,tplu))...
             +(-0.5*phi(xplu,y,z,t)-0.5*phi(xmin,y,z,t) ...
               -0.5*phi(x,yplu,z,t)-0.5*phi(x,ymin,z,t) ...
               -0.5*phi(x,y,zplu,t)-0.5*phi(x,y,zmin,t) ...
               -0.5*phi(x,y,z,tplu)-0.5*phi(x,y,z,tmin))*phiold;
         snew=(4+0.5*msq)*phinew^2+0.25*lam*phinew^4 ...
             +phinew*(-0.5*phi(xmin,y,z,t)-0.5*phi(xplu,y,z,t) ...
                      -0.5*phi(x,ymin,z,t)-0.5*phi(x,yplu,z,t) ...
                      -0.5*phi(x,y,zmin,t)-0.5*phi(x,y,zplu,t) ...
                      -0.5*phi(x,y,z,tmin)-0.5*phi(x,y,z,tplu))...
             +(-0.5*phi(xplu,y,z,t)-0.5*phi(xmin,y,z,t) ...
               -0.5*phi(x,yplu,z,t)-0.5*phi(x,ymin,z,t) ...
               -0.5*phi(x,y,zplu,t)-0.5*phi(x,y,zmin,t) ...
               -0.5*phi(x,y,z,tplu)-0.5*phi(x,y,z,tmin))*phinew;
       case(3)
         neighborsum=phi(xmin,y,z,t)+phi(xplu,y,z,t) ...
                    +phi(x,ymin,z,t)+phi(x,yplu,z,t) ...
                    +phi(x,y,zmin,t)+phi(x,y,zplu,t) ...
                    +phi(x,y,z,tmin)+phi(x,y,z,tplu);
         tmp=0.5*phiold^2; sold=tmp*(8+msq+lam*tmp)-phiold*neighborsum;
         tmp=0.5*phinew^2; snew=tmp*(8+msq+lam*tmp)-phinew*neighborsum;
      end

      if (table_r(2,touch)<exp(-snew+sold))
         phi(x,y,z,t)=phinew;
         acc=acc+1;
      end
   end; clear touch
   acc=acc/(Nx*Ny*Nz*Nt);   %%% note: turn into acceptancerate
   hitsize=acc/0.8*hitsize; %%% note: adjustment for next call
end

function [phi,acc]=phiwils_overrelax(phi,msq,lam)
   [Nx,Ny,Nz,Nt]=size(phi);
   table_x=randi(Nx,[1,Nx*Ny*Nz*Nt]);
   table_y=randi(Ny,[1,Nx*Ny*Nz*Nt]);
   table_z=randi(Nz,[1,Nx*Ny*Nz*Nt]);
   table_t=randi(Nt,[1,Nx*Ny*Nz*Nt]);
   table_r=rand([1,Nx*Ny*Nz*Nt]);
   acc=0;
   for touch=1:Nx*Ny*Nz*Nt %%% on average 1 OR update with randomized access
      x=table_x(touch); xplu=mod(x-0,Nx)+1; xmin=mod(x-2,Nx)+1;
      y=table_y(touch); yplu=mod(y-0,Ny)+1; ymin=mod(y-2,Ny)+1;
      z=table_z(touch); zplu=mod(z-0,Nz)+1; zmin=mod(z-2,Nz)+1;
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      phiold=phi(x,y,z,t);
      phinew=0.25*(phi(xplu,y,z,t)+phi(xmin,y,z,t) ...
                  +phi(x,yplu,z,t)+phi(x,ymin,z,t) ...
                  +phi(x,y,zplu,t)+phi(x,y,zmin,t) ...
                  +phi(x,y,z,tplu)+phi(x,y,z,tmin))...
            -phiold;
      %%% option one: complete s_old and s_new as in case(3) of MH-routine
      %neighborsum=phi(xmin,y,z,t)+phi(xplu,y,z,t) ...
      %           +phi(x,ymin,z,t)+phi(x,yplu,z,t) ...
      %           +phi(x,y,zmin,t)+phi(x,y,zplu,t) ...
      %           +phi(x,y,z,tmin)+phi(x,y,z,tplu);
      %sold=(4+0.5*msq)*phiold^2+0.25*lam*phiold^4-phiold*neighborsum;
      %snew=(4+0.5*msq)*phinew^2+0.25*lam*phinew^4-phinew*neighborsum;
      %%% option two: omit skin since this part is invariant
      tmp=0.5*phiold^2; sold=tmp*(msq+lam*tmp); %%% omit skin
      tmp=0.5*phinew^2; snew=tmp*(msq+lam*tmp); %%% omit skin
      %%% end options
      if (table_r(touch)<exp(-snew+sold))
         phi(t)=phinew;
         acc=acc+1;
      end
   end; clear touch
   acc=acc/(Nx*Ny*Nz*Nt);
end

function [corr]=calc_corr(phi)
   [Nx,Ny,Nz,Nt]=size(phi);
   corr=nan(1,Nt);
   avg=mean(reshape(phi,Nx*Ny*Nz,Nt)); %%% note: time-slice average
   avg_shift=avg;
   for delta=1:Nt
      avg_shift=circshift(avg_shift,-1); %%% note: shift one to the left
      tmp=avg.*avg_shift;
      corr(delta)=mean(tmp);
   end
end

