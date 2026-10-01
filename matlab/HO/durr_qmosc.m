function [x]=durr_qmosc(Nt,mass,omega,n_ther,n_meas,n_sepa,n_over,x)
% usage: [~]=durr_qmosc(Nt,mass,omega,n_ther,n_meas,n_sepa,n_over,x);

   if nargin<1, Nt=    1200, end
   if nargin<2, mass=120/Nt, end
   if nargin<3, omega=mass , end
   if nargin<4, n_ther= 100, end %%% number of configs to be skipped
   if nargin<5, n_meas=1000, end %%% number of configs for measurements
   if nargin<6, n_sepa= 100, end %%% number of packages between configs
   if nargin<7, n_over=   4, end %%% OR-to-MH sweep ratio in one package
   if nargin<8, x=zeros(1,Nt); elseif Nt~=length(x), error('durr_qmosc: Nt~=length(x)'), end
   figure(01); plot([1:Nt],x(:),'k.:'); title('path after initialization')

   hitsize=1;
   for loop_ther=1:n_ther
      for loop_sepa=1:n_sepa
        %if (n_over>0)
        %   [x,acc_or]=qmosc_overrelax(x,mass,omega);
        %end
         [x,acc_mh,hitsize]=qmosc_metropolis(x,mass,omega,hitsize);
      end; clear loop_sepa
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f\n',loop_ther-n_ther,hitsize,acc_mh)
   end; clear loop_ther acc_mh
   figure(02); plot([1:Nt],x(:),'k.:'); title('path after thermalization')

   x1=nan(n_meas,1); x2=nan(n_meas,1); x3=nan(n_meas,1); x4=nan(n_meas,1);
   skin=nan(n_meas,1); skfh=nan(n_meas,1); spot=nan(n_meas,1);
   corr=nan(n_meas,Nt); xxyy=nan(n_meas,Nt);
   for loop_meas=1:n_meas
      %%% generate new configuration by n_sepa packages of n_over OR + 1 MH sweeps
      for loop_sepa=1:n_sepa; acc_or=nan(1,n_over);
         for loop_over=1:n_over
            [x,acc_or(loop_over)]=qmosc_overrelax(x,mass,omega);
         end
         [x,acc_mh,hitsize]=qmosc_metropolis(x,mass,omega,hitsize);
      end; clear loop_sepa loop_over
      fprintf('%6i: hitsize=%4.2f, acc_mh=%6.4f, acc_or=%6.4f\n',loop_meas,hitsize,acc_mh,mean(acc_or))
      %%% start measurements on configuration number loop_meas
      x1(loop_meas)=mean(x.^1);
      x2(loop_meas)=mean(x.^2);
      x3(loop_meas)=mean(x.^3);
      x4(loop_meas)=mean(x.^4);
      [skin(loop_meas),skfh(loop_meas),spot(loop_meas)]=calc_sss(x,mass,omega);
      corr(loop_meas,:)=calc_corr(x);
      xxyy(loop_meas,:)=calc_xxyy(x);
   end; clear loop_meas acc acc_or
   figure(03); plot([1:Nt],x(:),'k.:'); title('path after last measurement')

   fprintf('x^1 : %7.4g +/- %6.4g\n',mean(x1),std(x1)/sqrt(n_meas))
   fprintf('x^2 : %7.4g +/- %6.4g\n',mean(x2),std(x2)/sqrt(n_meas))
   fprintf('x^3 : %7.4g +/- %6.4g\n',mean(x3),std(x3)/sqrt(n_meas))
   fprintf('x^4 : %7.4g +/- %6.4g\n',mean(x4),std(x4)/sqrt(n_meas))
   fprintf('skin: %7.4g +/- %6.4g\n',mean(skin),std(skin)/sqrt(n_meas))
   fprintf('skfh: %7.4g +/- %6.4g\n',mean(skfh),std(skfh)/sqrt(n_meas))
   fprintf('spot: %7.4g +/- %6.4g\n',mean(spot),std(spot)/sqrt(n_meas))

   figure(04); plot(x1); title('history of x^1')
   figure(05); plot(x2); title('history of x^2')
   figure(06); plot(x3); title('history of x^3')
   figure(07); plot(x4); title('history of x^4')
   figure(08); plot(skin); title('history of skin')
   figure(09); plot(skfh); title('history of skfh')
   figure(10); plot(spot); title('history of spot')

   corr_cen=mean(corr);              %%% note: avg along 1st direction --> row-vector
   corr_err=std (corr)/sqrt(n_meas); %%% note: std along 1st direction --> row-vector
   corr_cen(find(corr_cen<0))=nan;   %%% note: overwrite negative entries with nan
   corr_err( isnan(corr_cen))=nan;   %%% note: overwrite respective error with nan
  %find    (~isnan(corr_cen))
  %corr_cen(~isnan(corr_cen))
   corr_cen([1:5,Nt-5:Nt])
   corr_err([1:5,Nt-5:Nt])
   figure(11); errorbar([1:100],corr_cen(1:100),corr_err(1:100),'x:'); title('correlator x(t)^1*x(0)^1 - <x^1>^2')
   figure(12); semilogy(nan,nan); hold on; errorbar([1:100],corr_cen(1:100),...
               corr_err(1:100).*(corr_cen(1:100)>corr_err(1:100)),corr_err(1:100),'+:'); hold off
   figure(13); plot([1:99],0.5*log(corr_cen([Nt,1:98])./corr_cen(2:100)),'.'); axis([0,100,0,2*omega])

   xxyy_cen=mean(xxyy)-mean(x2)^2;   %%% note: avg along 1st direction --> row-vector
   xxyy_err=std (xxyy)/sqrt(n_meas); %%% note: std along 1st direction --> row-vector
   xxyy_cen(find(xxyy_cen<0))=nan;   %%% note: overwrite negative entries with nan
   xxyy_err( isnan(xxyy_cen))=nan;   %%% note: overwrite respective error with nan
  %find    (~isnan(xxyy_cen))
  %xxyy_cen(~isnan(xxyy_cen))
   xxyy_cen([1:5,Nt-5:Nt])
   xxyy_err([1:5,Nt-5:Nt])
   figure(14); errorbar([1:100],xxyy_cen(1:100),xxyy_err(1:100),'x:'); title('correlator x(t)^2*x(0)^2 - <x^2>^2')
   figure(15); semilogy(nan,nan); hold on; errorbar([1:100],xxyy_cen(1:100),...
               xxyy_err(1:100).*(xxyy_cen(1:100)>xxyy_err(1:100)),xxyy_err(1:100),'+:'); hold off
   figure(16); plot([1:99],0.5*log(xxyy_cen([Nt,1:98])./xxyy_cen(2:100)),'.'); axis([0,100,0,4*omega])

return

   moni=[x1,x2,x3,x4,skin,skfh,spot];
   save(['durr_qmosc.dat_',num2str(Nt),'_moni'],'moni','-ascii');
   save(['durr_qmosc.dat_',num2str(Nt),'_corr'],'corr','-ascii');
   save(['durr_qmosc.dat_',num2str(Nt),'_xxyy'],'xxyy','-ascii');

return

   system('rm -f durr_qmosc.ps; umask 022')
   for n=1:16
      figure(n); refresh; pause(1)
      eval(sprintf('print -dpsc -append durr_qmosc.ps'))
   end; clear n

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,acc,hitsize]=qmosc_metropolis(x,mass,omega,hitsize)
   Nt=length(x);
   table_t=randi(Nt,[1,Nt]);
   table_r=rand([2,Nt]);
   acc=0;
   for touch=1:Nt %%% in total 1 MH sweep with randomized access
      t=table_t(touch);
      tplu=mod(t-0,Nt)+1;
      tmin=mod(t-2,Nt)+1;
      xnew=x(t)+hitsize*(table_r(1,touch)-0.5);
      sold=0.5*mass*((x(tplu)-x(t))^2+(x(t)-x(tmin))^2+omega^2*x(t)^2); %%% skin+spot
      snew=0.5*mass*((x(tplu)-xnew)^2+(xnew-x(tmin))^2+omega^2*xnew^2); %%% skin+spot
      if (table_r(2,touch)<exp(-snew+sold))
         x(t)=xnew;
         acc=acc+1;
      end
   end; clear touch
   acc=acc/Nt;
   hitsize=acc/0.8*hitsize; %%% note: adjustment for next call
end

function [x,acc]=qmosc_overrelax(x,mass,omega)
   Nt=length(x);
   table_t=randi(Nt,[1,Nt]);
   table_r=rand([1,Nt]);
   acc=0;
   for touch=1:Nt %%% on average 1 approximate OR sweep with randomized access
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      xnew=x(tplu)+x(tmin)-x(t); %%% note: simplified version, only skin is invariant
     %sold=0.5*mass*((x(tplu)-x(t))^2+(x(t)-x(tmin))^2+omega^2*x(t)^2); %%% skin+spot
     %snew=0.5*mass*((x(tplu)-xnew)^2+(xnew-x(tmin))^2+omega^2*xnew^2); %%% skin+spot
      sold=0.5*mass*(0.0000000000000000000000000000000+omega^2*x(t)^2); %%% omit skin
      snew=0.5*mass*(0.0000000000000000000000000000000+omega^2*xnew^2); %%% omit skin
      if (table_r(touch)<exp(-snew+sold))
         x(t)=xnew;
         acc=acc+1;
      end
   end; clear touch
   acc=acc/Nt;
end

function [x,acc]=qmosc_mikecreutz(x,mass,omega)
   % > solve((2+omega^2)*xold^2-2*xplu*xold-2*xold*xmin=(2+omega^2)*xnew^2-2*xplu*xnew-2*xnew*xmin,xnew);
   %   { xold , -(omega^2*xold+2*xold-2*xmin-2*xplu)/(omega^2+2) }
   Nt=length(x);
   table_t=randi(Nt,[1,Nt]);
  %table_r=rand([1,Nt]);
   acc=0;
   for touch=1:Nt %%% on agerage 1 exact OR sweep with randomized access
      t=table_t(touch); tplu=mod(t-0,Nt)+1; tmin=mod(t-2,Nt)+1;
      xnew=(x(tplu)+x(tmin))/(1.0+0.5*omega^2)-x(t); %%% note: version by Mike Creutz
     %sold=0.5*mass*((x(tplu)-x(t))^2+(x(t)-x(tmin))^2+omega^2*x(t)^2); %%% skin+spot
     %snew=0.5*mass*((x(tplu)-xnew)^2+(xnew-x(tmin))^2+omega^2*xnew^2); %%% skin+spot
     %if (table_r(touch)<exp(-snew+sold))
         x(t)=xnew;
         acc=acc+1;
     %end
   end; clear touch
   acc=acc/Nt;
end

function [skin,skfh,spot]=calc_sss(x,mass,omega)
   Nt=length(x);
   skin=0;
   skfh=0;
   spot=0;
   for t=1:Nt
      tplu=mod(t-0,Nt)+1;
      tmin=mod(t-2,Nt)+1;
      skin=skin+(x(tplu)-x(t))^2;              %%% note: naive kinetic term of action
      skfh=skfh-(x(tplu)-x(t))*(x(t)-x(tmin)); %%% note: point-split Feynmann / Hibbs
      spot=spot+x(t)^2;                        %%% note: potential w/o factor omega^2
   end; clear t
   skin=0.5*mass*skin/Nt;
   skfh=0.5*mass*skfh/Nt;
   spot=0.5*mass*omega^2*spot/Nt;
end

function [corr]=calc_corr(x)
   Nt=length(x);
   corr=nan(1,Nt);
  %for delta=[1:min(Nt/4,100),                    Nt]
   for delta=[1:min(Nt/4,100),max(Nt-Nt/4,Nt-100):Nt]
      corr(delta)=mean(x.*circshift(x,delta));
   end
end

function [corr]=calc_xxyy(x)
   Nt=length(x);
   corr=nan(1,Nt);
   xx=x.^2;
  %for delta=[1:min(Nt/4,100),                    Nt]
   for delta=[1:min(Nt/4,100),max(Nt-Nt/4,Nt-100):Nt]
      corr(delta)=mean(xx.*circshift(xx,delta));
   end
end

% rng('default'); tic, durr_qmosc(1200,0.1,0.1,100,1000,120,0), toc
% x^1 : -0.2364 +/- 0.03233
% x^2 :   50.1 +/-  0.207
% x^3 : -30.99 +/-  5.434
% x^4 :   7568 +/-  65.87
% skin: 0.4748 +/- 0.000625
% skfh: 0.02277 +/- 0.000428
% spot: 0.02505 +/- 0.0001035
% Elapsed time is 20.267805 seconds.

% rng('default'); tic, durr_qmosc(1200,0.1,0.1,100,1000,120,2), toc
% x^1 : 0.03191 +/- 0.02932
% x^2 :   49.9 +/- 0.2075
% x^3 :  3.296 +/-  4.884
% x^4 :   7497 +/-  67.41
% skin:  0.475 +/- 0.0006266
% skfh: 0.02259 +/- 0.000425
% spot: 0.02495 +/- 0.0001037
% Elapsed time is 43.909264 seconds.

% rng('default'); tic, durr_qmosc(1200,0.1,0.1,100,1000,120,6), toc
% x^1 : -0.03026 +/- 0.02853
% x^2 :  50.11 +/- 0.2001
% x^3 : 0.7225 +/-  4.823
% x^4 :   7544 +/-  65.51
% skin: 0.4757 +/- 0.0006153
% skfh: 0.02233 +/- 0.0004168
% spot: 0.02506 +/- 0.0001001
% Elapsed time is 95.311556 seconds.

% simulation plan for continuum limit
%    Nt  mass=omega  n_ther  n_meas n_sepa n_over
%   120        1.0      100   10000     12      0 <== 11.7 sec
%   150        0.8      100   10000     15      1 <== 18.9 sec
%   200        0.6      100   10000     20      2 <== 32.2 sec
%   240        0.5      100   10000     24      3 <== 48.4 sec
%   400        0.3      100   10000     40      4 <==  118 sec
%   600        0.2      100   10000     60      5 <==  246 sec
%  1200        0.1      100   10000    120      6 <==  915 sec
%  1500        0.08     100   10000    150      7 <== 1558 sec
%  2000        0.06     100   10000    200      8 <== 2953 sec
%  2400        0.05     100   10000    240      9 <== 4612 sec
%  4000        0.03     100   10000    400     10 <==  233 min
%  6000        0.02     100   10000    600     11 <==  554 min
% 12000        0.01     100   10000   1200     12 <== 2327 min

% -rw-r--r-- 1 durr durr  64645 Jun 28 12:03 durr_qmosc.dat_1200_corr.gz
% -rw-r--r-- 1 durr durr   4184 Jun 28 12:03 durr_qmosc.dat_1200_moni.gz
% -rw-r--r-- 1 durr durr  56684 Jun 28 12:03 durr_qmosc.dat_1200_xxyy.gz
