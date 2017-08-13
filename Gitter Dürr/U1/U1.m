U = complex(ones(4,Nx,Ny,Nz,Nt))

function [U,acc,hitsize] = ulcompact_metropolis(U,be,hitsize)
	[~,Nx,Ny,Nz,Nt] = size(U);
	table_x/y/z/t/m
	for touch = 1:4*Nx*Ny*Nz*Nt
		...
		switch m
			case 1; staple = ...
			case 2; staple = ...
			case 3; staple = ...
			case 4; staple = ...

		end

% alle Wege nach x+mu (außer durch die Mitte)
staplex = U(2,x,y,z,t) * U(1,x,yplus,z,t) * U(2,xplus,y,z,t)' +
	U(2,x,ymin,z,t)' * U(1,x,ymin,z,t) * U(2,xplus,ymin,z,t) +
	U(3,x,y,z,t) * U(1,x,y,zplus,t) * U(3,xplus,y,z,t)' + 
	U(3,x,y,zmin,t)' * U(1,x,y,zmin,t) * U(3,xplus,y,zmin,t) +
	U(4,x,y,z,t) * U(1,x,y,z,tplus) * U(4,xplus,y,z,t)' +
	U(4,x,y,z,tmin)' + U(1,x,y,z,tmin) * U(4,xplus,y,z,tmin);

% alle Wege nach y+mu
stapley = U(1,x,y,z,t) * U(2,xplus,y,z,t) * U(1,x,yplus,z,t)' +
	U(1,xmin,y,z,t)' * U(2,xmin,y,z,t) * U(1,xmin,yplus,z,t) +
	U(3,x,y,z,t) * ...

...

r = rand(2,n_hits);
link_old = U(m,x,y,z,t);
s_old = 1 - real(link_old * staple');

% multiple hits - more efficient - converges against heat bath algorithm
for hit = 1:n_hits

	link_new = exp( complex(0,hitsize*(r(1,hit) - 0.5)) ) * link_old
	s_old = 1 - real(link_new * staple');

	if (r(2,hit) < exp(-be * (s_new - s_old)))
		link_old = link_new
		s_old = s_new
		acc = acc + 1
	end
end

U(m,x,y,z,t) = link_old

%%% Overrelaxation:
%%% Wie MH mit n_hits = 1,

staple = ...
% staple projected
staple_proj = staple / abs(staple);
link_new = staple_proj * link_old' * staple_proj
% accept rate should be 1 - no change in action


% beachte:
acc = acc / (4*Nx*Ny*Nz*Nt*nhits)
U = U_/abs(U);

www.durr.itp.unibe.ch/BUW/2017_lattice/durr_u1compact.m
