function FDTDslab

% 2D propagation in a parallel plate waveguide w/ a dielectric slab
% 
% short circuit walls at z=0, x=0, x=W;  Mur ABC termination at z=L
%
% one period sin line source excitation at user specified location
%
% dielectric slab at user specified location
%
% September 3, 2018   A. F. Peterson
 
  L = input('Give guide length ');
  A = input('Give guide width  ');
  Nz = input('Number of cells along length ');
  Nx = input('Number of cells along width ');
  delz = L/Nz;
  delx = A/Nx;
  Nsamp = round((Nx+1)/2);
  Nzm1=Nz-1;

  Sz = input('Cell index of source along length? ');
  Sx = input('Cell index of source along width?  ');

% initialize variables
  
  BigT = 1.0e-6;       % transient excitation of 1 microsecond duration
  mu = pi*4.0e-7;      % permeability
  epsilon = 8.854e-12; % permittivity
  c = 2.998e8;         % speed of light
  eta=sqrt(mu/epsilon); % eta
  Hy=zeros(Nx,Nz);      % Hy, size is Nx x Nz
  Ex=zeros(Nx,Nz);      % Ex, size is Nx x Nz
  Ez=zeros(Nx,Nz);      % Ez, size is Nx x Nz
  Jx=zeros(Nx,Nz);      % Jx, size is Nx x Nz
  Hymem(Nx)=0;          % vector of length Hx

% Ex(jj,ii) is used to denote Ex(jj,ii-1/2)
% Ez(jj,ii) is used to denote Ez(jj-1/2,ii)
% erx(jj,ii) samples of permittivity at Ex-field locations
% erz(jj,ii) samples of permittivity at Ez-field locations

  epsx=ones(Nx,Nz);
  epsz=ones(Nx,Nz);
  zstart = input('starting index in z for slab? ');
    zend = input('  ending index in z for slab? ');
   epslb = input('          epsilon-r for slab? ');

  for ii=zstart:zend % starting and ending index of slab in z
    if(ii == zstart) % values of epsx get average, epsz in slab
      for jj=1:Nx    % all values of x
        epsx(jj,ii)=(epslb+1)/2; % average epsilon at x field locations
        epsz(jj,ii)=epslb;       % epsz doesn't get averaged
      end      
    elseif(ii == zend) % values of epsx get average, epsz out of slab
      for jj=1:Nx
        epsx(jj,ii)=(epslb+1)/2; % average epsilon at x, epsilon at z is out of the slab
      end   
    else   % values of epsx and epsz are in slab
      for jj=1:Nx
        epsx(jj,ii)=epslb; % epsilon at both field locations are in the field
        epsz(jj,ii)=epslb; % epsilon at both field locations are in the field
      end
    end
  end
  
% determine time step based on background permittivity

  deltmax = 1/c/sqrt((1/delz)^2+(1/delx)^2) % timestep
  str = ['CFL time step = ',num2str(deltmax)];  disp(str); % maximum time step
  delt = input('Give time step  '); % chosen time step
  
  dtomudz = delt/(mu*delz);      % compute necessary constants
  dtomudx = delt/(mu*delx);      % compute necessary constants
  dtoepdz = delt/(epsilon*delz); % compute necessary constants
  dtoepdx = delt/(epsilon*delx); % compute necessary consta
  beta=1/(2*delz/c/delt + 1);
  alpha=beta*(2*delz/c/delt - 1);

% march in time
  
  fid = fopen('snapshots.txt', 'wt');
  
   for kk = 1:50 % 500 steps?

%    update source   ------------------------------------------------------

     t = kk*delt; % time is the step index times the time step
     
         if(( t >= 0) && (t <= BigT))  % source excitation if we are in the period
          Jx(Sx,Sz)=sin(2*pi*t/BigT)/eta/delx/delz; % sine wave
         else
          Jx(Sx,Sz)=0; % after the excitation time make the source zero
         end
     
%    update Hy at time kk+1/2  (Nx by Nz samples)  ------------------------

     for jj = 1:Nx
         Hymem(jj)=Hy(jj,Nzm1); % store next-to-last values for ABC (Nzm1=Nz-1)
     end
     
     % update equations for Hy, including boundaries at x = 1,Nx
     for ii = 1:Nzm1 % iterate over z    
       for jj = 1:Nx % iterate over x
         
         if(jj == 1) % zero Ez field at wall of guide jj % top PEC
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(Ez(jj+1,ii)          );
                           
         elseif(jj == Nx)  % zero Ez field at wall of guide jj+1 % BOTTOM PEC
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(           -Ez(jj,ii));
                           
         else % out in the middle of the mesh somewhere
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(Ez(jj+1,ii)-Ez(jj,ii));
                           
         end
       end
     end

%    update Hy at absorbing boundary Nz

       for jj = 1:Nx
         
         if(jj == 1) % zero Ez field at wall of guide jj
                
            Hy(jj,Nz)=alpha*Hy(jj,Nz)...
                     +beta*(Hy(jj,Nzm1)+Hymem(jj))...
                     +beta*(Ez(jj+1,Nz)          )/eta;
                           
         elseif(jj == Nx)  % zero Ez field at wall of guide jj+1
                
            Hy(jj,Nz)=alpha*Hy(jj,Nz)...
                     +beta*(Hy(jj,Nzm1)+Hymem(jj))...
                     +beta*(           -Ez(jj,Nz))/eta;
                           
         else % out in the middle of the mesh somewhere
                
            Hy(jj,Nz)=alpha*Hy(jj,Nz)...
                     +beta*(Hy(jj,Nzm1)+Hymem(jj))...
                     +beta*(Ez(jj+1,Nz)-Ez(jj,Nz))/eta;
                           
         end
       end

%    update Ex at time kk  (Nx by Nz-1 samples)  --------------------------

     for ii = 2:Nz  % don't update Ex(x,1), keep at zero for PEC
       for jj = 1:Nx
                         
         Ex(jj,ii)=Ex(jj,ii)-dtoepdz*(Hy(jj,ii)-Hy(jj,ii-1))/epsx(jj,ii)...
                            -(delt/epsilon)*Jx(jj,ii);
                           
       end 
     end
       
%    update Ez at time kk (Nx-1 by Nz samples)  ---------------------------

     for ii = 1:Nz     
       for jj = 2:Nx
                 
         Ez(jj,ii)=Ez(jj,ii)+dtoepdx*(Hy(jj,ii)-Hy(jj-1,ii))/epsz(jj,ii);  
            
       end
         
     end

% store snapshot down centerline in file 'fid'  ---------------------------

      str = ['result at time = ',num2str(t)];
      fprintf(fid,'%s \n',str);
      for ii=1:Nz
         z = delz * (ii-1);
         fprintf(fid,'%15.14g %15.14g\n',z,Ex(Nsamp,ii));
      end
      fprintf(fid,'%15.14g %15.14g\n',L,0.0);
      fprintf(fid,'\n\n');

   end
  
end

