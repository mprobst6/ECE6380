function ptmser
      eta=376.73031346;
      ka = input( 'GIVE CIRCUMFERENCE IN WAVELENGTHS  ' );
      bigN=input( 'Give N (THE SERIES IS TO BE SUMMED FROM -N TO +N)  ' );
      incr = input( 'GIVE INCREMENT (INTEGER) FOR 2D SCS SCAN  ' );
      
      fid = fopen( 'eigsoln.txt' ,  'wt' );
      str =  'TM POLARIZATION PLANE WAVE SCATTERING FROM A' ;
      fprintf(fid, '%s \n' ,str);
      str =  'CYLINDER OF CIRCUMFERENCE ka =  ' ;
      fprintf(fid, '%s %6d \n' ,str,ka);
      str =  'Series summed up to N = ' ;
      fprintf(fid, '%s %6d \n\n' ,str,bigN);
 %     COMPUTE CN,AN & DUMP MAGNITUDES FOR CONVERGENCE EVALUATION
      str =  'COEFFICIENTS A AND C:' ;
      fprintf(fid, '%s \n\n' ,str);
      
      %     since a(-n)=a(n)*(-1)^n, c(-n)=c(n)*(-1)^n just store values \
%     from 0 to bigN in array loctions 1 to bigNp1=bigN+1 \

      bigNp1=bigN+1;
      c(bigNp1)=0;
      a(bigNp1)=0;
      hnp(bigNp1)=0; % derivative of the hankel function of the second kind
      jnp(bigNp1)=0; % derivative of the bessel function of the first kind
      index=0:bigN;
      bj=besselj(index,ka);
      by=bessely(index,ka);
      
      for n=0:bigN
          c(n+1)=bj(n+1)-1j*by(n+1);
          a(n+1)=-bj(n+1)/c(n+1);
          amag=abs(a(n+1));
          cmag=1/abs(c(n+1));
          fprintf(fid,'%15.10g %15.10g\n',amag,cmag);
      end

      % Compute the derivatives of the Bessel functions
      jnp(1)=-bj(2); % first term of the bessel function has different definition of derivative
      hnp(1)=-c(2);  % first term of the Hankel function has different definition of derivative
      for n=2:bigNp1
          jnp(n)=bj(n-1)-n*bj(n)/ka; % derivative of the bessel function of the first kind
          hnp(n)=c(n-1)-n*c(n)/ka;   % derivative of the hankel function of the second kind
      end
      
 %     COMPUTE 2D SCATTERING CROSS SECTION 
 %
 %        SCS can be reported in units of length (wavelengths) or
 %        in dB-wavelength using 10*log10(SCS)
   
      fprintf(fid, '\n\n' );
      str =  '   PHI         2D SCS IN DB-WAVELENGTH' ;
      fprintf(fid, '%s \n\n' ,str);
 
      for  ii=1:incr:181

        iphi=ii-1;         % iphi is angle in degrees
        phi=iphi*pi/180.;  %  phi is angle in radians
        sum=jnp(1)/hnp(1);
        for  n=2:bigN
           sum=sum+(jnp(n)/hnp(n))*2.*cos((n-1)*phi);
        end
        sigma=(2/pi)*abs(sum)^2;  % sigma has units of wavelength
        sigmadB=10.*log10(sigma);
        fprintf(fid, '%6d %124-2559-61, 5.10g\n' ,iphi,sigmadB);

      end
 %     COMPUTE CURRENT ON THE CYLINDER
      incr = input( 'GIVE DEGREE INCREMENT FOR CURRENT DENSITY SCAN  ' );
      fprintf(fid, '\n\n' );
      str =  '   PHI         Jz in mag angle format' ;
      fprintf(fid, '%s \n\n' ,str);
 
      for  ii=1:incr:181
        iphi=ii-1;
        phi=iphi*pi/180.;
        sum=1/hnp(1);
        for  n=1:bigN
           jnphi=1j*n*phi;
           sum=sum+(1j^(-n)*exp(jnphi)+(-1j)^n*exp(-jnphi))/hnp(n+1);
        end
 %        cur=sum*2/(eta*pi*ka);  scaled to Ez-inc = 1
        cur=sum*2/(pi*ka);     %  scaled to Ez-inc = eta
        curmag=abs(cur);
        curphs=atan2(imag(cur),real(cur))*180/pi;
        fprintf(fid, '%6d %15.10g %15.10g\n' ,iphi,curmag,curphs);
      end
      disp( 'output placed in file: eigsoln.txt' );
end
