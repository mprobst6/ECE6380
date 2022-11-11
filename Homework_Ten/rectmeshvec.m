function rectmeshvec
%
% Generate triangular-cell model of rectangular cavity with boundary
% edges at the end of the list
%
% January 19, 2013, adapted from quadmesh version April 4, 2012
% A. F. Peterson
% pointer consistency issues Oct 2020

% initialization

 global nedgs pedtop pedton;

 nsize = 5000;
 x = zeros(nsize);
 y = zeros(nsize);
 ppaton = zeros(nsize,3);
 ppatoe = zeros(nsize,3);
 pedtop = zeros(nsize,2);
 pedton = zeros(nsize,2);
 ppaton = zeros(nsize,6);
 point2 = zeros(nsize);

%  input description of cavity

 bigx = input('give cavity x dimension  ');
 bigy = input('give cavity y dimension  ');

% err = input('give real part of epsilon-r  ');
% eri = input('give imaginary part of epsilon-r  ');
% er(1)=err+1j*eri;

 ncx = input('give number of cells along x  ');
 ncy = input('give number of cells along y  ');

%--------------------------------------------------------------------
% part one: generate nodes and pointer 'ppaton'
%--------------------------------------------------------------------

	  xlow=-0.5*bigx;
	  ylow=-0.5*bigy;
	  delx=bigx/ncx;
	  dely=bigy/ncy;
%     initialize the cell counter

      ncells=0;

%     initialize the nodes along left edge of plate

      for iy=0:ncy
	     x(iy+1)=xlow;
	     y(iy+1)=ylow+iy*dely;
      end
      nodes=ncy+1;

	  nstart=1;
%	  nend=nodes;

%     loop through plate by columns of cells

      for ix=1:ncx

%     add a new column of cell nodes

	     for iy=0:ncy
	        in=nodes+1+iy;
	        x(in)=xlow+delx;
	        y(in)=ylow+iy*dely;
         end
	     nright=nodes+1;
	     nodes=nodes+ncy+1;

%     build connectivity -- cell to node
%
%     -- nodes on left run from 'nstart' to 'nend'
%     -- nodes on right run from 'nright' to 'nodes'

         for iy=1:ncy
	        ncells=ncells+1;
	        ppaton(ncells,1)=nstart;
	        ppaton(ncells,2)=nright;
	        ppaton(ncells,3)=nright+1;

	        ncells=ncells+1;
	        ppaton(ncells,1)=nstart;
	        ppaton(ncells,2)=nright+1;
	        ppaton(ncells,3)=nstart+1;

	        nstart=nstart+1;
	        nright=nright+1;
         end

         nstart=nodes-ncy;
		 xlow=xlow+delx;
      end
      
%--------------------------------------------------------------------
% part two:  generate pointers associated with the edges:
%            pedton & pedtop
%--------------------------------------------------------------------

      nedgs=0;

%   initialize pointers to -1 to flag edge-of-model

      for ii=1:nsize
         pedtop(ii,2) = -1;
      end

%   each patch has three associated edges -- add edges to database
%   if necessary and construct pointer 'ppatoe'

      for ii=1:ncells
         edgenum = adedg(ii,ppaton(ii,1),ppaton(ii,2));
         edgenum = adedg(ii,ppaton(ii,2),ppaton(ii,3));
         edgenum = adedg(ii,ppaton(ii,3),ppaton(ii,1));
      end
 
      str=['number of edges in model = ',num2str(nedgs)];
      disp(str);
      if(nedgs > nsize) 
	     disp('WARNING:  nedgs exceeds limit of arrays ')
      end

%--------------------------------------------------------------------
%  part three -- sort edges to put boundary edges at end of list
%--------------------------------------------------------------------

 nbound=0;
 ninter=0;
 for ii=1:nedgs
     if(pedtop(ii,2) == -1)
         nbound=nbound+1;
         point2(ii)=3000+nbound;
     else
         ninter=ninter+1;
         point2(ii)=ninter;
     end
 end
 for ii=1:nedgs
     if(point2(ii) > 3000)
     point2(ii)=point2(ii)-3000+ninter;
     end
 end
 
%   point2 points from original edge # to new edge # (interior edges first)
%
%   rearrange other pointer arrays

% for ii=1:ncells
%     for jj=1:3
%     index = ppatoe(ii,jj);
%     ppatoe(ii,jj)=point2(index);
%     end
% end
 
 pedton2=pedton;
 pedtop2=pedtop;
 for ii=1:nedgs
     for jj=1:2
         index=point2(ii); % index is the new edge number
         pedton(index,jj)=pedton2(ii,jj);
         pedtop(index,jj)=pedtop2(ii,jj);
     end
 end

%--------------------------------------------------------------------
% part five: clean up 'ppaton' so nodes appear in CCW orientation
%--------------------------------------------------------------------

  for ii=1:ncells
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);
      xmid=(x(n1)+x(n2)+x(n3))/3;
      ymid=(y(n1)+y(n2)+y(n3))/3;
      
      phi1=atan2(y(n1)-ymid,x(n1)-xmid);
      if(phi1 < 0) 
          phi1=phi1+2*pi;
      end
      
      phi2=atan2(y(n2)-ymid,x(n2)-xmid);
      if(phi2 < 0) 
          phi2=phi2+2*pi;
      end
      
      phi3=atan2(y(n3)-ymid,x(n3)-xmid);
      if(phi3 < 0) 
          phi3=phi3+2*pi;
      end
      
      if(phi1 < phi2)
        if(phi1 < phi3)
          if(phi2 > phi3)
%              disp('changing order in PPATON');
              ntemp=n2;
              n2=n3;
              n3=ntemp;
          end
        end
      else
          if(phi2 > phi3)
%              disp('changing order in PPATON');
              ntemp=n1;
              n1=n3;
              n3=ntemp;
          elseif(phi3 > phi1)
%              disp('changing order in PPATON');
              ntemp=n1;
              n1=n2;
              n2=ntemp;
          end
      end
      ppaton(ii,1)=n1;
      ppaton(ii,2)=n2;
      ppaton(ii,3)=n3;      
  end

%--------------------------------------------------------------------
%  part six -- compute cell to edge pointer
%--------------------------------------------------------------------

  ppatoe(ncells,nedgs)=0;

       for ii=1:ncells
         iedg=0;

%        find all occurrances of patch 'i' in pointer 'pedtop'
 
         for jj=1:nedgs
           if(pedtop(jj,1) == ii)
             iedg=iedg+1;
             ppatoe(ii,iedg)=jj;
           elseif(pedtop(jj,2) == ii)
             iedg=iedg+1;
             ppatoe(ii,iedg)=jj;
           end
         end
         if(iedg ~= 3)
          disp('error in patch to edge connectivity for patch');
          disp(ii);
         end
       end

% clean up 'ppatoe' so edges are synchronized with nodes in 'ppaton'

  for ii=1:ncells
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);

      pstart(1:3)=ppatoe(ii,1:3);
      pnew=zeros(3,1);
      
%  find edge 1 between n2 and n3

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n2) || (nd1 == n3))
          if((nd2 == n3) || (nd2 == n2))
              pnew(1)=ppatoe(ii,jj);
          end
      end
    end
%  find edge 2 between n3 and n1

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n3) || (nd1 == n1))
          if((nd2 == n1) || (nd2 == n3))
              pnew(2)=ppatoe(ii,jj);
          end
      end
    end
      
%  find edge 3 between n1 and n2

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n1) || (nd1 == n2))
          if((nd2 == n2) || (nd2 == n1))
              pnew(3)=ppatoe(ii,jj);
          end
      end
    end
    
    for jj=1:3
        if(pnew(jj) == 0)
            disp('error in edge sorting for patch ');
            disp(ii);
            break
        else
            ppatoe(ii,jj)=pnew(jj);   
        end
    end
  end

%--------------------------------------------------------------------
%  part seven -- compute edge lengths
%--------------------------------------------------------------------

  min=99999;
  max=0;
  sum=0;
  length(nedgs)=0;
  for ii=1:nedgs
      n1=pedton(ii,1);
      n2=pedton(ii,2);
      length(ii)=sqrt((x(n2)-x(n1))^2+(y(n2)-y(n1))^2);
      if(length(ii) > max)
          max=length(ii);
      end
      if(length(ii) < min)
          min=length(ii);
      end
      sum=sum+length(ii);
  end
  average=sum/nedgs;
  
  disp('min, max, and average edge lengths:');
  disp(min);  disp(max);  disp(average);
  disp('ratio of max length to min length:');
  sum=max/min;
  disp(sum);

%--------------------------------------------------------------------
%  part eight -- compute interior angles
%--------------------------------------------------------------------

  nptchs=ncells;
  ang(nptchs,3)=0;
  for ii=1:nptchs
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);
      
      v12x=x(n2)-x(n1);
      v12y=y(n2)-y(n1);
      mag12=sqrt(v12x^2+v12y^2);
      
      v13x=x(n3)-x(n1);
      v13y=y(n3)-y(n1);
      mag13=sqrt(v13x^2+v13y^2);
      
      v32x=x(n2)-x(n3);
      v32y=y(n2)-y(n3);
      mag32=sqrt(v32x^2+v32y^2);
      
%     ang1 is between edge 12 and edge 13

      dot=v12x*v13x + v12y*v13y;
      ang(ii,1)=acos(dot/(mag12*mag13))*180/pi;
      
%     ang2 is between edge 21 and edge 23

      dot=v12x*v32x + v12y*v32y;
      ang(ii,2)=acos(dot/(mag12*mag32))*180/pi;
      
%     ang3 is between edge 32 and edge 31

      dot=-(v32x*v13x + v32y*v13y);
      ang(ii,3)=acos(dot/(mag32*mag13))*180/pi;
      
%      str=[num2str(ang(ii,1)),'  ',num2str(ang(ii,2)),'  ',...
%           num2str(ang(ii,3))];
%      disp(str);
  end

  min=99999;
  max=0;
  for ii=1:nptchs
      for jj=1:3
          angle=ang(ii,jj);
          if(angle < min)
              min=angle;
          end
          if(angle > max)
              max=angle;
          end
      end
  end
  
  disp('min and max interior angles:');
  disp(min);  disp(max);
  disp('ratio of max angle to min angle:');
  sum=max/min;
  disp(sum);

%--------------------------------------------------------------------
%  part four -- dump data to file
%--------------------------------------------------------------------

 fid = fopen('cylfil.txt', 'wt');

 disp('CYLFIL.TXT contains nnods,nptchs,nedgs,ninn');
 disp('                    x,y coordinates');
 disp('                    ppaton(1:nptchs,3)');
 disp('                    ppatoe(1:nptchs,3)');
 disp('                    er(1:nptchs)');
 disp('                    pedtop(1:nedgs,2)');
 disp('                    pedton(1:nedgs,2)');
 disp('  ');
 disp('Run PLOTTRIMESH.M to generate a plot of the mesh');

 fprintf(fid,'%6d %6d %6d %6d\n',nodes,ncells,nedgs,ninter);

 for ii=1:nodes
   fprintf(fid,'%6d %15.14g %15.14g\n',ii, x(ii),y(ii));
 end

 for ii=1:ncells
   fprintf(fid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 for ii=1:ncells
   fprintf(fid,'%6d %6d %6d %6d\n',ii,ppatoe(ii,1),ppatoe(ii,2),...
       ppatoe(ii,3));
 end
 
 err = 1;
 for ii=1:ncells
  fprintf(fid,'%6d %15.14g\n',ii,err);
 end

 for ii=1:nedgs
  fprintf(fid,'%6d %6d %6d\n',ii,pedtop(ii,1),pedtop(ii,2));
 end
 
 for ii=1:nedgs
   fprintf(fid,'%6d %6d %6d\n',ii,pedton(ii,1),pedton(ii,2));
 end

% make copy in file 'trifile.txt' for plotter
 
 gid = fopen('trifile.txt', 'wt');
 
  fprintf(gid,'%6d %6d\n',nodes,ncells);

 for ii=1:nodes
   fprintf(gid,'%6d %15.14g %15.14g\n',ii,x(ii),y(ii));
 end
  
 for ii=1:ncells
   fprintf(gid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 disp('data placed in file: cylfil.txt');
 disp(' ');
 disp('number of nodes:   ');
 disp(nodes);
 disp('number of edges:   ');
 disp(nedgs);
 disp('number of patches: ');
 disp(ncells);

end

%
%--------------------------------------------------------------------
%     subroutine to create edge-related database
%--------------------------------------------------------------------
%

function[edgenum] = adedg(npat,n1,n2)

%  add the edge between nodes 'n1' and 'n2' to the database
%
%     npat: one of the patches associated with this edge
%     n1,n2: nodes associated with this edge
%     nedgs: number of edges already in the pointers

%     for new edge, add new row to pointers 'pedtop' and 'pedton'

%     if this edge is already in the database, add new info
%     about patch 'np' to the existing pointers

      global nedgs pedtop pedton;

%    find out if the edge is already in the database

      edgenum = 0;
      for k=1:nedgs
        if pedton(k,1) == n1 || pedton(k,2) == n1
          if pedton(k,1) == n2 || pedton(k,2) == n2
          edgenum=k;
          end
        end
      end

      if(edgenum == 0)           %  new edge

          nedgs=nedgs+1;
          pedton(nedgs,1)=n1;
          pedton(nedgs,2)=n2;
          pedtop(nedgs,1)=npat;
          edgenum=nedgs;
          
         else               %  edge already in database

          pedtop(edgenum,2)=npat;

      end
end
      

