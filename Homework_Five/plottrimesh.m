function plottrimesh

% load a triangular cell mesh of a 3D object and plot it
%
% March 30, 2016   A. F. Peterson

% read mesh from file 'trifile.txt'

 nnodes = dlmread('trifile.txt','', [0,0,0,0]);
 ncells = dlmread('trifile.txt','', [0,1,0,1]);

 Node = dlmread('trifile.txt','', [1,1,nnodes,2]);
 
 nstart=nnodes + 1;
 nend=nstart + ncells - 1;
 Element_to_Node = dlmread('trifile.txt','', [nstart,1,nend,3]);
 
 % disp(Element_to_Node);
 
 figure
 trimesh( Element_to_Node, Node(:,1), Node(:,2) )

end

