function plottrimeshQ

% load a triangular cell mesh of a 2D object and plot it
%
% Oct 25 2018   A. F. Peterson

% read mesh from file 'cylfil.txt.txt'

 nnodes = dlmread('cylfil.txt','', [0,0,0,0]);
 ncells = dlmread('cylfil.txt','', [0,1,0,1]);

 Node = dlmread('cylfil.txt','', [1,0,nnodes,1]);
 
 nstart=nnodes + 1;
 nend=nstart + ncells - 1;
 Element_to_Node = dlmread('cylfil.txt','', [nstart,0,nend,2]);
 
 % disp(Element_to_Node);
 
 figure
 trimesh( Element_to_Node, Node(:,1), Node(:,2) )

end

