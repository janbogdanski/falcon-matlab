% beam.m
%
% Solves a linear elastic 2D beam problem ( plane stress or strain )
% with several element types.
%
% ^ y
% |
% ---------------------------------------------% | |
% | |
% ---------> x | 2c
% | |
% | L |
% ---------------------------------------------%
% with the boundary following conditions:
%
% u_x = 0 at (0,0), (0,-c) and (0,c)
% u_y = 0 at (0,0)
%
% t_x = y along the edge x=0
% t_y = P*(x^2-c^2) along the edge x=L
%
% ******************************************************************************
%
% This file and the supporting matlab files can be found at
% http://www.tam.northwestern.edu/jfc795/Matlab
%
% by Jack Chessa
% Northwestern University
% ******************************************************************************
clear
colordef black
state = 0;
% ******************************************************************************
% *** I N P U T ***
% ******************************************************************************
tic;
disp('************************************************')
disp('*** S T A R T I N G R U N ***')
disp('************************************************')
disp([num2str(toc),' START'])
% MATERIAL PROPERTIES
E0 = 10e7; % Young's modulus
nu0 = 0.30; % Poisson's ratio
% BEAM PROPERTIES
L = 16; % length of the beam
c = 2; % the distance of the outer fiber of the beam from the mid-line
% MESH PROPERTIES
elemType = 'Q9'; % the element type used in the FEM simulation; 'T3' is for a
% three node constant strain triangular element, 'Q4' is for
% a four node quadrilateral element, and 'Q9' is for a nine
% node quadrilateral element.
numy = 4; % the number of elements in the x-direction (beam length)
numx = 18; % and in the y-direciton.
plotMesh = 1; % A flag that if set to 1 plots the initial mesh (to make sure
% that the mesh is correct)
% TIP LOAD
P = -1; % the peak magnitude of the traction at the right edge
% STRESS ASSUMPTION
stressState='PLANE_STRESS'; % set to either 'PLANE_STRAIN' or "PLANE_STRESS'
% nuff said.
% ******************************************************************************
% *** P R E - P R O C E S S I N G ***
% ******************************************************************************
I0=2*c^3/3; % the second polar moment of inertia of the beam cross-section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELASTICITY MATRIX
if ( strcmp(stressState,'PLANE_STRESS') ) % Plane Strain case
C=E0/(1-nu0^2)*[ 1 nu0 0;
nu0 1 0;
0 0 (1-nu0)/2 ];
else % Plane Strain case
C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 0;
nu0 1-nu0 0;
0 0 1/2-nu0 ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
% Here we gnerate the finte element mesh (using the approriate elements).
% I won't go into too much detail about how to use these functions. If
% one is interested one can type - help 'function name' at the matlab comand
% line to find out more about it.
%
% The folowing data structures are used to describe the finite element
% discretization:
%
% node - is a matrix of the node coordinates, i.e. node(I,j) -> x_Ij
% element - is a matrix of element connectivities, i.e. the connectivity
% of element e is given by > element(e,:) -> [n1 n2 n3 ...];
%
% To apply boundary conditions a description of the boundaries is needed. To
% accomplish this we use a separate finite element discretization for each
% boundary. For a 2D problem the boundary discretization is a set of 1D elements.
%
% rightEdge - a element connectivity matrix for the right edge
% leftEdge - I'll give you three guesses
%
% These connectivity matricies refer to the node numbers defined in the
% coordinate matrix node.
disp([num2str(toc),' GENERATING MESH'])
switch elemType
case 'Q4' % here we generate the mesh of Q4 elements
nnx=numx+1;
nny=numy+1;
node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
case 'Q9' % here we generate a mehs of Q9 elements
nnx=2*numx+1;
nny=2*numy+1;
node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
inc_u=2;
inc_v=2*nnx;
node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];
element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
otherwise %'T3' % and last but not least T3 elements
nnx=numx+1;
nny=numy+1;
node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
node_pattern1=[ 1 2 nnx+1 ];
node_pattern2=[ 2 nnx+2 nnx+1 ];
inc_u=1;
inc_v=nnx;
element=[make_elem(node_pattern1,numx,numy,inc_u,inc_v);
make_elem(node_pattern2,numx,numy,inc_u,inc_v) ];
end
% DEFINE BOUNDARIES
% Here we define the boundary discretizations.
uln=nnx*(nny-1)+1; % upper left node number
urn=nnx*nny; % upper right node number
lrn=nnx; % lower right node number
lln=1; % lower left node number
cln=nnx*(nny-1)/2+1; % node number at (0,0)
switch elemType
case 'Q9'
rightEdge=[ lrn:2*nnx:(uln-1); (lrn+2*nnx):2*nnx:urn; (lrn+nnx):2*nnx:urn ]';
leftEdge =[ uln:-2*nnx:(lrn+1); (uln-2*nnx):-2*nnx:1; (uln-nnx):-2*nnx:1 ]';
edgeElemType='L3';
otherwise % same discretizations for Q4 and T3 meshes
rightEdge=[ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
leftEdge =[ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
edgeElemType='L2';
end
% GET NODES ON DISPLACEMENT BOUNDARY
% Here we get the nodes on the essential boundaries
fixedNodeX=[uln lln cln]'; % a vector of the node numbers which are fixed in
% the x direction
fixedNodeY=[cln]'; % a vector of node numbers which are fixed in
% the y-direction
uFixed=zeros(size(fixedNodeX)); % a vector of the x-displacement for the nodes
% in fixedNodeX ( in this case just zeros )
vFixed=zeros(size(fixedNodeY)); % and the y-displacements for fixedNodeY
numnode=size(node,1); % number of nodes
numelem=size(element,1); % number of elements
% PLOT MESH
if ( plotMesh ) % if plotMesh==1 we will plot the mesh
clf
plot_mesh(node,element,elemType,'g.-');
hold on
plot_mesh(node,rightEdge,edgeElemType,'bo-');
plot_mesh(node,leftEdge,edgeElemType,'bo-');
plot(node(fixedNodeX,1),node(fixedNodeX,2),'r>');
plot(node(fixedNodeY,1),node(fixedNodeY,2),'r^');
axis off
axis([0 L -c c])
disp('(paused)')
pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES
%
% Here we define the system data structures
% U - is vector of the nodal displacements it is of length 2*numnode. The
% displacements in the x-direction are in the top half of U and the
% y-displacements are in the lower half of U, for example the displacement
% in the y-direction for node number I is at U(I+numnode)
% f - is the nodal force vector. It's structure is the same as U,
% i.e. f(I+numnode) is the force in the y direction at node I
% K - is the global stiffness matrix and is structured the same as with U and f
% so that K_IiJj is at K(I+(i-1)*numnode,J+(j-1)*numnode)
disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])
U=zeros(2*numnode,1); % nodal displacement vector

f=zeros(2*numnode,1); % external load vector
K=sparse(2*numnode,2*numnode); % stiffness matrix
% a vector of indicies that quickly address the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal x-displacements
xs=1:numnode; % x portion of u and v vectors
ys=(numnode+1):2*numnode; % y portion of u and v vectors
% ******************************************************************************
% *** P R O C E S S I N G ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE EXTERNAL FORCES
% integrate the tractions on the left and right edges
disp([num2str(toc),' COMPUTING EXTERNAL LOADS'])
switch elemType % define quadrature rule
case 'Q9'
[W,Q]=quadrature( 4, 'GAUSS', 1 ); % four point quadrature
otherwise
[W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature
end
% RIGHT EDGE
for e=1:size(rightEdge,1) % loop over the elements in the right edge
sctr=rightEdge(e,:); % scatter vector for the element
sctrx=sctr; % x scatter vector
sctry=sctrx+numnode; % y scatter vector
for q=1:size(W,1) % quadrature loop
pt=Q(q,:); % quadrature point
wt=W(q); % quadrature weight
[N,dNdxi]=lagrange_basis(edgeElemType,pt); % element shape functions
J0=dNdxi'*node(sctr,:); % element Jacobian
detJ0=norm(J0); % determiniat of jacobian
yPt=N'*node(sctr,2); % y coordinate at quadrature point
fyPt=P*(c^2-yPt^2)/(2*I0); % y traction at quadrature point
f(sctry)=f(sctry)+N*fyPt*detJ0*wt; % scatter force into global force vector
end % of quadrature loop
end % of element loop
% LEFT EDGE
for e=1:size(leftEdge,1) % loop over the elements in the left edge
sctr=rightEdge(e,:);
sctrx=sctr;
sctry=sctrx+numnode;
for q=1:size(W,1) % quadrature loop
pt=Q(q,:); % quadrature point
wt=W(q); % quadrature weight
[N,dNdxi]=lagrange_basis(edgeElemType,pt); % element shape functions
J0=dNdxi'*node(sctr,:); % element Jacobian
detJ0=norm(J0); % determiniat of jacobian
yPt=N'*node(sctr,2);
fyPt=-P*(c^2-yPt^2)/(2*I0); % y traction at quadrature point
fxPt=P*L*yPt/I0; % x traction at quadrature poin

f(sctry)=f(sctry)+N*fyPt*detJ0*wt;
f(sctrx)=f(sctrx)+N*fxPt*detJ0*wt;
end % of quadrature loop
end % of element loop
% set the force at the nodes on the top and bottom edges to zero (traction free)
% TOP EDGE
topEdgeNodes = find(node(:,2)==c); % finds nodes on the top edge
f(topEdgeNodes)=0;
f(topEdgeNodes+numnode)=0;
% BOTTOM EDGE
bottomEdgeNodes = find(node(:,2)==-c); % finds nodes on the bottom edge
f(bottomEdgeNodes)=0;
f(bottomEdgeNodes+numnode)=0;
%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),' COMPUTING STIFFNESS MATRIX'])
switch elemType % define quadrature rule
case 'Q9'
[W,Q]=quadrature( 4, 'GAUSS', 2 ); % 4x4 Gaussian quadrature
case 'Q4'
[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
otherwise
[W,Q]=quadrature( 1, 'TRIANGULAR', 2 ); % 1 point triangural quadrature
end
for e=1:numelem % start of element loop
sctr=element(e,:); % element scatter vector
sctrB=[ sctr sctr+numnode ]; % vector that scatters a B matrix
nn=length(sctr);
for q=1:size(W,1) % quadrature loop
pt=Q(q,:); % quadrature point
wt=W(q); % quadrature weight
[N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
J0=node(sctr,:)'*dNdxi; % element Jacobian matrix
invJ0=inv(J0);
dNdx=dNdxi*invJ0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE B MATRIX
% _ _
% | N_1,x N_2,x ... 0 0 ... |
% B = | 0 0 ... N_1,y N_2,y ... |
% | N_1,y N_2,y ... N_1,x N_2,x ... |
% - -B=zeros(3,2*nn);
B(1,1:nn) = dNdx(:,1)';
B(2,nn+1:2*nn) = dNdx(:,2)';
B(3,1:nn) = dNdx(:,2)';
B(3,nn+1:2*nn) = dNdx(:,1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
end % of quadrature loop

end % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%
% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),' APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=fixedNodeX; % global indecies of the fixed x displacements
vdofs=fixedNodeY+numnode; % global indecies of the fixed y displacements
f=f-K(:,udofs)*uFixed; % modify the force vector
f=f-K(:,vdofs)*vFixed;
f(udofs)=uFixed;
f(vdofs)=vFixed;
K(udofs,:)=0; % zero out the rows and columns of the K matrix
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
% SOLVE SYSTEM
disp([num2str(toc),' SOLVING SYSTEM'])
U=K\f;
%******************************************************************************
%*** P O S T - P R O C E S S I N G ***
%******************************************************************************
%
% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don't go into too much detail - use help
% 'function name' to get more details.
disp([num2str(toc),' POST-PROCESSING'])
dispNorm=L/max(sqrt(U(xs).^2+U(ys).^2));
scaleFact=0.1*dispNorm;
fn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
figure(fn)
clf
plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,U(ys));
hold on
plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS
stress=zeros(numelem,size(element,2),3);
switch elemType % define quadrature rule
case 'Q9'
stressPoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0 ];
case 'Q4'
stressPoints=[-1 -1;1 -1;1 1;-1 1];
    otherwise
    
        stressPoints=[0 0;1 0;0 1];
end
for e=1:numelem % start of element loop
sctr=element(e,:);
sctrB=[sctr sctr+numnode];
nn=length(sctr);
for q=1:nn
pt=stressPoints(q,:); % stress point
[N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
J0=node(sctr,:)'*dNdxi; % element Jacobian matrix
invJ0=inv(J0);
dNdx=dNdxi*invJ0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE B MATRIX
B=zeros(3,2*nn);
B(1,1:nn) = dNdx(:,1)';
B(2,nn+1:2*nn) = dNdx(:,2)';
B(3,1:nn) = dNdx(:,2)';
B(3,nn+1:2*nn) = dNdx(:,1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
strain=B*U(sctrB);
stress(e,q,:)=C*strain;
end
end % of element loop
stressComp=1;
figure(fn)
clf
plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,stress(:,:,stressComp));
hold on
plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED STRESS PLOT, BENDING COMPONENT')
%print(fn,'-djpeg90',['beam_',elemType,'_sigma',num2str(stressComp),'.jpg'])
disp([num2str(toc),' RUN FINISHED'])
% ***************************************************************************
% *** E N D O F P R O G R A M ***
% ***************************************************************************
disp('************************************************')
disp('*** E N D O F R U N ***')
disp('************************************************')