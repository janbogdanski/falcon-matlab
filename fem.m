function fem_1D
% This is a simple 1D FEM program. The FEM
% solution is based on linear elements also called hat functions.The
% problem addressed is the extension of a bar under the action of applied
% forces. Use mesh parameters under the heading mesh of this code to change
% values. For example change the number of nodes to 2 to really see the
% difference between the exact and the FEM solution. In general, change
% constants the way you like.
%
% The second plot of stresses in the bar suggests that for each of the
% finite elements in the bar the solution (that is the slope of the
% extension) is a constant. This is illustrated with the help of horizontal
% lines. Each line also represents one element. It would be shown in later
% codes that the slope of the approximate solution (stresses) coincides
% with the slope of the exact solution exactly at the mid point of the
% element.

% This code is inspired from a manual written by Jack Chessa on FEM
% implementation in MATLAB.

% Written 26-Dec-2010
% Author: Husnain Inayat Hussain (husnain.inayat@gmail.com)

% FEM_1D

% 1. PreProcessor
% 2. Processor
% 3. PostProcessor

% ----------------   PREPROCESSOR   -------------------

% Preprocessor consists of defining the problem, material, elemental
% discretization, geometry, applied conditions such as forces and
% displacements.

% PROBLEM



% MATERIAL
E = 1; % Young's Modulus
A = 1; % Area of the bar
L = 1; % Length of the bar
f_o = 1;
u_o =0;
P = 10;

% MESH
num_nodes = 3;
num_elems = num_nodes - 1;
nodes = linspace (0,L,num_nodes);
elm_cnc = [1:1:num_elems ; 2:1:num_nodes]';
h_e = nodes (2) - nodes (1);

% ----------------   PROCESSOR   -------------------

% Stiffness Matrix
K_e = E*A/h_e * [1 -1;
                -1  1];
K = zeros (num_nodes, num_nodes);

F_e = f_o*h_e/2;
F = zeros (num_nodes, 1);

for i = 1 : size(elm_cnc,1)
    asm = elm_cnc (i,:);
    K(asm , asm) = K(asm , asm) + K_e;
    F(asm) = F(asm) + F_e;
end

F = F - K(:,1) * u_o;
F (1) = u_o;
F (end) = F (end) + P;

K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;

u = K\F;

% ----------------   POSTPROCESSOR   -------------------

% Plot the displacements and the stressed of the analytical and FEM
% solutions.

% FEM Solution & Plot
plot(nodes, u,'Marker','o','LineWidth',1,'LineStyle','--',...
    'Color',[1 0 0],...
    'DisplayName','FEM');

% Analytical Solution & Plot
x = linspace(0,L,100);
u1 = - f_o/(2*E*A)*x.^2 + (P + f_o*L)*x;
hold
plot(x, u1,'LineWidth',1,'DisplayName','Analytical');
h = legend('FEM','Analytical',2);
set(h,'Interpreter','none')
grid on

% Create title
title({'Extension of a Longitudianl Bar','Analytical vs.FEM Solution'},...
    'FontWeight','normal',...
    'FontSize',12);

% Create xlabel
xlabel({'Lenght (m)'},'FontWeight','normal','FontSize',11);

% Create ylabel
ylabel({'Extension (m)'},'FontWeight','normal','FontSize',11);


% Plot the stresses of the analytical and FEM solutions

% FEM Solution & Plot
figure
du = zeros (5,1);
plot(nan,nan)
hold

h1 = zeros (size(elm_cnc,1),1);

for i = 1 : size(elm_cnc,1)
    du(i) = (u(i+1) - u (i))/ (nodes(i+1) - nodes(i));
    asm = elm_cnc (i,:);
    h1(i) = plot (nodes(asm), repmat(du(i),1,length(nodes(asm))),...
        'Marker','o','LineWidth',1,'LineStyle','--','Color',[1 0 0],...
        'DisplayName','FEM');
end

% Analytical Solution & Plot
du1 = - f_o/(E*A)*x + (P + f_o*L);
h2 = plot (x, du1, 'LineWidth',1,'LineStyle','-',...
    'Color',[0 0 1], 'DisplayName','Analytical');
h = legend([h1(end) h2], 'FEM','Analytical', 1);
set(h,'Interpreter','none')
grid on

% Create title
title({'Stresses in a Longitudianl Bar','Analytical vs.FEM Solution'},...
    'FontWeight','normal',...
    'FontSize',12);

% Create xlabel
xlabel({'Lenght (m)'},'FontWeight','normal','FontSize',11);

% Create ylabel
ylabel({'Stress (N/m^2)'},'FontWeight','normal','FontSize',11);