
n_nodes = 10;
n_dofs = n_nodes*2;
n_elements = n_nodes-1;
x = linspace(0, 1, n_nodes);



tspan = linspace(0,10, 1000);
y0 = zeros(2*n_dofs-4,1);

[t, y] = ode45(@f, tspan, y0); 


function dydt = f(t,y)


n_nodes = 10;
n_dofs = n_nodes*2;
n_elements = n_nodes-1;
x = linspace(0, 1, n_nodes);


L = 1000; %mm
rho = 0.006; %kg/mm
E = 210000; %N/mm2
I = 0.801*10^6; %mm4
A = 764; %mm^2

%Element stiffnes matrix
L = L/n_elements;
k = [12, -6*L, -12, -6*L; -6*L, 4*L^2, 6*L, 2*L^2; -12, 6*L, 12, 6*L; -6*L, 2*L^2, 6*L, 4*L^2;]*E*I/L^3;

%Element mass matrix
m = [156, 22*L, 54, -13*L; 22*L, 4*L^2, 13*L, -3*L^2; 54, 13*L, 156, -22*L; -13*L, -3*L^2, -22*L, 4*L^2;]*rho*A*L/420;

%Global stiffness matrix
n_dofs = n_nodes*2;
K = zeros(n_dofs, n_dofs);
for i=1:n_elements
    K(1+2*(i-1):1+2*(i-1)+3,1+2*(i-1):1+2*(i-1)+3) = K(1+2*(i-1):1+2*(i-1)+3,1+2*(i-1):1+2*(i-1)+3) + k;

end

%Global mass matrix
M = zeros(n_dofs, n_dofs);
for i=1:n_elements
    M(1+2*(i-1):1+2*(i-1)+3,1+2*(i-1):1+2*(i-1)+3) = M(1+2*(i-1):1+2*(i-1)+3,1+2*(i-1):1+2*(i-1)+3) + m;
end

%Force vector
F = zeros(n_dofs, 1); 
F(end-1) = 10000;

%Initial conditions
dz0 = zeros(n_nodes, 1);
dzt0 = zeros(n_nodes, 1);

d0 = zeros(n_dofs, 1);
dt0 = zeros(n_dofs, 1);
d0(1:2:end) = dz0;

F = F(3:end);
K = K(3:end,3:end);
M = M(3:end,3:end);

%Static solution
d = K\F;
d = [0; 0; d];
%plot(d(1:2:end))
    
    n_dofs = 20;
    dydt = zeros(2*n_dofs-4,1);
    dydt(1:n_dofs-2) = y(n_dofs-1:2*n_dofs-4);
    dydt(n_dofs-1:end) = M\(F - K *y(1:n_dofs-2));
end

