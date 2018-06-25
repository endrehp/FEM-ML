clear all
close all

n_nodes = 10;
n_dofs = n_nodes*2;
n_elements = n_nodes-1;
xs = linspace(0, 1, n_nodes);
xs4 = linspace(0,1, n_nodes*4);
dx = 1/n_elements;

L = 1000; %mm
rho = 0.006; %kg/mm
E = 210000; %N/mm2
I = 0.801*10^6; %mm4
A = 764; %mm^2
alpha1= 0.064; alpha2= 4.72e-5; %proportional damping coefficients
gamma= 0.5; beta= 0.25; %Newmark coefficients 
dt = 0.001; %timestep length
tf = 3;
TOL = 0.01;
error = 1;

%Shape functions
PSI(:,1)= 1 - (3/h^2)*xs.^2 + (2/h^3)*xs.^3;
PSI(:,2)= h*((xs/h) -(2/h^2)*xs.^2 + (1/h^3)*xs.^3);
PSI(:,3)= (3/h^2)*xs.^2 - (2/h^3)*xs.^3;
PSI(:,4)= h*(-(1/h^2)*xs.^2 + (1/h^3)*xs.^3);

PSIP(:,1)= -6*(xs/h^2) + (6/h^3)*xs.^2;
PSIP(:,2)= 1 - 4*(xs/h) + (3/h^2)*xs.^2;
PSIP(:,3)= 6*(xs/h^2) - (6/h^3)*xs.^2;
PSIP(:,4)= - 2*(xs/h) + (3/h^2)*xs.^2;

PSIP_f2_4(:,1)= -6*(xs4/h^2) + (6/h^3)*xs4.^2; %FOR EVALUATION OF f2 ONLY
PSIP_f2_4(:,2)= 1 - 4*(xs4/h) + (3/h^2)*xs4.^2;
PSIP_f2_4(:,3)= 6*(xs4/h^2) - (6/h^3)*xs4.^2;
PSIP_f2_4(:,4)= - 2*(xs4/h) + (3/h^2)*xs4.^2;

PSIPP(:,1)= -6/h^2 + (12/h^3)*xs;
PSIPP(:,2)= -4/h + (6/h^2)*xs;
PSIPP(:,3)= 6/h^2 - (12/h^3)*xs;
PSIPP(:,4)= -2/h + (6/h^2)*xs;

PSIPPP(:,1)= (12/h^3)*ones(1,(h/deltax)+1);
PSIPPP(:,2)= (6/h^2)*ones(1,(h/deltax)+1);
PSIPPP(:,3)= -(12/h^3)*ones(1,(h/deltax)+1);
PSIPPP(:,4)= (6/h^2)*ones(1,(h/deltax)+1); 

%Element stiffnes matrix
L = L/n_elements;
ke = [12, -6*L, -12, -6*L; -6*L, 4*L^2, 6*L, 2*L^2; -12, 6*L, 12, 6*L; -6*L, 2*L^2, 6*L, 4*L^2;]*E*I/L^3;

%Element mass matrix
me = [156, 22*L, 54, -13*L; 22*L, 4*L^2, 13*L, -3*L^2; 54, 13*L, 156, -22*L; -13*L, -3*L^2, -22*L, 4*L^2;]*rho*A*L/420;
ce = alpha1*me + alpha2*ke;

%Global linear matrices
n_dofs = n_nodes*2;
K = zeros(n_dofs, n_dofs);
M = zeros(n_dofs, n_dofs);
C = zeros(n_dofs, n_dofs);
for j=1:n_elements
    K(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) = K(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) + ke;
    M(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) = M(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) + me;
    C(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) = C(1+2*(j-1):1+2*(j-1)+3,1+2*(j-1):1+2*(j-1)+3) + ce;
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
C = C(3:end,3:end);

%Linear Newmark Matrices
A1 = M + gamma*dt*C + beta*dt^2*K;
A2 = -2*M + (1-2*gamma)*dt*C + (0.5-2*beta+gamma)*dt^2*K;
A3 = M -(1-gamma)*dt*C + (0.5+eta-gamma)*dt^2*K;


%Main loop
UL = zeros(tf/dt,n_dofs);
ULr = zeros(tf/dt, n_dofs-2);
UNLr = zeros(tf/dt, n_dofs-2);

for j=1:tf/dt

    if j == 1 
       
        U0 = zeros(n_dofs-2,1);
        U1 = zeros(n_dofs-2,1);
        F0 = F;
        F1 = F;
        
    else
        U0= ULr(:,j-1); 
        U1= ULr(:,j);
        F0= F; 
        F1= F; 
    end
    
    %FLin(:,j+1)= PHI*FD*cos(THETA*time(j+1));
    F2= F;
    F= dt^2*(beta*F2 + (0.5-2*beta+gamma)*F1 + (0.5+beta-gamma)*F0; 
    U2 = inv(A1)*(A2*U1-A3*U0 + dt^2*F);
    
    ULr(j+1,:) = U2;
    UL(j+1,3:end) = U2;

    
    %Nolinear loop
    eps= 10^5; 
    counter= 0; 
    
    while TOL < error
        

    %first nonlinear stiffness matrix
     a1= 1; 
     b0= 0;
        for o=1:N
            for p= 1:(1/dx)+1
                 WWP(p)= PSIP(p,1)*UL(a1,j+1) + PSIP(p,2)*UL(a1+1,j+1)+...
                 PSIP(p,3)*UL(a1+2,j+1) + PSIP(p,4)*UL(a1+3,j+1);
                 WWPP(p)= PSIPP(p,1)*UL(a1,j+1)+ PSIPP(p,2)*UL(a1+1,j+1)+...
                 PSIPP(p,3)*UL(a1+2,j+1) + PSIPP(p,4)*UL(a1+3,j+1);
                 WWPPP(p)= PSIPPP(p,1)*UL(a1,j+1)+ PSIPPP(p,2)*...
                 UL(a1+1,j+1)+PSIPPP(p,3)*UL(a1+2,j+1) + ...
                 PSIPPP(p,4)*UL(a1+3,j+1);
                 WP(p+b0)= WWP(p);
                 WPP(p+b0)= WWPP(p);
                 WPPP(p+b0)= WWPPP(p);
            end
            
            a1= a1+2;
            b0= b0 + 1/dx;
        end
     end
%Calculation of ki


end
