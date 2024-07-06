% Iteration 1. In the following code, a function is defined which has
% equations written for every single node. But its not necessary to write
% equations for boundary points since their values are known for all time.
l = input("Enter value of side length of plate : ");
dx = l/floor(l/input("Enter value of delta x : "));
dy = l/floor(l/input("Enter value of delta y : "));
dt = input("Enter value of delta t : ");
theta = input("Enter value of theta : ");
a=(l/dx)+1;              % number of points in a row
b=(l/dy)+1;              % number of points in a column
% l = 4;
% dx =1;
% dy =1;
% dt =.2;
% theta = .5;
% k=[1 1 2 3 1; 1 3 4 7 1; 1 8 6 4 1; 1 1 1 1 1; 2 3 4 2 1];
% xsq=(dx)^2;
% ysq=(dy)^2;
prompt = sprintf("Set values of k in %dx%d matrix: ", a, b);
k = input(prompt);

% phi = sym('x', [a b]);  % Creates a 10x10 symbolic matrix with variables x1_1, x1_2, ..., x10_10
old_phi = ones(a,b);
function F = matrixEquations(vars,t,dx,dy,dt,theta,k,l,old_phi)
a=(l/dx)+1;              % number of points in a row
b=(l/dy)+1;              % number of points in a column
xsq=(dx)^2;
ysq=(dy)^2;

phi = reshape(vars, [a, b]);  % Reshape vars into a axb matrix
F = zeros(a*b, 1);            % Preallocate the function result
for i =2:b-1
    for j = 2:a-1
        F(a*(i-1)+j) = (phi(i,j)- old_phi(i,j))/dt -theta*((k(i,j+1)*(phi(i,j+1))^2-2*k(i,j)*(phi(i,j))^2 + k(i,j-1)*(phi(i,j-1))^2)/xsq + (k(i+1,j)*(phi(i+1,j))^2-2*k(i,j)*(phi(i,j))^2 + k(i-1,j)*(phi(i-1,j))^2)/ysq) + (theta-1)*((k(i,j+1)*(old_phi(i,j+1))^2-2*k(i,j)*(old_phi(i,j))^2 + k(i,j-1)*(old_phi(i,j-1))^2)/xsq + (k(i+1,j)*(old_phi(i+1,j))^2-2*k(i,j)*(old_phi(i,j))^2 + k(i-1,j)*(old_phi(i-1,j))^2)/ysq);
    end
end
for i=[1 b]
    for j=1:a
        F(a*(i-1)+j) = phi(i,j)-1;
    end
end
for i=1:b
    j=a;
    F(a*(i-1)+j) = phi(i,j)-1;
end
for i=1:b
    j=1;
    F(a*(i-1)+j) = phi(i,j)-sin(pi/2+t*dt);
end
end

initialGuess = ones(a*b, 1);  % Initial guess as a vector
options = optimoptions('fsolve', 'Display', 'iter');  % Solve the equations numerically using fsolve. Option to display iteration details

for t=1:10
    [sol, fval, exitflag] = fsolve(@(vars)matrixEquations(vars,t,dx,dy,dt,theta,k,l,old_phi), initialGuess, options);
    old_phi = reshape(sol, [a, b]);
end
disp(reshape(sol, [a, b]));  % Reshape the solution to a 10x10 matrix
