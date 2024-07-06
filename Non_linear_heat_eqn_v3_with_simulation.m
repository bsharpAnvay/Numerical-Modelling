% Iteration 2. In the following code, a function is defined which has
% equations written only for unknown nodes and not boundary nodes.
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
prompt = sprintf("Set values of k in %dx%d matrix: ", b, a);
k = input(prompt);
old_phi = ones(b,a);
% Create a VideoWriter object
video_writer = VideoWriter('heatmap_video.avi', 'Uncompressed AVI');
open(video_writer);

% Create figure for heatmap
figure;

% Define color axis limits
caxis_limits = [-1, 1];

function F = matrixEquations(vars,t,dx,dy,dt,theta,k,l,old_phi)
a=(l/dx)+1;              % number of points in a row
b=(l/dy)+1;              % number of points in a column
xsq=(dx)^2;
ysq=(dy)^2;
phi = zeros(b,a);
for i=[1 b]
    for j=1:a
        phi(i,j)=1;
    end
end
for i=1:b
    j=a;
    phi(i,j)=1;
end
for i=1:b
    j=1;
    phi(i,j)=5*sin(pi/2+t*dt);
end
phi(2:b-1,2:a-1) = reshape(vars, [b-2, a-2]);  % Reshape vars into a (b-2)x(a-2) matrix
F = zeros((a-2)*(b-2), 1);            % initialising F
for i =2:b-1
    for j = 2:a-1
        F((a-2)*(i-2)+j-1) = (phi(i,j)- old_phi(i,j))/dt -theta*((k(i,j+1)*(phi(i,j+1))^2-2*k(i,j)*(phi(i,j))^2 + k(i,j-1)*(phi(i,j-1))^2)/xsq + (k(i+1,j)*(phi(i+1,j))^2-2*k(i,j)*(phi(i,j))^2 + k(i-1,j)*(phi(i-1,j))^2)/ysq) + (theta-1)*((k(i,j+1)*(old_phi(i,j+1))^2-2*k(i,j)*(old_phi(i,j))^2 + k(i,j-1)*(old_phi(i,j-1))^2)/xsq + (k(i+1,j)*(old_phi(i+1,j))^2-2*k(i,j)*(old_phi(i,j))^2 + k(i-1,j)*(old_phi(i-1,j))^2)/ysq);
    end
end
end

initialGuess = ones((a-2)*(b-2), 1);  % Initial guess as a vector
options = optimoptions('fsolve', 'Display', 'iter');  % Solve the equations numerically using fsolve. Option to display iteration details

for t=1:15
      % Plot the heatmap
    heatmap(old_phi);  % Plot the heatmap
    colormap('jet');  % Set colormap (e.g., 'jet', 'hot', 'cool', etc.)
    colorbar;         % Add a colorbar
   

     % Record the frame to the video multiple times
    for frameRepeat = 1:100  % Repeat each frame 5 times (adjust as needed)
        frame = getframe(gcf);  % Capture current figure
        writeVideo(video_writer, frame);  % Write frame to video
    end
    
    % Clear the figure for the next frame
    clf;
    [sol, fval, exitflag] = fsolve(@(vars)matrixEquations(vars,t,dx,dy,dt,theta,k,l,old_phi), initialGuess, options);
    old_phi(2:b-1,2:a-1) = reshape(sol, [b-2, a-2]);
    for i=1:b
        j=1;
        old_phi(i,j)=5*sin(pi/2+t*dt);
    end
end
phi=old_phi;
% disp(reshape(sol, [b-2, a-2]));  % Reshape the solution to a 10x10 matrix
% disp(phi);