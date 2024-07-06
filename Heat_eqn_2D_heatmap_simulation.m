% Consider a square plate of length l. The temperature of its 3 edges is set to 1 at all times. The y-axis edge
% is varying sinusoidally with time.

l = input("Enter value of side length of plate : ");
dx = l/floor(l/input("Enter value of delta x : "));
dy = l/floor(l/input("Enter value of delta y : "));
dt = input("Enter value of delta t : ");
theta = input("Enter value of theta : ");
% l = 4;
% dx =1;
% dy =1;
% dt =.2;
% theta = .5;
xsq=1/((dx)^2);
ysq=1/((dy)^2);
a=(l/dx)-1;              % number of points in a row
b=(l/dy)-1;              % number of points in a column
A=zeros(a*b,a*b);        % initialising a matrix of the order a*b for finding temperatures of a*b number of points.

phi=zeros(a*b,1);         % setting values of phi for all points to be 1.
new_phi=zeros(a*b,1);     % initialising new_phi column matri
B=zeros(a*b,1);          % initialising the matrix B of A*T(t+1)=B

prompt = sprintf("Set values of k in %dx%d matrix: ", a, b);
k = input(prompt);
% k=[1 2 3; 3 4 7; 8 6 4];

% Creating VideoWriter object
video_writer = VideoWriter('heatmap_video.avi', 'Uncompressed AVI');
open(video_writer);

% Creating figure for heatmap
figure;

% Defining color axis limits
caxis_limits = [-1, 1];

for t=1:15

    final_temp=zeros(b+2,a+2);
    for p=1:b+2
    final_temp(p,1)=sin(t*dt);
    end
    for w=2:b+1
    A=phi((w-1-1)*a+1:((w-1)*a),1);
    final_temp(w,2:a+1)=A';
    end

    % Plot the heatmap
    heatmap(final_temp);  % Ploting the heatmap
    colormap('jet');  
    colorbar;         
   
    for frameRepeat = 1:10 
        frame = getframe(gcf);  
        writeVideo(video_writer, frame);  
    end
    
    % Clear the figure for the next frame
    clf;

    % making the matrix A and B of A*T(t+1)=B
    for i=2:b-1          % considering equations of all points except 1st and last row
        for j=1:a
            A((i-1)*a+j,(i-1)*a+j) = -2*k(i,j)*theta*dt*(xsq+ysq) - 1;
            A((i-1)*a+j,(i-2)*a+j) = k(i,j)*theta*dt*ysq;  % upper point
            A((i-1)*a+j,(i)*a+j) = k(i,j)*theta*dt*ysq;    % lower point
            if j==a
                A((i-1)*a+j,(i-1)*a+j-1) = k(i,j)*theta*dt*xsq; % left point
                B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+0) + ysq*(phi((i-2)*a+j)+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - theta*k(i,j)*dt*xsq;
                % subtracted the constant term from LHS
            elseif j==1
                A((i-1)*a+j,(i-1)*a+j+1) = k(i,j)*theta*dt*xsq; % right point
                B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j)+sin(((t-1)*dt)+0)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - theta*k(i,j)*dt*xsq*sin(((t)*dt)+0);
                % subtracted time dependent term from LHS
            else
                A((i-1)*a+j,(i-1)*a+j-1) = k(i,j)*theta*dt*xsq;
                A((i-1)*a+j,(i-1)*a+j+1) = k(i,j)*theta*dt*xsq;
                B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
            end
        end
    end
    % filling matrix A for 1st and last row of heat plate
    for j=2:a-1
        for i=[1,b]
            A((i-1)*a+j,(i-1)*a+j) = -2*k(i,j)*theta*dt*(xsq+ysq) - 1;
            A((i-1)*a+j,(i-1)*a+j-1) = k(i,j)*theta*dt*xsq;
            A((i-1)*a+j,(i-1)*a+j+1) = k(i,j)*theta*dt*xsq;
            if i==1
                A((i-1)*a+j,(i)*a+j) = k(i,j)*theta*dt*ysq;
                B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(1+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
            elseif i==b
                A((i-1)*a+j,(i-2)*a+j) = k(i,j)*theta*dt*ysq;
                B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+1)-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
            end
        end
    end
    % filling matrix A for corners of heat plate
    for j=[1,a]
        for i=[1,b]
            A((i-1)*a+j,(i-1)*a+j) = -2*k(i,j)*theta*dt*(xsq+ysq) - 1;
            if j==a
                A((i-1)*a+j,(i-1)*a+j-1) = k(i,j)*theta*dt*xsq;
                % B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-k(i,j)*theta*dt*ysq- k(i,j)*theta*dt*xsq;
                if i==1
                    A((i-1)*a+j,(i)*a+j) = k(i,j)*theta*dt*ysq;
                    B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+0) + ysq*(0+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
                elseif i==b
                    A((i-1)*a+j,(i-2)*a+j) = k(i,j)*theta*dt*ysq;
                    B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+0) + ysq*(phi((i-2)*a+j)+0)-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
                end
            end
            if j==1
                A((i-1)*a+j,(i-1)*a+j+1) = k(i,j)*theta*dt*xsq;
                % B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - theta*k(i,j)*dt*xsq*side_1_temp(t)-k(i,j)*theta*dt*ysq;
                if i==1
                    A((i-1)*a+j,(i)*a+j) = k(i,j)*theta*dt*ysq;
                    B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(sin(((t-1)*dt)+0)+phi((i-1)*a+j+1)) + ysq*(0+phi(i*a+j))-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - theta*k(i,j)*dt*xsq*sin((t*dt)+0);
                elseif i==b
                    A((i-1)*a+j,(i-2)*a+j) = k(i,j)*theta*dt*ysq;
                    B((i-1)*a+j,1) = (theta-1)*k(i,j)*dt*(xsq*(sin(((t-1)*dt)+0)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+0)-2*(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - theta*k(i,j)*dt*xsq*sin((t*dt)+0);
                end
            end
        end
    end
    new_phi = A\B;
    phi=new_phi;
end