% Consider a square plate of length l. The temperature of its 3 edges is set to 1 at all times. The y-axis edge
% is varying sinusoidally with time. 
dx = input("Enter value of delta x : ");
dy = input("Enter value of delta y : ");
l = input("Enter value of side length of plate : ");
dt = input("Enter value of delta t : ");
tao = input("Enter value of theta : ");
xsq=1/((dx)^2);
ysq=1/((dy)^2);
a=(l/dx)-1;       %number of points in a row
b=(l/dy)-1;       %number of points in a column
A=zeros(a*b,a*b); %initialising a matrix of the order a*b for finding temperatures of a*b number of points.

phi=ones(a*b,1);         %setting values of phi for all points to be 1.
new_phi=ones(a*b,1);     %initialising new_phi column matrix.

side_1_temp=ones(10,1);
for i=1:10
    side_1_temp(i)=sin(((i-1)*dt)+pi/2);
end
B=zeros(a*b,1);   % initialising the matrix B of A*T(t+1)=B

% for i=2:b-1       % making the matrix A and B of A*T(t+1)=B
%     for j=1:a
%         A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*((1/(dx)^2)+(1/(dy)^2)) - 1;
%         A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2; % upar
%         A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;   % neeche
%         if j==1 || j==a
%             if j==a
%               A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;% left
%               B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx;
%             end % subtracted the constant term from LHS
%             if j==1
%               A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx; %right
%               B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*side_1_temp(t);
%             end % subtracted time dependent term from LHS
%         else
%             A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
%             A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx;
%             B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
%         end
%     end
% end
% % filling matrix A for 1st and last row of heat plate
% for j=2:a-1
%     for i=[1,b]
%         A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*((1/(dx)^2)+(1/(dy)^2)) - 1;
%         A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
%         A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx; 
%         B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2;
%         if i==1
%             A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
%         elseif i==b
%             A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2; 
%         end
%     end   
% end
% for j=[1,a]
%     for i=[1,b]
%         A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*((1/(dx)^2)+(1/(dy)^2)) - 1;
%         if j==a
%             A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
%             B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2- tao*j*dt/dx;
%             if i==1
%                 A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
%             elseif i==b
%                 A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2;
%             end
%         end
%         if j==1
%             A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx;
%            B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*side_1_temp(t)-tao*j*dx*dt/(dy)^2;
%             if i==1
%                 A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
%             elseif i==b
%                 A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2;
%             end
%         end
%     end
% end

for t=1:5
    for i=2:b-1       % making the matrix A and B of A*T(t+1)=B
    for j=1:a
        A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*(xsq+ysq) - 1;
        A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2; % upar
        A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;   % neeche
        if j==1 || j==a
            if j==a
              A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;% left
              B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+1) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx;
            end % subtracted the constant term from LHS
            if j==1
              A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx; %right
              B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j)+sin(((t-1)*dt)+pi/2)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*sin(((t)*dt)+pi/2);
            end % subtracted time dependent term from LHS
        else
            A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
            A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx;
            B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j);
        end
    end
    end
% filling matrix A for 1st and last row of heat plate
for j=2:a-1
    for i=[1,b]
        A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*((1/(dx)^2)+(1/(dy)^2)) - 1;
        A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
        A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx; 
        if i==1
            A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
            B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(1+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2;
        elseif i==b
            A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2; 
            B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+1)+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2;
        end
    end   
end
for j=[1,a]
    for i=[1,b]
        A((i-1)*a+j,(i-1)*a+j) = tao*j*dx*dt*((1/(dx)^2)+(1/(dy)^2)) - 1;
        if j==a
            A((i-1)*a+j,(i-1)*a+j-1) = tao*j*dt/dx;
            % B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2- tao*j*dt/dx;
            if i==1
                A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
                B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+1) + ysq*(1+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2- tao*j*dt/dx;
            elseif i==b
                A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2;
                B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+1) + ysq*(phi((i-2)*a+j)+1)+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j)-tao*j*dx*dt/(dy)^2- tao*j*dt/dx;
            end
        end
        if j==1
           A((i-1)*a+j,(i-1)*a+j+1) = tao*j*dt/dx;
          % B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*side_1_temp(t)-tao*j*dx*dt/(dy)^2;
            if i==1
                A((i-1)*a+j,(i)*a+j) = tao*j*dx*dt/(dy)^2;
                B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(sin(((i-1)*dt)+pi/2)+phi((i-1)*a+j+1)) + ysq*(1+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*sin((t*dt)+pi/2)-tao*j*dx*dt/(dy)^2;
            elseif i==b
                A((i-1)*a+j,(i-2)*a+j) = tao*j*dx*dt/(dy)^2;
                B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(sin(((i-1)*dt)+pi/2)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+1)+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*sin((t*dt)+pi/2)-tao*j*dx*dt/(dy)^2;
            end
        end
    end
end
    % j=1;
    % for i=2:b-1
    %     B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*side_1_temp(t);
    % end
    % for i=[1,b]
    %      B((i-1)*a+j,1) = (tao-1)*j*dx*dt*(xsq*(phi((i-1)*a+j-1)+phi((i-1)*a+j+1)) + ysq*(phi((i-2)*a+j)+phi(i*a+j))+(xsq+ysq)*phi((i-1)*a+j)) - phi((i-1)*a+j) - tao*a*dt/dx*side_1_temp(t)-tao*j*dx*dt/(dy)^2;
    % end
    new_phi = A\B;
    phi=new_phi;
end