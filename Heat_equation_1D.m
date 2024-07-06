%Say we have 10 discrete points on rod 1 to 10.
dx = input("Enter value of delta x : ");
dt = input("Enter value of delta t : ");
tao = input("Enter value of theta : ");
c=dt/(dx)*(dx);
% rows 1-10 of phi represent discrete points, columns represent successive time intervals

phi=ones(10,10); % initialising phi matrix

for i=1:10       % assigning values to phi(1,t) at discrete time intervals
    phi(1,i)=sin(((i-1)*dt)+pi/2);
end

A=zeros(9,9);
for i=1:9
    A(i,i)=-2;
    if i<9
    A(i,i+1)=1;
    A(i+1,i)=1;
    end
end
A(9,8)=2; % using ghost node technique

B=zeros(9,10);
B(1,:)=phi(1,:);
I=eye(9);

M=(I-c*(1-tao)*A);
for i=1:9
    N=(I+c*tao*A)*phi(2:10,i) + c*(tao*B(:,i)+(1-tao)*B(:,i+1));
    R=[M N];
    RA=rref(R);
    for j=9:-1:1
        sum=0;
        for k=0:1:8-j
            sum=sum+RA(j,9-k)*phi(10-k,i+1);
        end
        phi(j+1,i+1)= (RA(j,10)-sum)/RA(j,j); % backward substitution
    end
end
 % phi(10,i+1)=RA(9,10)/RA(9,9);
 % phi(2:10,1+i)=phi(2:10,i)+(dt/(dx)*(dx))*(A*phi(2:10,i) + B(:,i));