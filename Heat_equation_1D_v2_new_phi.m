%Say we have 10 discrete points on rod 1 to 10.
dx = input("Enter value of delta x : ");
dt = input("Enter value of delta t : ");
tao= input("Enter value of theta : ");
c=dt/(dx)*(dx);

phi=ones(10,1);         %setting values of phi for points 1 to 10 to be 1.
new_phi=ones(10,1);     %initialising new_phi column matrix

point_1_temp=ones(10,1);
for i=1:10
    point_1_temp(i)=sin(((i-1)*dt)+pi/2);
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
B(1,:)=point_1_temp;
I=eye(9);

M=(I-c*(1-tao)*A);
for i=1:9
    new_phi(1)=point_1_temp(i+1);
    N=(I+c*tao*A)*phi(2:10) + c*(tao*B(:,i)+(1-tao)*B(:,i+1));
    R=[M N];
    RA=rref(R);
    for j=9:-1:1
        sum=0;
        for k=0:1:8-j
            sum=sum+RA(j,9-k)*new_phi(10-k);
        end
        new_phi(j+1)= (RA(j,10)-sum)/RA(j,j); % backward substitution
    end
    phi=new_phi;
end