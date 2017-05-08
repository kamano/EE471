clear
frequency=1000; %(Hz)
lambda=(3*10^8)/frequency;
a=0.3*lambda; %outer radius (cm)
b=0.25*lambda; %inner radius (cm)
N=20; %number of divisions
theta0=360/N; %angle of each division
Er=4; %dielectric constant 

k=2*pi/lambda; %k (or beta from the wave equation)

r=sqrt((a^2-b^2)/N); %radius of equivalent area circle
c=(a+b)/2; %distance from the center to the middle of a and b
theta=theta0/2;
f=1;

while f<N+1
    if f==1
        theta=0;
    end
    if f==2
        theta=theta0/2;
    end
    
xi(1,f)=c*cosd(theta); %calculating x's
yi(1,f)=c*sind(theta); %calculating y's
thetab(f,1)=theta; %keeping track of thetas
theta=theta+theta0; 
Ei(1,f)=exp(-1i*k*xi(1,f)); %initial electric field

f=f+1;
end


A=ones(N,N);
b=ones(N,1);
m=1;
while m<(N+1)
    n=1;
    while n<(N+1)
        
        if m==n %when m==n, (diagonals of matrix)
        A(m,n)=1+((Er-1)*(1i/2)*(pi*k*r*besselh(1,2,k*r)-(2*1i)));
        
        elseif m~=n %when m=!n for everything else
        ro=sqrt((xi(1,m)-xi(1,n))^2+(yi(1,m)-yi(1,n))^2); %the ro between the different sections
        A(m,n)=(Er-1)*(1i*pi*k*r/2)*besselj(1,k*r)*besselh(0,2,k*ro); 
        end
        
        n=n+1;
    end
    b(m,1)=Ei(1,m); %the Ei that goes along with that distance
    m=m+1;
end

ans=inv(A)*b; %getting the answer
W=ans(1:(N/2),:); %taking only the first half
RAWR=thetab(1:(N/2),:); %from 0 to 180 theta

plot(RAWR,abs(W)) %plotting 
