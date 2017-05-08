clear
length= 1; %1meter
N= 10; %the number of divisions
a=0.001; %1mm 
l=length/N;

epsilon=8.854*10^(-12);

b=4*pi*epsilon*ones(N,1);
A=ones(N,N);

i=1;
while i<(N+1)
    j=1;
    while j<(N+1)
       if (i*l-j*l)~=0
           A(i,j)=A(i,j)*2*pi*a*l/(abs(i*l-j*l));
       else 
           A(i,j)=A(i,j)*4*pi*log(l/a);
       end
       j=j+1;
    end
    i=i+1;
end

result=A\b;

phi=reshape(result,[1,N])
phi=[phi phi(:,end)]
x=linspace(0,length,N+1);
stairs(x,phi)

