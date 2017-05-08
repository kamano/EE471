N=70; %the number of divisions

c=3*10^8; %speed of light
frequency= 2*10^9; %1 GHz
wavelength= 2*pi*c/frequency;
a=0.001; %radius

prompt = 'Length of Antenna in terms of lambda  '; %request user input
len = input(prompt);

total_length=len*wavelength; %total length of the antenna
length=total_length/2; %length of one half of thing
l=total_length/N; %length of a segment

z=-1*length-l/2:l:length-l/2;
z(N+1)=length-l;

b=4*pi*ones(N+1,1);
A=ones(N+1,N+1);
B0=(2*pi)/wavelength;
Vg=1; 
eta_naught=120*pi;

i=1;
%need N+1 unks and equations
while i<(N+2)
    j=1;
    
    while j<(N+2)
       R=sqrt(a^2+(z(i)-z(j))^2);
       
       A(i,j)=A(i,j)*l*(exp(-1i*B0*R)/R);
      
       if j==N+1 %for the last colum
           A(i,j)=cos(B0*z(i)); 
       end
       
       j=j+1;
    end
    b(i,1)=b(i,1)*-1i*Vg*sin(B0*abs(z(i)))/(2*eta_naught);
    i=i+1;
end

result=A\b;
current=result(1:N);
current_middle=result(round(N/2))
impedance=Vg/current_middle

cur=abs(current);

phi=reshape(abs(result),[1,N+1]);
phi(:, N+1)=[];
phi(1,1)=0;

z(:, N+1)=[];
plot(phi,z)
