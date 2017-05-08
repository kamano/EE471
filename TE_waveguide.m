clc
clear all
width=20; %(20cm)
height=10; %(10cm)

prompt = 'Enter the number of divisions, N  '; %request user input
N = input(prompt)
h=width/N;
M=floor(height/h);


unks=ones(M,N); %set a matrix for the number of unknowns +1 in each dimension
top_voltage=2*ones(1,N); % set the top potential
bottom_voltage=2*ones(1,N); %set the bottom potential
left_voltage=2*ones(M+2,1);  %set left potential
right_voltage=2*ones(M+2,1); %set right potential

problem=[top_voltage;unks;bottom_voltage];
x=[left_voltage problem right_voltage];

A=-4*eye(N*M,N*M);

A_col=1;
A_row=1;

i=2;

while i<M+2

    j=2;
    while j<(N+2)
        i;
        if x(i,j+1)==1 % for adding rightside ones
            A(A_row,A_col+1)=1;
        end
        
        if x(i,j-1)==1      %adding left hand ones
            A(A_row,A_col-1)=1;
        end
        
        if x(i-1,j)==1%if above
            A(A_row,A_col-N)=1;
        
        end
        
        if x(i+1,j)==1  %if below
            A(A_row,(i-1)*N+j-1)=1;
        end
        
        if x(i,j-1)==2 %waveguide left
            A(A_row,A_col+1)=2;
        end
        
        if x(i-1,j)==2 %waveguide below
            A(A_row,(i-1)*N+j-1)=2;
        end
        
        if x(i+1,j)==2  %waveguide above
            A(A_row,A_col-N)=2;
        end
        
        if x(i,j+1)==2 %waveguide right
            A(A_row,A_col-1)=2;
        end
        
        j=j+1;
        A_row=A_row+1;
        A_col=A_col+1;
    end
    i=i+1;
end


[v c]=eig(A); %get the eigen vector and values

c=eig(A); %just eigen values
C=abs(c); %take absolute value of eigen values

[order,Ind]=sort(C,'ascend'); %taking the indices of c
vector=zeros(N*M,N*M);
counter=1;

while counter<(N*M+1)
ind=Ind(counter);
vector(:,counter)=v(:,ind); %need to match eigen values with vectors
counter=counter+1;
end


mu_naught= 4*pi*10^(-9); %for cm
epsilon_naught=8.854*10^(-14); %for cm

n=0;
counter=1;
while n<(N*M)+1
     m=0;
    while m<(N*M)+1
       analytical_cutoff(counter,1)=(1/(2*pi*sqrt(mu_naught*epsilon_naught)))*sqrt((m*pi/width)^2+(n*pi/height)^2);
       m_value(counter,1)=m;
       n_value(counter,1)=n;
       
       counter=counter+1;
       m=m+1;
    end
    n=n+1;
end
n_and_m_values=[m_value n_value];

[ordered_cutoff,index_cutoff]=sort(analytical_cutoff,'ascend');

values=zeros(N*M*N*M,2);

counter=1;
while counter<(N*M*N*M-1)
    index=index_cutoff(counter);
    values(counter,:)=n_and_m_values(index,:); %need to match eigen values with vectors
    counter=counter+1;
end

prompt2 = 'Enter m '; %request user input
m = input(prompt2)

prompt3 = 'Enter n '; %request user input
n = input(prompt3)

i=1;
j=1;
while i<(N*M*N*M)
    if values(i,j)==m && values(i,j+1)==n
        which_eigen_value=i;
    end
    i=i+1;
end

vec=-1*vector(:,which_eigen_value);

phi=reshape(vec,[N,M]);

lamda=order(which_eigen_value); 

experimental_cutoff_answer= sqrt(lamda/(mu_naught*epsilon_naught*h^2))/(2*pi)
analytical_cutoff_answer=(1/(2*pi*sqrt(mu_naught*epsilon_naught)))*sqrt((m*pi/width)^2+(n*pi/height)^2)



%plotting
[Ex,Ey]=gradient(phi,0.1,0.1);
figure;
imagesc(phi)
hold on
quiver(Ex,Ey,'black')
hold off
