%The unshielded microstrip
clc
clear all
width=0.12;
metal_width=0.01;
potential_on_metal=100; %in volts
height_metal=0.002;
s=0.01; %height of the dielectric
height_freespace=0.1;
e1=1; %dielectric for freespace
e2=4; %dielectric for bottom

N=3;

h=metal_width/N;
M1=floor(s/h); %number of nodes in dielectric
M2=round(height_freespace/h); %number of nodes freespace
N1=round(width/h);

dielectric=e2*ones(M1,N1);
freespace=e1*ones(M2,N1);

midpoint=floor(N1/2)-floor(N/2);
metal=e1*ones(1,N1);
counter=0;

while counter<N
    metal(1,midpoint+counter)=100;
    counter=counter+1;
end

top_bottom=0*ones(1,N1);
left_right=0*ones(M1+M2+3,1);
problem=[freespace; metal; dielectric];

problem=[top_bottom; problem; top_bottom];

problem=[left_right problem left_right];

A_col=1;
A_row=1;

b=0*ones((M1+M2+1)*N1,1);

i=2;
imagesc(problem)
% with dielectric
while i<M1+M2+3

    j=2;
    while j<(N1+2)
        
        
        if problem(i,j)==e1 % for freespace self
            A(A_row,A_col)=-4*e1;
        end
        
        if problem(i,j)==e2 % for dielectric self
            A(A_row,A_col)=-4*e2;
        end
      
        if problem(i,j+1)==e1 %freespace right of x
            A(A_row,A_col+1)=e1;
        end
        
        if problem(i,j+1)==e2 %dielectric right of x
            A(A_row,A_col+1)=e2;
        end
        
        if problem(i,j-1)==e1 %freespace left of x
            A(A_row,A_col-1)=e1;
        end
        
        if problem(i,j-1)==e2 %dielectric left of x
            A(A_row,A_col-1)=e2;
        end
        
        if problem(i-1,j)==e1 && problem(i+1,j)==e1 %freespace - if not border of freespace and dielectric
            A(A_row,A_col-N1)=e1; %above
            A(A_row,(i-1)*N1+j-1)=e1; %below
        end
        
        if problem(i-1,j)==e2 && problem(i+1,j)==e2 %dielectric - if not border of freespace and dielectric
            A(A_row,A_col-N1)=e2; %above
            A(A_row,(i-1)*N1+j-1)=e2; %below
        end
        
        if problem(i-1,j)==100 || problem(i+1,j)==100
            A(A_row,A_col)=-4*(e1+e2); % for itself
            A(A_row,A_col-N1)=2*e1; %for above
            A(A_row,(i-1)*N1+j-1)=2*e2; %for below
            A(A_row,A_col-1)=e2+e1;%for left
            A(A_row,A_col+1)=e2+e1;%for right
        end
        
        if problem(i-1,j)==e1 && problem(i+1,j)==e2  %if at border of freespace and dielectric
            A(A_row,A_col)=-4*(e1+e2); % for itself
            A(A_row,A_col-N1)=2*e1; %for above
            A(A_row,(i-1)*N1+j-1)=2*e2; %for below
            A(A_row,A_col-1)=e2+e1;%for left
            A(A_row,A_col+1)=e2+e1;%for right
        end
        
        if problem(i,j)==100
            b(j,1)=100;
            A(A_row,A_col)=1;
            A(A_row,A_col-1)=0;
            A(A_row,A_col+1)=0;
            A(A_row,(i-1)*N1+j-1)=0;
            A(A_row,A_col-N1)=0;
        end
        
        j=j+1;
        A_row=A_row+1;
        A_col=A_col+1;
    end
    i=i+1;
end

phi=A\b;

k=M1+M2+1;

result1=reshape(phi,k,N1);

imagesc(problem)
colorbar
figure;
imagesc(result1)
colorbar

%again
%without the dielectric
b2=0*ones((M1+M2+1)*N1,1);
A_row=1;
A_col=1;
i=2;
while i<M1+M2+3

    j=2;
    while j<(N1+2)
        
        
        if problem(i,j)==e1 % for freespace self
            A2(A_row,A_col)=e1;
        end
        
        if problem(i,j)==e2 % for dielectric self
            A2(A_row,A_col)=e1;
        end
      
        if problem(i,j+1)==e1 %freespace right of x
            A2(A_row,A_col+1)=e1;
        end
        
        if problem(i,j+1)==e2 %dielectric right of x
            A2(A_row,A_col+1)=e1;
        end
        
        if problem(i,j-1)==e1 %freespace left of x
            A2(A_row,A_col-1)=e1;
        end
        
        if problem(i,j-1)==e2 %dielectric left of x
            A2(A_row,A_col-1)=e1;
        end
        
        if problem(i-1,j)==e1 && problem(i+1,j)==e1 %freespace - if not border of freespace and dielectric
            A2(A_row,A_col-N1)=e1; %above
            A2(A_row,(i-1)*N1+j-1)=e1; %below
        end
        
        if problem(i-1,j)==100 || problem(i+1,j)==100
            A2(A_row,A_col)=-4*(e1+e2); % for itself
            A2(A_row,A_col-N1)=2*e1; %for above
            A2(A_row,(i-1)*N1+j-1)=2*e2; %for below
            A2(A_row,A_col-1)=e2+e1;%for left
            A2(A_row,A_col+1)=e2+e1;%for right
        end
        
        
        if problem(i,j)==100
            b2(j,1)=100;
            A2(A_row,A_col)=1;
            A2(A_row,A_col-1)=0;
            A2(A_row,A_col+1)=0;
            A2(A_row,(i-1)*N1+j-1)=0;
            A2(A_row,A_col-N1)=0;
        end
        
        j=j+1;
        A_row=A_row+1;
        A_col=A_col+1;
    end
    i=i+1;
end

phi2=A2\b2;

result2=reshape(phi2,k,N1)';

figure;
imagesc(result2)
colorbar


