width=20; %(20cm)
height=10; %(10cm)

prompt = 'Enter the number of divisions, N  '; %request user input
N = input(prompt)
h=width/N;
M=floor(height/h);

prompt1 = 'Enter top potential  '; %request user input
top = input(prompt1);

prompt2 = 'Enter bottom potential  '; %request user input
bottom = input(prompt2);

prompt3 = 'Enter left potential   '; %request user input
left = input(prompt3);

prompt4 = 'Enter right potential   '; %request user input
right = input(prompt4);

unks=ones(M,N); %set a matrix for the number of unknowns

top_voltage=top*ones(1,N+1); % set the top potential
bottom_voltage=bottom*ones(1,N+1); %set the bottom potential
left_voltage=left*ones(M,1);   %set left potential
right_voltage=right*ones(M+2,1);  %set right potential

problem=[top_voltage;left_voltage unks;bottom_voltage];
x=[problem right_voltage]; %complete image of what problem looks like

symmetry_about_y_axis=0; %set as 0 for no symmetry
symmetry_about_x_axis=0; %set as 0 for no symmetry

%determining if there is symmetry 

if mod(N,2)==1 && right==left  %symmetry_about_y_axis
     symmetry_about_y_axis=1
     N=floor(N/2)+1
end

if mod(M,2)==1 && top==bottom %symmetry_about_x_axi
     symmetry_about_x_axis=1  
     M=floor(M/2)+1
end

    if bottom==0 && top==0 && right==0 && left==0
        b=-2*ones(N*M,1);
    else
        b=zeros(N*M,1);
    end



A=-4*eye(N*M,N*M);


A_col=1;
A_row=1;

i=2;


while i<M+2

    j=2;
    while j<(N+2)
        if x(i,j-1)==1  %adding left hand ones
            A(A_row,A_col-1)=1;
        end
        
        if x(i,j+1)==1 % for adding rightside ones

            if symmetry_about_y_axis==1 && j==N+1
                A(A_row,A_col-1)=2;
            else
                A(A_row,A_col+1)=1;
            end
        end
        
        
        if x(i-1,j)==1 %if unknown phi above
            A(A_row,A_col-N)=1;
        end
        
        if x(i+1,j)==1  %if unknown phi below
            if symmetry_about_x_axis==1 && i==M+1
                A(A_row,A_col-N)=2;
            else
                A(A_row,(i-1)*N+j-1)=1;
            end
        end
        
            if x(i,j+1)==right  %determining b matrix entries
                b(A_row,1)=-right;

            elseif x(i-1,j)==top
                 b(A_row,1)=-top;
        
        
            elseif x(i+1,j)==bottom
                 b(A_row,1)=-bottom;
        
        
            elseif x(i,j-1)==left
                 b(A_row,1)=-left;
            end

        if symmetry_about_y_axis==1 && j==N+1
            j=2*N
        end
        j=j+1;
        A_row=A_row+1;
        A_col=A_col+1;
        
        
    end %while j end
    
    i=i+1;
    
end %while i end

if bottom==0 && top ==0 && right==0 && left==0
       b=-2*ones(N*M,1); 
end
        

A; %the complete A matrix
figure;
b; %complete b matrix

imagesc(A) %display A 
colorbar
figure;
imagesc(b)
colorbar

ANS=inv(A)*b; %calculating phi's


phi=reshape(ANS,[N,M])';  %reshaping

figure;
imagesc(phi)

if symmetry_about_x_axis==1 && symmetry_about_y_axis==0%flip up down
Ans2=flipud(phi);

phi1=[phi; Ans2(2:end,:)];
end

if symmetry_about_y_axis==1 %flip left right
Ans1=fliplr(phi);

phi1=[phi Ans1(:,2:end)];
end


if symmetry_about_x_axis==1 && symmetry_about_y_axis==1
    Ans2=flipud(phi1);
    phi1=[phi1; Ans2(2:end,:)];
end

if symmetry_about_x_axis==0 && symmetry_about_y_axis==0
    phi1=phi;
end

figure;
imagesc(phi1)
 
