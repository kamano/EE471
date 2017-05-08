width=20; %(20cm)
height=10; %(10cm)

prompt = 'Enter the number of divisions, N  '; %request user input
N = input(prompt)
h=width/N;
M=floor(height/h);

prompt1 = 'Enter top potential  '; %request user input
top = input(prompt1)

prompt2 = 'Enter bottom potential  '; %request user input
bottom = input(prompt2)

prompt3 = 'Enter left potential   '; %request user input
left = input(prompt3)

prompt4 = 'Enter right potential   '; %request user input
right = input(prompt4)

unks=ones(M,N); %set a matrix for the number of unknowns
top_voltage=top*ones(1,N+1); % set the top potential
bottom_voltage=bottom*ones(1,N+1); %set the bottom potential
left_voltage=left*ones(M,1);   %set left potential
right_voltage=right*ones(M+2,1);  %set right potential

problem=[top_voltage;left_voltage unks;bottom_voltage];
x=[problem right_voltage]; %complete image of what problem looks like

A=-4*eye(N*M,N*M);

b=zeros(N*M,1);

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
        
        if x(i,j+1)==right  %for the side
            b(A_row,1)=-right;

        elseif x(i-1,j)==top
             b(A_row,1)=-top;
        
        
        elseif x(i+1,j)==bottom
            b(A_row,1)=-bottom;
        
        
        elseif x(i,j-1)==left
            b(A_row,1)=-left;
        end
        

        j=j+1;
        A_row=A_row+1;
        A_col=A_col+1;
    end
    i=i+1;
end


A;
figure;
b;

imagesc(A) %display A
colorbar

ANS=A\b;

phi=reshape(ANS,[N,M])';  %reshaping

figure;
imagesc(phi) %display phi






