clc
clear all
%solving only for N terms
length=1; %1 meter, length of the plate
width=1; %1 meter, width of the plate
n=5; %the number of divisions

X=length/n; %length of division in x direction
Y=width/n; %length of division in y direction

area=(X/2)^2*pi; %area for equivalent circle
epsilon=8.854*10^(-12); %epsilon

N=n*n;%for a plate

b=epsilon*ones(N,1); %top plate is 1

A=ones(N,N); %creating A matrix

yi=1;
  
d=0; %0 and then 1 meter, distance between the plates

x=X/2:X:length; %the spacings for x
y=Y/2:Y:width;
%creating a matrix for x's
counter=1;
H=x';
while counter<n
    H=[H; x'];
    counter=counter+1;
end
%creating a matrix for y's
K=y;
counter=1;
while counter<n
    K=[K;y];
    counter=counter+1;
end

Testing_point=1;

while Testing_point<N+1
     point_2=1;
        while point_2<N+1
            
            
                R=sqrt((H(Testing_point)-H(point_2))^2+(K(Testing_point)-K(point_2))^2);
                
                    if R==0
                        A(Testing_point,point_2)=0.282*sqrt(area); %diagonal elements
                    else
                        A(Testing_point,point_2)=area/(4*pi*abs(R)); %non-diaginal elements
                    end
                point_2=point_2+1;     
            
        end            
     Testing_point=Testing_point+1;  
     
end
result=A\b;
Top=reshape(result,[n,n]);
Bottom=-Top;
subplot(2,1,1)
surf(Top)
subplot(2,1,2)
surf(Bottom)

figure;

subplot(2,1,1)
imagesc(Top)
colorbar
subplot(2,1,2)
imagesc(Bottom)
colorbar


S=sum(result);
v=2;
cap=S/v;