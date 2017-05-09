syms phi_e1_1 phi_e1_2 phi_e1_3 x y

x1=0;
x2=1;
x3=0;
x4=0.3;

y1=0;
y2=0;
y3=1;
y4=0.25;

%{
for element 1
phi_e1_1=a+b*0+c*0;
phi_e1_2=a+b*0.3+c*0.25;
phi_e1_3=a+b*0+c*1;
%}

P_1=[phi_e1_1;phi_e1_2;phi_e1_3];
B_1=[1 0 0;1 0.3 0.25; 1 0 1];
coeff_1=inv(B_1)*P_1;
phi_one=coeff_1(1)+coeff_1(2)*x+coeff_1(3)*y;
A_one=0.3/2;
alpha_one=coeffs(phi_one,flipud(P_1));

for n=1:3
    grad_one(n,1:2)=gradient(alpha_one(n),[x,y]);
end

for m=1:3
    for n=1:3
        c_one(m,n)=A_one*(grad_one(m,1)*grad_one(n,1)+grad_one(m,2)*grad_one(n,2));
    end
end

syms phi_e2_1 phi_e2_2 phi_e2_3
%{
for element 2
phi_e2_1=a+b*0+c*0;
phi_e2_2=a+b*1+c*0;
phi_e2_3=a+b*0.3+c*0.25;
%}

P_2=[phi_e2_1;phi_e2_2;phi_e2_3];
B_2=[1 0 0;1 1 0;1 0.3 0.25];
coeff_2=inv(B_2)*P_2;
phi_two=coeff_2(1)+coeff_2(2)*x+coeff_2(3)*y;
A_two=((0.25)+0.3*(0))/2;
alpha_two=coeffs(phi_two,flipud(P_2));

for n=1:3
    grad_two(n,1:2)=gradient(alpha_two(n),[x,y]);
end

for m=1:3
    for n=1:3
        c_two(m,n)=A_two*(grad_two(m,1)*grad_two(n,1)+grad_two(m,2)*grad_two(n,2));
    end
end

syms phi_e3_1 phi_e3_2 phi_e3_3

%{
for element 3
phi_e3_1=a+b*0.3+c*0.25;
phi_e3_2=a+b*1+c*0;
phi_e3_3=a+b*0+c*1;
%}

P_3=[phi_e3_1;phi_e3_2;phi_e3_3];
B_3=[1 0.3 0.25;1 1 0;1 0 1];
coeff_3=inv(B_3)*P_3;
phi_three=coeff_3(1)+coeff_3(2)*x+coeff_3(3)*y;
A_3=(0.3*(-1)+(1-0.25))/2;
alpha_three=coeffs(phi_three,flipud(P_3));

for n=1:3
    grad_three(n,1:2)=gradient(alpha_three(n),[x,y]);
end

for m=1:3
    for n=1:3
        c_three(m,n)=A_3*(grad_three(m,1)*grad_three(n,1)+grad_three(m,2)*grad_three(n,2));
    end
end

C=[c_one(1,1)+c_two(1,1) c_three(1,2) c_one(1,3) c_one(1,2)+c_two(1,3);
    c_three(1,2) c_three(2,2)+c_two(2,2) c_three(2,3) c_two(2,3)+c_three(1,2);
    c_one(1,3) c_three(2,3) c_one(3,3)+c_three(3,3) c_three(1,3)+c_one(2,3);
    c_one(1,2)+c_two(1,3) c_two(2,3)+c_three(1,2) c_three(1,3)+c_one(2,3) c_one(2,2)+c_two(3,3)+c_three(1,1)];

syms phi4

phi=[0 50 50 phi4];
L=phi*C*transpose(phi);

phi4=solve(diff(L,phi4)==0,phi4);
phi4=double(phi4)
phi=[50 0 0 0; 0 100 0 0; 0 phi4 100 0; 0 0 0 50];
imagesc(phi)
colorbar
