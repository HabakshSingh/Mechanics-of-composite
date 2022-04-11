E1=138; %gpa
E2=10.4;G12=4.84;v12=0.25;
S=[1/E1 -v12/E1 0; -v12/E1 1/E2 0;0 0 1/G12];
Q=inv(S);%for 0 degree
a=4; %number of plates
thet=[-45 45 -45 45];   %change it as per the given question
for f=1:a
    theta=thet(f);
    m=cosd(theta);n=sind(theta);
    Tsigma=[m^2 n^2 2*m*n;n^2 m^2 -2*m*n;-m*n m*n (m^2-n^2)];
    Tepsilon=[m^2 n^2 m*n;n^2 m^2 -m*n;-2*m*n 2*m*n (m^2-n^2)];
    Qtheta{f}=inv(Tsigma)*Q*Tepsilon;
end
%now we will do calculation for A,B,C and D matrix.
t=0.25;% thickness of each plate
z=t*[-2 -1 0 1 2];  %to be changed as given in the question
A=zeros(3);
B=zeros(3);
D=zeros(3);
for i=1:3
    for j=1:3        
        for k=1:a
            A(i,j)=A(i,j)+(Qtheta{k}(i,j)*(z(k+1)-z(k)));
            B(i,j)=B(i,j)+Qtheta{k}(i,j)*((z(k+1))^2-(z(k))^2)/2;
            D(i,j)=D(i,j)+Qtheta{k}(i,j)*((z(k+1))^3-(z(k))^3)/3;
        end
    end
end
E=[A B;B D];  %full stiffness bending
%now some external load is applied and we are asked to find the stress.
%then [e0;k]=inv(E)*[applied load]
strain0=inv(E)*[0;0;0;20;0;0]; %for a moment Mx .where stain0=[ex0 ey0 exy0 kx ky kz]
syms z
strain=strain0(1:3,1)+z*strain0(4:6,1);  %strain=[ex;ey;exy]
%now we will calculate stress in each plate
for g=1:a
    stress{g}=((Qtheta{g}*strain));
end