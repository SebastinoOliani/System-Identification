clear all
clc

%% 
load('C-2017-K2.mat')
[row1,col1]=size(xy);
figure (1)
xlabel('x position')
ylabel('y position')
plot(xy(:,1),xy(:,2),'r.','DisplayName','C-2017-K2')
hold on;
load('P-2011-NO1.mat')
[row2,col2]=size(xy);
plot(xy(:,1),xy(:,2),'b.','DisplayName','R-2011-NO1')
legend
hold off;
%% C-2017-K2
load('C-2017-K2.mat')
Y1=-ones(row1,1);
X1=zeros(row1,5);
for i=1:1:row1
    X1(i,1)=xy(i,1)^2;
    X1(i,2)=xy(i,1)*xy(i,2);
    X1(i,3)=xy(i,2)^2;
    X1(i,4)=xy(i,1);
    X1(i,5)=xy(i,2);
end
theta1=inv(transpose(X1)*X1)*transpose(X1)*Y1;
e1=eccentricity(theta1(1),theta1(2),theta1(3),theta1(4),theta1(5),1);
syms x y
fimplicit(theta1(1)*x^2+theta1(2)*x*y+theta1(3)*y^2+theta1(4)*x+theta1(5)*y+1==0,[-4e-3 1e-3 -7e-4 7e-4])
title('C-2017-K2')
hold on;
plot(xy(:,1),xy(:,2),'rx')
legend('Least square solution', 'Experimental data')
hold off;
%% P-2011-NO1
load('P-2011-NO1.mat')
[row2,col2]=size(xy);
Y2=-ones(row2,1);
X2=zeros(row2,5);
for i=1:1:row2
    X2(i,1)=xy(i,1)^2;
    X2(i,2)=xy(i,1)*xy(i,2);
    X2(i,3)=xy(i,2)^2;
    X2(i,4)=xy(i,1);
    X2(i,5)=xy(i,2);
end
theta2=inv(transpose(X2)*X2)*transpose(X2)*Y2;
e2=eccentricity(theta2(1),theta2(2),theta2(3),theta2(4),theta2(5),1);
syms x y
fimplicit(theta2(1)*x^2+theta2(2)*x*y+theta2(3)*y^2+theta2(4)*x+theta2(5)*y+1==0,[-16e-5 4e-5 -6e-5 6e-5])
title('P-2011-NO1')
hold on;
plot(xy(:,1),xy(:,2),'rx','DisplayName','R-2011-NO1')
legend('Least square solution','Experimental data')
hold off;