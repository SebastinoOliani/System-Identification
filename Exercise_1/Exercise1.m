clc
clear all
load('SysID_Exercise_1.mat');
%data size
dim=size(x);
lenght=dim(1,1);
%noise
mu=0.2;
var=0.1;
%wcomputation
w=zeros(lenght,2);
for i=1:1:lenght
    w(i,1)=x(i);
    w(i,2)=2*x(i)^2-1;
end
%maximum likelihood
theta_ML=inv(transpose(w)*w)*transpose(w)*(y-mu);
%theta as gaussian random variable
mu_theta=[1;0.4];
var_theta=0.01;
%maximum a posteriori
theta_MAP=inv(transpose(w)*w/var+eye(2)/var_theta)*(transpose(w)*(y-mu)/var+mu_theta/var_theta);

%%
%additional measurements
dim_v=size(x_v);
lenght_v=dim_v(1,1);
w_v=zeros(lenght_v,2);

for i=1:1:lenght_v
    w_v(i,1)=x_v(i);
    w_v(i,2)=2*x_v(i)^2-1;
end

%predicted output with maximum likelihood
yML=zeros(lenght_v,1);
for i=1:1:lenght_v
    yML(i,1)=[w_v(i,1) w_v(i,2)]*theta_ML+mu;
end

countML=0;
for i=1:1:lenght_v
    countML=countML+(yML(i)-y_v(i))^2;
end
%error with ML
E_ML=1/lenght_v*countML;

%predicted output with MAP
yMAP=zeros(lenght_v,1);
for i=1:1:lenght_v
    yMAP(i,1)=[w_v(i,1) w_v(i,2)]*theta_MAP+mu;
end

countMAP=0;
for i=1:1:lenght_v
    countMAP=countMAP+(yMAP(i)-y_v(i))^2;
end
%error with MAP
E_MAP=1/lenght_v*countMAP;

if E_MAP<E_ML
    disp('MAP estimate is more accurate');
else disp('ML estimate is more accurate');
end