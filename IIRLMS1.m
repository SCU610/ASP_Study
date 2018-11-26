%% This is a test .m file

clear;
clc;
%IIR-LMS
k=800;
wk=zeros(3,k);
wk(:,1)=[0 0 0]';% ��ʼȨֵ,���е�һ�б�ʾȨֵa0���ڶ��б�ʾȨֵb1���ڶ��б�ʾȨֵb2
alpha0=zeros(1,k);  %�˴���ʼ��alpha0 beita1��beita2
beita1=zeros(1,k);
beita2=zeros(1,k);
xk=sqrt(0.5)*randn(1,k);%���������
M=diag([0.05 0.005 0.0025]);%M��ֵ
dk=zeros(1,k);%����dk������
dk(1)=xk(1);  %�ȼ���dk��ǰ����ֵ
dk(2)=xk(2)+dk(1)*1.2;
yk=zeros(1,k);
yk(1)=wk(1,1)*xk(1);%�ȼ���yk��ǰ����ֵ
yk(2)=wk(1,1)*xk(2)+yk(1)*wk(2,2);
Uk=zeros(3,k);
Uk(:,1)=[xk(1) 0 0]';%�ȼ���U��ǰ����ֵ
Uk(:,2)=[xk(2) yk(1) 0]';
err=zeros(1,k);err(1)=dk(1)-yk(1);err(2)=dk(2)-yk(2);
for i=3:k-1
Uk(:,i)=[xk(i) yk(i-1) yk(i-2)]';
yk(i)=wk(:,i)'*Uk(:,i);%����ź�
dk(i)=xk(i)+1.2*dk(i-1)-0.6*dk(i-2);%�����ź�
err(i)=dk(i)-yk(i);%���
alpha0(i)=xk(i)+wk(2,i)*alpha0(i-1)+wk(3,i)*alpha0(i-2);
beita1(i)=yk(i-1)+wk(2,i)*beita1(i-1)+wk(3,i)*beita1(i-2);
beita2(i)=yk(i-2)+wk(2,i)*beita2(i-1)+wk(3,i)*beita2(i-2);
gra=-2*err(i)*[alpha0(i) beita1(i) beita2(i)]';
wk(:,i+1)=wk(:,i)-M*gra;
end

%���۲���ͼ
ite=1:k;
figure
plot(ite,err.^2);
title('lms��ѧϰ����');
xlabel('������');
ylabel('�����');
hold on 
%lmsʧ��ͼ
b1=wk(2,:);b2=wk(3,:);
z1=0.5*(b1+sqrt(b1.^2+4*b2));%����z1
z2=0.5*(b1-sqrt(b1.^2+4*b2));%����z2
%����z1��������
R1=(1./(1-b1.*z1-b2.*z1.^2)-2./(0.6*z1.^2-1.2*z1+1)).*z1./(z1-z2);
%����z2��������
R2=(1./(1-b1.*z2-b2.*z2.^2)-2./(0.6*z2.^2-1.2*z2+1)).*z2./(z2-z1);
MSEzhen1=25/7+R1+R2;
MSEzhen1=real(MSEzhen1);
fenmu=1/min(MSEzhen1);
shitiao=MSEzhen1.*fenmu-min(MSEzhen1)*fenmu;
figure
plot(ite,shitiao);
title('lmsʧ��');
xlabel('������');ylabel('ʧ��');
hold on;


figure
plot3(wk(2,:),wk(3,:),MSEzhen1,'r');
hold on 


%IIR-SER
k=650;
lamda=0.51;%��av
a=0.93;
q0=1;
%������ʼֵ
Q=q0*eye(3);%Q0����������ʼֵ
wk1=zeros(3,k);
wk1(:,1)=[0 0 0]';% ��ʼȨֵ,���е�һ�б�ʾȨֵa0���ڶ��б�ʾȨֵb1���ڶ��б�ʾȨֵb2
k=800;
for i=3:k-1
Uk(:,i)=[xk(i) yk(i-1) yk(i-2)]';
yk(i)=wk1(:,i)'*Uk(:,i);%����ź�
dk(i)=xk(i)+1.2*dk(i-1)-0.6*dk(i-2);%�����ź�
err(i)=dk(i)-yk(i);%���
alpha0(i)=xk(i)+wk1(2,i)*alpha0(i-1)+wk1(3,i)*alpha0(i-2);
beita1(i)=yk(i-1)+wk1(2,i)*beita1(i-1)+wk1(3,i)*beita1(i-2);
beita2(i)=yk(i-2)+wk1(2,i)*beita2(i-1)+wk1(3,i)*beita2(i-2);
gra=-2*err(i)*[alpha0(i) beita1(i) beita2(i)]';
S=Q*Uk(:,i);
v=a+Uk(:,i)'*S;
Q=(1/a)*(Q-S*S'/v);
wk1(:,i+1)=wk1(:,i)-M*lamda*(1-a^(i+1))*Q*gra/(1-a);
end


%%�ȸ�ƽ��ͼ
[b1,b2]=meshgrid(-2:0.02:2,-1:0.01:1);
z1=0.5*(b1+sqrt(b1.^2+4*b2));%����z1
z2=0.5*(b1-sqrt(b1.^2+4*b2));%����z2
%����z1��������
R1=(1./(1-b1.*z1-b2.*z1.^2)-2./(0.6*z1.^2-1.2*z1+1)).*z1./(z1-z2);
%����z2��������
R2=(1./(1-b1.*z2-b2.*z2.^2)-2./(0.6*z2.^2-1.2*z2+1)).*z2./(z2-z1);
MSE=25/7+R1+R2;
MSE=real(MSE);
figure,contour(b1,b2,MSE,[0.1  0.3 0.5 0.7 0.9 0.99]*5);%�ȸ���
axis([-2 2 -1 1]);
hold on;
plot(wk(2,:),wk(3,:),'r');
hold on
title('IIR LMS/IIR SER')
xlabel('Ȩ��b1')
ylabel('Ȩ��b2')
hold on;
plot(1.2,-0.6,'*');                    %��ע���Ȩֵ��λ��
hold on;plot(wk1(2,:),wk1(3,:),'b');
legend('','IIR-LMS', '','IIR-SER');

%��άͼ��
figure
plot3(1.2,-0.6,1/fenmu,'*');     %��ѵ� 
title('IIR��άͼʾ');
hold on;
surf(b1,b2,MSE);
hold on;
plot3(wk(2,:),wk(3,:),MSEzhen1,'r');
hold on


