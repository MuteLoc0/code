clear;
clc;
hold off
%-----输出图像预备-----%
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

ht  = 10;       %行数
wth = 10;       %列数
sel = [ht wth]; %方便索引，使用sub2ind(sel,x,y)
n   = ht*wth;   %方便索引，总粒子数
t0  = [0 10];   %求解区间
m0  = zeros(n*4,1);    %初始位置与速度. 分配规则为：按列依次对应，分成四大块
%分别为x轴位置、y轴位置、x速度、y速度.例如：
%   1   3;  5   7; 9    11; 13  14
%   2   4;  6   8; 10   12; 14  16
%分别为粒子对应x位置、y位置、x速度、y速度四个矩阵
jj_max=1;%把0.5*pi细分成40份
v0=0.2;
beta0=0.05;%以上两参量用来判断何时静止以及取向角度相同
k=100;%力的均分份数
F_s=zeros(jj_max,1);
F_c=zeros(jj_max,1);
kapa=0.5;
T=0;
%生成初始速度
v_x=randn(n,1);
v_x=(v_x-mean(v_x))/std(v_x)*sqrt(0.5*T);
v_y=randn(n,1);
v_y=(v_y-mean(v_y))/std(v_y)*sqrt(0.5*T);
alpha_l=linspace(0,0.5*pi,jj_max);%生成角度梯度
%%
for jj=1:jj_max
    alpha=0;
for ii = 1:wth %该循环设置粒子初始速度
      m0(2*n+ii)=v_x(ii);
      m0(3*n+ii)=v_y(ii);
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v1=zeros(k,1);
   beta1=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.01,1.5,k);%设置力梯度
   parfor ii=1:k
    tic
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    T_t=(mean(m(:,2*n+1:3*n).^2,2)+mean(m(:,3*n+1:4*n).^2,2))-(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2);
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v1(ii)=mean(v_t(round(size(v_t)/2):size(v_t)),'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta1(ii)=beta_t(size(beta_t,1));
    if (abs(beta1(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v1(ii))>v0) 
        F_s1(ii)=A_e;
    end
    toc
   end
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk1=ii;
       end
   end
   for ii=kk1:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
%%
F_s=zeros(jj_max,1);
F_c=zeros(jj_max,1);
for jj=1:jj_max
    alpha=pi/3;
for ii = 1:wth %该循环设置粒子初始速度
      m0(2*n+ii)=v_x(ii);
      m0(3*n+ii)=v_y(ii);
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v2=zeros(k,1);
   beta2=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.01,1.5,k);%设置力梯度
   parfor ii=1:k
    tic
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    T_t=(mean(m(:,2*n+1:3*n).^2,2)+mean(m(:,3*n+1:4*n).^2,2))-(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2);
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v2(ii)=mean(v_t(round(size(v_t)/2):size(v_t)),'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta2(ii)=beta_t(size(beta_t,1));
    if (abs(beta2(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v2(ii))>v0) 
        F_s1(ii)=A_e;
    end
    toc
   end
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk2=ii;
       end
   end
   for ii=kk2:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
%%
F_s=zeros(jj_max,1);
F_c=zeros(jj_max,1);
for jj=1:jj_max
    alpha=pi/6;
for ii = 1:wth %该循环设置粒子初始速度
      m0(2*n+ii)=v_x(ii);
      m0(3*n+ii)=v_y(ii);
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v3=zeros(k,1);
   beta3=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.01,1.5,k);%设置力梯度
   parfor ii=1:k
    tic
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v1_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    T_t=(mean(m(:,2*n+1:3*n).^2,2)+mean(m(:,3*n+1:4*n).^2,2))-(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2);
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v3(ii)=mean(v1_t(round(size(v1_t)/2):size(v1_t)),'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta3(ii)=beta_t(size(beta_t,1));
    if (abs(beta3(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v3(ii))>v0) 
        F_s1(ii)=A_e;
    end
    toc
   end
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk3=ii;
       end
   end
   for ii=kk3:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
%%

F_s=zeros(jj_max,1);
F_c=zeros(jj_max,1);
for jj=1:jj_max
    alpha=pi/2;
for ii = 1:wth %该循环设置粒子初始速度
      m0(2*n+ii)=v_x(ii);
      m0(3*n+ii)=v_y(ii);
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v4=zeros(k,1);
   beta4=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.01,1.5,k);%设置力梯度
   parfor ii=1:k
    tic
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    T_t=(mean(m(:,2*n+1:3*n).^2,2)+mean(m(:,3*n+1:4*n).^2,2))-(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2);
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v4(ii)=mean(v_t(round(size(v_t)/2):size(v_t)),'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta4(ii)=beta_t(size(beta_t,1));
    if (abs(beta4(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v4(ii))>v0) 
        F_s1(ii)=A_e;
    end
    toc
   end
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk4=ii;
       end
   end
   for ii=kk4:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
%%
F_s=zeros(jj_max,1);
F_c=zeros(jj_max,1);
for jj=1:jj_max
    alpha=pi/4;
for ii = 1:wth %该循环设置粒子初始速度
      m0(2*n+ii)=v_x(ii);
      m0(3*n+ii)=v_y(ii);
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v5=zeros(k,1);
   beta5=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.01,1.5,k);%设置力梯度
   parfor ii=1:k
    tic
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    T_t=(mean(m(:,2*n+1:3*n).^2,2)+mean(m(:,3*n+1:4*n).^2,2))-(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2);
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v5(ii)=mean(v_t(round(size(v_t)/2):size(v_t)),'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta5(ii)=beta_t(size(beta_t,1));
    if (abs(beta5(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v5(ii))>v0) 
        F_s1(ii)=A_e;
    end
    toc
   end
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk5=ii;
       end
   end
   for ii=kk5:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
%%
beta1=smooth(beta1);
beta2=smooth(beta2);
beta3=smooth(beta3);
beta4=smooth(beta4);
beta5=smooth(beta5);
plot(A_l(kk1:size(beta1,1)),beta1(kk1:size(beta1,1)),'-^r','MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r');
ylim([-0.1,1.6])
xlabel('$F$','Interpreter', 'latex');
ylabel('$\beta$','Interpreter', 'latex');
hold on
plot(A_l(kk2+4:size(beta2,1)),beta2(kk2+4:size(beta2,1)),'-vk','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(A_l(kk3+4:size(beta3,1)),beta3(kk3+4:size(beta3,1)),'-ob','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot(A_l(kk4:size(beta4,1)),beta4(kk4:size(beta4,1)),'-squareg','MarkerSize',3,'MarkerEdgeColor','g','MarkerFaceColor','g');
plot(A_l(kk5:size(beta5,1)),beta5(kk5:size(beta5,1)),'-dc','MarkerSize',3,'MarkerEdgeColor','c','MarkerFaceColor','c');
legend('$\alpha=0$','$\alpha=\frac{\pi}{3}$','$\alpha=\frac{\pi}{6}$','$\alpha=\frac{\pi}{2}$','$\alpha=\frac{\pi}{4}$')
%%
% v11=[v1(1:16);smooth(v1(17:size(v1,1)))];
% 
% v31=[v3(1:18);smooth(v3(19:size(v3,1)))];
% 
% v51=[v5(1:21);smooth(v5(22:size(v5,1)))];
% plot(A_l,v11,'-^r','MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r');
% xlabel('$F$','Interpreter', 'latex');
% ylabel('$v$','Interpreter', 'latex');
% hold on
% % plot(A_l,v2,'-vk','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
% plot(A_l,v31,'-ob','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b');
% % plot(A_l,v4,'-squareg','MarkerSize',3,'MarkerEdgeColor','g','MarkerFaceColor','g');
% plot(A_l,v51,'-dk','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
% legend('$\alpha=0$','$\alpha=\frac{\pi}{6}$','$\alpha=\frac{\pi}{4}$')

function dmdt = odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa)
g   = 1;        %弹性系数
mass= 1;      %质量
a0=1;         %原长
dmdt = zeros(n*4,1);
F_e=A_e*[cos(alpha),sin(alpha)];
gamma=0;%阻尼系数（嘎掉）
omg=sqrt(T);
for ii = 1:n
    dmdt(ii)    = m(2*n+ii);
    dmdt(n+ii)  = m(3*n+ii);
end
for ii = 1:n
    [y,x] = ind2sub(sel,ii);
    if x==1&&y==1
        F_i=g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==1&&y==ht
        F_i=g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==ht
        F_i=g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==1
        F_i=g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif x==1
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==wth
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif y==1
        F_i = g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif y==ht
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    else
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    end
        [fx,fy] = basepotential(m(ii),m(n+ii));
        dmdt(2*n+ii) = (fx*(1+kapa*cos(omg*t)) + F_i(1)+F_e(1)-gamma*m(2*n+ii)) / mass;
        dmdt(3*n+ii) = (fy*(1+kapa*cos(omg*t)) + F_i(2)+F_e(2)-gamma*m(3*n+ii)) / mass;
end
end

function [fx,fy] = basepotential(x,y)
A   = 0.35;     %外场强度
as  = 0.5;        %外场周期
fx  = -A*sin(x*2*pi/as);%*cos(y*2*pi/as);
fy  = -A*sin(y*2*pi/as);%cos(x*2*pi/as);
end