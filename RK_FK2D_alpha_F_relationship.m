clear;
clc;

ht  = 10;       %行数
wth = 10;       %列数
sel = [ht wth]; %方便索引，使用sub2ind(sel,x,y)
n   = ht*wth;   %方便索引，总粒子数
t0  = [0 10];   %求解区间
m0  = zeros(n*4,1);    %初始位置与速度. 分配规则为：按列依次对应，分成四大块
%分别为x轴位置、y轴位置、x速度、y速度.例如：
%   1   3;  5   7; 9    11; 13  15
%   2   4;  6   8; 10   12; 14  16
%分别为粒子对应x位置、y位置、x速度、y速度四个矩阵
A_e=1;%外力大小

j=40;%把0.5*pi细分成40份
F_s=zeros(j+1,1);%每个角度对应的静摩擦
F_c=zeros(j+1,1);
v0=1e-1;
beta0=1e-2;%以上两参量用来判断何时静止以及取向角度相同
k=200;%力的均分份数
T=0;

for alpha=0:(1/j)*0.5*pi:0.5*pi
    fprintf('The process has reached %.3f percent\n',100*alpha/(0.5*pi));  
   for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
   end
   v_x=normrnd(0, sqrt(0.5*T), sel);
   v_y=normrnd(0, sqrt(0.5*T), sel);
   for ii = 1:wth %该循环设置粒子初始速度
     for tt = 1:ht
       m0(2*n+sub2ind(sel,ii,tt))=v_x(ii,tt);
       m0(3*n+sub2ind(sel,ii,tt))=v_y(ii,tt);
     end
   end
   %计算平均速度
   v=zeros(k,1);
   beta=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   parfor ii=1:k
    %tic
    A_e=ii*0.01;%力可以再变得细一些
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，再对时间求平均
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v(ii)=mean(v_t(round(size(v_t)/2):size(v_t)),'all');
    beta(ii)=beta_t(size(beta_t,1));
    if (abs(beta(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v(ii))>v0) 
       F_s1(ii)=A_e;
    end
    %fprintf('Elapsed time: %.3f seconds\n', toc);
   end
   jj=round(alpha/((1/j)*0.5*pi))+1;
   for ii=1:k
       if (F_s1(ii)>0)&&(F_s(jj)==0)%静摩擦临界点
           F_s(jj)=F_s1(ii);
           kk=ii;
       end
   end
   for ii=kk:k%从静摩擦临界点往后找
       if  (F_c1(ii)>0)&&(F_c(jj)==0)
           F_c(jj)=F_c1(ii);
       end
   end
end
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
plot(linspace(0,pi/2,j+1),F_c,'-^r','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(linspace(0,pi/2,j+1),F_s,'-.vb','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('$\alpha$','Interpreter', 'latex');
ylabel('$F$','Interpreter', 'latex');
legend('$F_{c}$','$F_{s}$')
function dmdt = odefun(t,m,n,sel,ht,wth,A_e,alpha,T)
g   = 1;        %弹性系数
mass= 1;      %质量
a0=1;         %原长
dmdt = zeros(n*4,1);
F_e=A_e*[cos(alpha),sin(alpha)];
gamma=0;%阻尼系数（欠阻尼）
kapa=0;
omega=sqrt(T)*0.2;
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
        dmdt(2*n+ii) = (fx*(1+kapa*cos(omega*t)) + F_i(1)+F_e(1)-gamma*m(2*n+ii)) / mass;
        dmdt(3*n+ii) = (fy*(1+kapa*cos(omega*t)) + F_i(2)+F_e(2)-gamma*m(3*n+ii)) / mass;
end
end

function [fx,fy] = basepotential(x,y)
A   = 0.1;     %外场强度
as  = 0.5;        %外场周期
fx  = -A*sin(x*2*pi/as);%*cos(y*2*pi/as);
fy  = -A*sin(y*2*pi/as);%cos(x*2*pi/as);
end