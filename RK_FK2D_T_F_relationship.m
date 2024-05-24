clear;
clc;
hold off
%-----输出图像预备-----%
set(0, 'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

ht  = 5;       %行数
wth = 5;       %列数
sel = [ht wth]; %方便索引，使用sub2ind(sel,x,y)
n   = ht*wth;   %方便索引，总粒子数
t0  = [0 10];   %求解区间
m0  = zeros(n*4,1);    %初始位置与速度. 分配规则为：按列依次对应，分成四大块
j=30;%(现在是温度分成10份)
v0=0.1;
beta0=0.05;%以上两参量用来判断何时静止以及取向角度相同
k=101;%力的均分份数
alpha=(0/j)*0.5*pi;
F_ans=zeros(j,5);

for h=1:1:5
kapa=(h-1)*0.05;
F_ss=zeros(j,1,4);
F_cc=zeros(j,1,4);%存每个kappa点四次的结果
F_s=zeros(j,1);
F_c=zeros(j,1);
for jj=1:j
    fprintf('The process has reached %.2f percent\n',(h-1)*20+20*jj/j)
    T=0.1*jj;
%生成初始速度
v_x=sqrt(0.5*T)*normalize(randn(ht));
v_y=sqrt(0.5*T)*normalize(randn(wth));
for ii = 1:wth %该循环设置粒子初始速度
    for tt = 1:ht
      m0(2*n+sub2ind(sel,ii,tt))=v_x(ii,tt);
      m0(3*n+sub2ind(sel,ii,tt))=v_y(ii,tt);
    end
end
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
   %计算平均速度
   v=zeros(k,1);
   beta=zeros(k,1);
   F_c1=zeros(k,1);
   F_s1=zeros(k,1);%用来找每个角度的F_s和F_c时用的中间量
   A_l = linspace(0.1,0.3,k);%设置力梯度
   parfor ii=1:k
    A_e=A_l(ii);
    [t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa),[1,15],m0);
    v_t=(mean(m(:,2*n+1:3*n),2).^2+mean(m(:,3*n+1:4*n),2).^2).^(0.5);%先求每时刻v平均，下面v(ii)再对时间求平均
    beta_t=abs(atan(mean(m(:,3*n+1:4*n),2)./mean(m(:,2*n+1:3*n),2)));
    v(ii)=mean(v_t,'all');
    %beta(ii)=mean(beta_t(round(size(beta_t)/2):size(beta_t)),'all');
    beta(ii)=beta_t(size(beta_t,1));
    if (abs(beta(ii)-alpha)<beta0) 
        F_c1(ii)=A_e;
    end
    if (abs(v(ii))>v0) 
        F_s1(ii)=A_e;
    end
   end
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
F_ans(:,h)=F_s;
end

figure('Color',[1 1 1])
T_set=[0.1:0.1:3]';
plot(T_set,smooth(F_ans(:,1)),'.r-','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
plot(T_set,smooth(F_ans(:,2)),'^b-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
plot(T_set,smooth(F_ans(:,3)),'vc-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
plot(T_set,smooth(F_ans(:,4)),'dg-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
plot(T_set,smooth(F_ans(:,5)),'sm-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('$T$','Interpreter', 'latex');
ylabel('$F_{s}$','Interpreter', 'latex');
legend('$\kappa=0$','$\kappa=0.05$','$\kappa=0.10$','$\kappa=0.15$','$\kappa=0.20$','Interpreter', 'latex');


function dmdt = odefun(t,m,n,sel,ht,wth,A_e,alpha,T,kapa)
g   = 1;        %弹性系数
mass= 0.04;      %质量
a0=1;         %原长
dmdt = zeros(n*4,1);
F_e=A_e*[cos(alpha),sin(alpha)];
gamma=0.2;%阻尼系数（欠阻尼）
omg=3*sqrt(T);
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
A   = 0.3;     %外场强度
as  = 0.5;        %外场周期
fx  = -A*sin(x*2*pi/as);%*cos(y*2*pi/as);
fy  = -A*sin(y*2*pi/as);%cos(x*2*pi/as);
end