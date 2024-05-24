clear;
clc;
ht  = 36;       %行数
wth = 36;       %列数
dett= 4;        %计算温度的方块边长
if mod(ht,dett)||mod(wth,dett)
    error('dett ERROR')
end
sel = [ht wth]; %方便索引，使用sub2ind(sel,x,y)
n   = ht*wth;   %方便索引，总粒子数
t0  = [0 10];   %求解区间
T   = 1;%不是温度，是用来调初始温度的一个量
T0=5; %表征基底温度的量，不是基底温度
m0  = zeros(n*4,1);    %初始位置与速度. 分配规则为：按列依次对应，分成四大块
%分别为x轴位置、y轴位置、x速度、y速度.例如：
%   1   3;  5   7; 9    11; 13  15
%   2   4;  6   8; 10   12; 14  16
%分别为粒子对应x位置、y位置、x速度、y速度四个矩阵
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end
v_x=normrnd(0, sqrt(0.5*T), sel);
v_y=normrnd(0, sqrt(0.5*T), sel);
% for ii = 1:wth %该循环设置粒子初始速度正态分布
%     for tt = 1:ht
%       m0(2*n+sub2ind(sel,ii,tt))=v_x(ii,tt);
%       m0(3*n+sub2ind(sel,ii,tt))=v_y(ii,tt);
%     end
% end
y0=round(ht/2)+1/2;
x0=round(wth/2)+1/2;
for xx = 1:wth
    for yy = 1:ht
        m0(2*n+sub2ind(sel,yy,xx))=sqrt(0.5*T)*(x0-abs(xx-x0))/10+sqrt(0.5*T0);
        m0(3*n+sub2ind(sel,yy,xx))=sqrt(0.5*T)*(y0-abs(yy-y0))/10+sqrt(0.5*T0);
    end
end

[t,m]=ode45(@(t,m) odefun(t,m,n,sel,ht,wth),[0,30],m0);

EE = zeros(ht,wth,size(t,1));%三维数组！
for ii = 1:ht %该循环计算各点能量
    for tt = 1:wth
      EE(ii,tt,:) = m(:,2*n+sub2ind(sel,ii,tt)).^2+m(:,3*n+sub2ind(sel,ii,tt)).^2;
    end
end

Temp = zeros(ht/dett,wth/dett,size(t,1));
for ii = 1:ht/dett  %该循环存储各区温度
    for tt = 1:wth/dett
        Temp(ii,tt,:) = squeeze(sum(sum(EE(dett*(ii-1)+1:ii*dett,dett*(tt-1)+1:tt*dett,:),1),2)./dett^2);
    end
end

for tt = 1:size(t,1)
figure(1)
h = heatmap(Temp(:,:,tt),'ColorLimits', [5, 15]);
colormap(h, 'hot');
h.CellLabelColor = 'none';
pause(0.01)
end

function dmdt = odefun(t,m,n,sel,ht,wth)
fprintf('The process has reached %.2f percent\n',100*t/30)
g   = 5;        %弹性系数
mass= 1;        %质量
a0  = 1;        %原长
dmdt = zeros(n*4,1);
for ii = 1:n
    dmdt(ii)    = m(2*n+ii);
    dmdt(n+ii)  = m(3*n+ii);
end
for ii = 1:n
    [y,x] = ind2sub(sel,ii);
    if x==1&&y==1
        F=g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==1&&y==ht
        F=g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==ht
        F=g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==1
        F=g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif x==1
        F = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==wth
        F = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif y==1
        F = g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif y==ht
        F = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    else
        F = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    end
        [fx,fy] = basepotential(m(ii),m(n+ii));
        dmdt(2*n+ii) = (fx + F(1)) / mass;
        dmdt(3*n+ii) = (fy + F(2)) / mass;
end
end
function [fx,fy] = basepotential(x,y)
A   = 0.1;     %外场强度
as  = 0.5;        %外场周期
fx  = -A*sin(x*2*pi/as);
fy  = -A*sin(y*2*pi/as);
end