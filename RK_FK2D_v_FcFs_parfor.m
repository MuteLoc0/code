clear;
clc;
hold off
ht  = 16;       %行数
wth = 16;       %列数
sel = [ht wth]; %方便索引，使用sub2ind(sel,x,y)
n   = ht*wth;   %方便索引，总粒子数
t0  = [0 10];   %求解区间
m0  = zeros(n*4,1);    %初始位置与速度. 分配规则为：按列依次对应，分成四大块
A_e=1;%外力大小
alpha=pi/6;%外力夹角
%分别为x轴位置、y轴位置、x速度、y速度.例如：
%   1   3;  5   7; 9    11; 13  15
%   2   4;  6   8; 10   12; 14  16
%分别为粒子对应x位置、y位置、x速度、y速度四个矩阵
for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0(ii)  = w;
    m0(n+ii)= l;
end

deltaa=0.08;
numn=24;
%计算平均速度

tic
% 预先定义共享数组用于存储 vi 和 vd
vi = zeros(2, numn); % 创建一个 2 行 numn 列的数组来存储每个 jj 的 vi
vd = zeros(2, numn); % 创建一个 2 行 numn 列的数组来存储每个 jj 的 vd

parfor jj = 1:2
    m0_cache = zeros(2*numn, n*4); % 缓存所有初始状态
    for ii = 1:n    %该循环设置粒子初始位置
    [l,w]   = ind2sub(sel,ii);
    m0_cache(1,ii)  = w;
    m0_cache(1,n+ii)= l;
    m0_cache(2*numn,ii)  = w;
    m0_cache(2*numn,n+ii)= l;
    end
    if jj == 1
        local_vi = zeros(numn, 1); % 每个并行工作者都有自己的 local_vi
        for ii = 1:numn
            A_e = ii * deltaa;
            % 调试输出初始条件
            fprintf('jj = %d, ii = %d, A_e = %f, m0_cache(ii, :) = %s\n', jj, ii, A_e, mat2str(m0_cache(ii, :)));
            [t, m] = ode45(@(t, m) odefun(t, m, n, sel, ht, wth, A_e, alpha), [1, 24], m0_cache(ii, :));
            if any(isnan(m(:)))
                warning('NaN detected in ODE solution for jj = %d, ii = %d', jj, ii);
            end
            m0_cache(ii+1, :) = m(end, :); % 缓存当前状态作为下一个的初始状态
            local_vi(ii) = mean((mean(m(:, 2*n+1:3*n), 2).^2 + mean(m(:, 3*n+1:4*n), 2).^2).^(0.5), 'all');
        end
        vi(jj, :) = local_vi; % 将 local_vi 的值存储到共享数组中
    end

    if jj == 2
        % 第二部分循环
        local_vd = zeros(numn, 1); % 每个并行工作者都有自己的 local_vd
        for ii = 1:numn
            A_e = deltaa * (numn + 1 - ii);
            % 调试输出初始条件
            fprintf('jj = %d, ii = %d, A_e = %f, m0_cache(numn + 1 - ii, :) = %s\n', jj, ii, A_e, mat2str(m0_cache(2*numn + 1 - ii, :)));
            [t, m] = ode45(@(t, m) odefun(t, m, n, sel, ht, wth, A_e, alpha), [1, 26], m0_cache(2*numn + 1 - ii, :));
            if any(isnan(m(:)))
                warning('NaN detected in ODE solution for jj = %d, ii = %d', jj, ii);
            end
            m0_cache(2*numn - ii, :) = m(end, :);
            local_vd(numn + 1 - ii) = mean((mean(m(:, 2*n+1:3*n), 2).^2 + mean(m(:, 3*n+1:4*n), 2).^2).^(0.5), 'all');
        end
        vd(jj, :) = local_vd;
    end
end

toc

% 通过 jj 的值来绘制 vi 和 vd
plot(vi(1, :));
hold on
plot(vd(2, :));

function dmdt = odefun(t,m,n,sel,ht,wth,A_e,alpha)
g   = 1;        %弹性系数
mass= 0.1;      %质量
a0=1;         %原长
A = 0.5;       %基底势震荡强度
as  = 0.5;        %外场周期
dmdt = zeros(n*4,1);
F_e=A_e*[cos(alpha),sin(alpha)];
gamma=0.5;%阻尼系数（欠阻尼）
for ii = 1:n
    dmdt(ii)    = m(2*n+ii);
    dmdt(n+ii)  = m(3*n+ii);
end
for ii = 1:n
    [y,x] = ind2sub(sel,ii);
    if x==1&&y==1
        F_i = g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==1&&y==ht
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)]-a0*normalize([m(ii+ht)-m(ii) m(n+ii+ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==ht
        F_i = g*([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)]-a0*normalize([m(ii-1)-m(ii) m(n+ii-1)-m(n+ii)],'norm')+ ...
            [m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)]-a0*normalize([m(ii-ht)-m(ii) m(n+ii-ht)-m(n+ii)],'norm'));
    elseif x==wth&&y==1
        F_i = g*([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)]-a0*normalize([m(ii+1)-m(ii) m(n+ii+1)-m(n+ii)],'norm')+ ...
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
        dmdt(2*n+ii) = (-A*sin(m(ii)*2*pi/as)...
            + F_i(1)+F_e(1)-gamma*m(2*n+ii)) / mass;
        dmdt(3*n+ii) = (-A*sin(m(n+ii)*2*pi/as)...
            + F_i(2)+F_e(2)-gamma*m(3*n+ii)) / mass;
end
end
