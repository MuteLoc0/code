clear;
n   = 80;       %粒子数
a0  = 1;        %平衡距离
t0  = [0 100];  %求解区间
y0  = zeros(2*n,1); %单数为速度，偶数为位形
y0(40*2-1)=0.5;
y0(40*2+1)=0.5;
y0(2:2:2*n) = 0:a0:(n-1)*a0;    %粒子初始位形与初始速度
sol = ode45(@(t,x) odefun(t,x,n,a0),t0,y0);
t = linspace(t0(1),t0(2),500)';
x = deval(sol,t);
filename = 'FK1D.gif';
f = figure;
f.set(Color='w')
subplot(2, 1, 1);
line1 = plot(NaN, NaN); %创建一个空的线条对象
axis([0, 1, -0.5, 0.5]);
title('相对位移');
subplot(2, 1, 2);
line2 = plot(NaN, NaN, 'ok', 'MarkerSize', 2); %创建一个空的线条对象
axis([-1, n+1, -0.1, 0.1]);
title('实际位置');
%循环生成每一帧图像并保存
frames  = cell(size(t, 1),1);
for ii = 1:size(t, 1)
    set(line1, 'XData', linspace(0, 1, n), 'YData', x(2:2:2*n, ii)-(0:a0:(n-1)*a0)');
    set(line2, 'XData', x(2:2:2*n, ii), 'YData', zeros(size(x(2:2:2*n, ii))));
    drawnow;

    frame = getframe(f);
    frames{ii} = frame;
    pause(0.05);
end
for ii = 1:numel(frames)
    im = frame2im(frames{ii});
    [imind, cm] = rgb2ind(im, 256);
    %在第一次循环时创建 GIF 文件
    if ii == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

function dxdv = odefun(~,x,n,a0)
A   = 0.08;     %外场强度
g   = 1;        %弹性系数
as  = 1;        %外场周期
m   = 1;        %质量
%dxdv的奇数是x''，偶数是x'. x的奇数是x'，偶数是x. 
dxdv = zeros(2*n,1);
dxdv(1) = (-A*sin(x(2)*2*pi/as)+g*(x(4)-x(2)-a0))/m;
dxdv(2) = x(1);
for i= 2:1:n-1
dxdv(2*i-1) = (-A*sin(x(2*i)*2*pi/as)+g*(x(2*(i+1))+x(2*(i-1))-x(2*i)*2))/m;
dxdv(2*i)   = x(2*i-1);
end
dxdv(2*n-1) = (-A*sin(x(2*n)*2*pi/as)+g*(x(2*n-2)-x(2*n)+a0))/m;
dxdv(2*n)   = x(2*n-1);
end