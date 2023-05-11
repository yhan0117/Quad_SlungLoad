% Simulation of dynamics with OL control 
clc;clear;close all
addpath Mats
addpath Animation
load("EOM.mat")
load("sys_param.mat")

% animation specs
fps = 20;
duration = 15;

%--------------------------------------System Inputs--------------------------------------%
% const input
input = param(8)*(param(4) + param(5))/(4*param(10));
input = input*ones(1, 4);

% input = input + 0.1*[sin(0.5*t) cos(0.8*t) -sin(0.5*t) -cos(0.8*t)];

%--------------------------------------State Space Transformation--------------------------------------%
ode = subs(ode, u, input);
ode = subs(ode, param, sys_parameters);

V = odeToVectorField(ode);
F = matlabFunction(V,'vars',{'t','Y'});

%--------------------------------------Numerical Solution--------------------------------------%
% alpha x beta y z phi theta psi
x0 = [0 0];
y0 = [0 0];
z0 = [3 0];
a0 = [0 1];
b0 = [pi/3 0];
ph0 = [0 0];
t0 = [0 0];
ps0 = [0 0];

sol = ode45(F, [0 duration], [a0 x0 b0 y0 z0 ph0 t0 ps0]);
t_ = linspace(0, duration, duration*fps);
q_ = deval(sol,t_)';

x = q_(:,3);
y = q_(:,7);
z = q_(:,9);
alpha = q_(:,1); 
beta = q_(:,5);
phi = q_(:,11);
theta = q_(:,13);
psi = q_(:,15);

r_po = subs(r_po, param, sys_parameters);
tmp = zeros(3, length(t_));
R = zeros(3, 3, length(t_));
for i = 1:length(t_)
    tmp(:,i) = subs(r_po, q([1:3,7:8]) , [x(i) y(i) z(i) alpha(i) beta(i)]);
    R(:,:,i) = subs(R_IB, q(4:6), [phi(i) theta(i) psi(i)]);
end

r_ = double(vpa(tmp));
x_po = r_(1,:);
y_po = r_(2,:);
z_po = r_(3,:);

R_IB = double(vpa(R));

clearvars -except t_ x y z x_po y_po z_po R_IB fps
save("Mats\input.mat")
disp("Numerical Solution Done")

%% --------------------------------------Angled Animation--------------------------------------%
clear;clc;clf
load("Mats\input.mat")

% recordning specs
%{
angled
topdown
right
left
%}
POV = 'left';
pb_speed = 1;

% gif name
filename = "Animation\Input\ConstInput\" + POV + ".gif";
delete(filename)

% Axis and labels
figure; hold on
title(sprintf('Trajectory\nTime: %0.2f sec', t_(1)), 'Interpreter', 'Latex');
xlabel('x', 'Interpreter', 'Latex')
ylabel('y', 'Interpreter', 'Latex')
zlabel('z', 'Interpreter', 'Latex')

if strcmp(POV,'angled')
    view(-27.5,25);
elseif strcmp(POV,'topdown')
    view(0,90);
elseif strcmp(POV,'left')
    view(90,0);
elseif strcmp(POV,'right')
    view(0,0);
end

grid minor; axis equal; rotate3d on;
xlim([-2 2]); ylim([-2 2]); zlim([1.5 3.5]); 
set(gcf,'WindowState','fullscreen')

% Plotting with no color to set axis limits
plot3(x,y,z,'Color','none');
plot3(x_po,y_po,z_po,'Color','none');


% Plotting the first iteration
p = plot3(x(1),y(1),z(1),'b');
m = scatter3(x(1),y(1),z(1),'filled','b','square');
p_ = plot3(x_po(1),y_po(1),z_po(1),'r');
m_ = scatter3(x_po(1),y_po(1),z_po(1),10,'filled','r');
L = plot3([x_po(1), x(1)], [y_po(1), y(1)], [z_po(1), z(1)],'k', "LineWidth", 0.1);

% Plot tht 3D box
d = 0.25;
w = 0.25;
h = 0.1/2;
vv          = [d -w -h; d w -h; -d w -h; -d -w -h;
               d -w h; d w h; -d w h; -d -w h];
face         = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
VV = R_IB(:,:,1)*vv.';
VV = VV + [x(1) y(1) z(1)].';
box        = patch('Vertices',VV.','Faces',face,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2);

% Iterating through the length of the time array
for k = 1:length(t_)
    % Updating the line
    p.XData = x(1:k);
    p.YData = y(1:k);
    p.ZData = z(1:k);
    p_.XData = x_po(1:k);
    p_.YData = y_po(1:k);
    p_.ZData = z_po(1:k);

    % Updating the point
    m.XData = x(k); 
    m.YData = y(k);
    m.ZData = z(k);
    m_.XData = x_po(k); 
    m_.YData = y_po(k);
    m_.ZData = z_po(k);

    % Updating the link
    L.XData = [x_po(k), x(k)];
    L.YData = [y_po(k), y(k)];
    L.ZData = [z_po(k), z(k)];
    
    % updating box
    VV = R_IB(:,:,k)*vv.';
    VV = VV + [x(k) y(k) z(k)].';
    set(box,'Faces',face,'Vertices',VV.','FaceColor','black');

    % Updating the title
    title(sprintf('Trajectory\nTime: %0.2f sec', t_(k)),'Interpreter','Latex');

    % Saving the figure
    set(gcf,'WindowState','fullscreen')
    frame = getframe(gcf);

    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime', 1/fps/pb_speed);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 1/fps/pb_speed);
    end
end

close all

