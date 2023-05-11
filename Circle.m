% Simulation of dynamics with OL control 
clc;clear;close all
addpath Mats
addpath Animation
load("EOM.mat")
load("sys_param.mat")

% animation specs
fps = 20;
duration = 20;

%--------------------------------------Trajectory--------------------------------------%
% Quad follows circular trajectory
% Assumes z0 is normalized to 1 and yaw is 0
T = 7;  % Period of flying around the circle
r = 2;  % Radius of the circle
tau = 0.6;    % Time constant from rest to expected circular traj (99% in 5 time constants)

x_ = r*(1-exp(-t/tau))*cos(2*pi/T*t);
y_ = r*(1-exp(-t/tau))*sin(2*pi/T*t);

% Constrained Equations of Motion
con_eqns = subs(int_ode, q(1:3), [x_ y_ 0]);
con_eqns = subs(con_eqns, param, sys_parameters);

% Numerical solution to alpha beta
V = odeToVectorField(con_eqns(7:8));
F = matlabFunction(V,'vars',{'t','Y'});

a0 = [0 0];
b0 = [eps 0];
sol = ode45(F, [0 duration], [b0 a0]);

t_ = linspace(0, duration, duration*fps);
y = deval(sol,t_)';


% Translational Traj
ddab = zeros(length(t_), 4);
for i = 1:length(t_)
    ddab(i,:) = F(t_(i), y(i,:));
end
beta = y(:,1);
d_beta = y(:,2);
dd_beta = ddab(:,2);
alpha = y(:,3);
d_alpha = y(:,4);
dd_alpha = ddab(:,4);

trans = matlabFunction([x_; y_]);
trans = trans(t_);
x = trans(1,:);
y = trans(2,:);
z = zeros(length(t_), 1);

% Payload traj wrt I
r_po = subs(r_po, param, sys_parameters);
tmp = zeros(length(t_), 3);
for i = 1:length(t_)
    tmp(i,:) = subs(r_po, q([1:3,7:8]) , [x(i) y(i) z(i) alpha(i) beta(i)]);
end

r_ = double(vpa(tmp));
x_po = r_(:,1);
y_po = r_(:,2);
z_po = r_(:,3);


% Rotational EOM
syms theta phi U
con_eqns = subs(con_eqns, [q(4:5) sum(u)], [phi theta U]);

rot = t*zeros(length(t_), 3);
dq = diff(q,t);
ddq = diff(q,t,t);
for i = 1:length(t_)
    rot(i,:) = subs(con_eqns(1:3), ddq(7:8), [dd_alpha(i) dd_beta(i)]);
    rot(i,:) = subs(rot(i,:), dq(7:8), [d_alpha(i) d_beta(i)]);
    rot(i,:) = subs(rot(i,:), q(7:8), [alpha(i) beta(i)]);
    rot(i,:) = subs(rot(i,:), t, t_(i));
end

% roll and pitch
phi_ = zeros(length(t_), 1);
theta_ = zeros(length(t_), 1);
psi_ = zeros(length(t_), 1);
U_ = zeros(length(t_), 1);
for i = 1:length(t_)
    fun = rot(i,:);
    X1 = U*sin(theta);
    X2 = U*sin(phi)*cos(theta);
    X3 = U*cos(phi)*cos(theta);
    X1 = vpa(rhs(isolate(fun(:,1), X1)));
    X2 = vpa(rhs(isolate(fun(:,2), X2)));
    X3 = vpa(rhs(isolate(fun(:,3), X3)));
    
    phi_(i) = atan2(X2, X3);
    theta_(i) = atan2(X1, X3/cos(phi_(i)));
    U_(i) = X1/sin(theta_(i));
end

% Rotation matrix at each time step
R = zeros(3, 3, length(t_));
for i = 1:length(t_)
    R(:,:,i) = subs(R_IB, q(4:6), [phi_(i) theta_(i) psi_(i)]);
end

R_IB = double(vpa(R));
clearvars -except t_ x y z x_po y_po z_po R_IB fps U_
save("Mats\circle.mat")
disp("Numerical Solution Done")
%%

%% --------------------------------------Animation--------------------------------------%
clear;clc;clf

load("Mats\circle.mat")

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
filename = "Animation\Circle\" + POV + ".mp4";
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
xlim([-3 3]); ylim([-3 3]); zlim([-1.5 0.5]); 
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

    % Delay
    pause(0.01)

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
