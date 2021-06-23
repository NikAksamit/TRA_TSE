%% Unsteady 2D Demo
clear
close all

%% Preparing Data and Path

% Determine where this script is.
folder = fileparts(which('Demo_Unsteady2D.m'));
% Add that folder and all subfolders to the path.
addpath(genpath(folder));
load('OceanVelocity.mat')

NCores=maxNumCompThreads-2;
if isempty(gcp)
    parpool(max([NCores,1]));
end

u_interp = griddedInterpolant({time,xc,yc},vx,'spline','none');
v_interp = griddedInterpolant({time,xc,yc},vy,'spline','none');

days_transport=30;
tsteps=days_transport*24;     %xxx steps per day

xmin=-2.5;
xmax=5;
ymin=-40;
ymax=-30;
pointsInX = 251;
pointsInY = 250;
x_span = linspace(xmin,xmax,pointsInX);
y_span = linspace(ymin,ymax,pointsInY);

%%%%%%%
[x_grid,y_grid] = ndgrid(x_span,y_span);

tic
t0 = time(1);
tf = t0+days_transport;
tspan=linspace(t0,tf,tsteps);
disp(' ')
[xt,yt,time_note,TSE_Bar_EPS,TRA_Bar_EPS,TSE_EPS] = Advect_and_Calculate_2DUnsteady(tspan,x_grid(:),y_grid(:),u_interp,v_interp,NCores);

TRA_Bar_EPS=real(reshape(TRA_Bar_EPS,size(x_grid)));
TSE_EPS=real(reshape(TSE_EPS,size(x_grid)));
TSE_Bar_EPS=reshape(TSE_Bar_EPS,size(x_grid));

%% Plotting
col=linspace(0.25,0,256)';
gray_col=colormap(gray(256));
figure
ax(1)=subplot(1,3,1);
surf(x_grid,y_grid,TSE_EPS)
view(0,90)
shading interp
title('$\mathrm{TSE}_0^{30}(\mathbf{x}_0,\mathbf{v}_0)$','Interpreter','latex')
axis tight
colormap(ax(1),[gray_col(:,1), gray_col(:,2)+col, gray_col(:,3)])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('Lon','Interpreter','latex')
ylabel('Lat','Interpreter','latex')
colorbar

ax(2)=subplot(1,3,3);
surf(x_grid,y_grid,TRA_Bar_EPS)
view(0,90)
shading interp
title('$\overline{\mathrm{TRA}}_0^{30}(\mathbf{x}_0,\mathbf{v}_0)$','Interpreter','latex')
axis tight
colormap(ax(2),[gray_col(:,1)+col, gray_col(:,2), gray_col(:,3)])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('Lon','Interpreter','latex')
ylabel('Lat','Interpreter','latex')
hold on
colorbar

ax(3)=subplot(1,3,2);
TSE_Bar_EPS=reshape(real(TSE_Bar_EPS),size(x_grid));
surf(x_grid,y_grid,TSE_Bar_EPS)
view(0,90)
shading interp
title('$\overline{\mathrm{TSE}}_0^{30}(\mathbf{x}_0,\mathbf{v}_0)$','Interpreter','latex')
axis tight
colormap(ax(3),[gray_col(:,1), gray_col(:,2), gray_col(:,3)+col])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('Lon','Interpreter','latex')
ylabel('Lat','Interpreter','latex')
colorbar