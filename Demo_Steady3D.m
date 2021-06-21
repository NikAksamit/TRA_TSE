%% Steady 3D Demo
clear all
close all

%% Preparing Data and Path

% Determine where this script is.
folder = fileparts(which('Demo_Steady3D.m'));
% Add that folder and all subfolders to the path.
addpath(genpath(folder));
load('MomentumBarrierField.mat')

NCores=maxNumCompThreads-2;
if isempty(gcp)
    parpool(max([NCores,1]));
end

u_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vx,'linear','none');
v_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vy,'linear','none');
w_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vz,'linear','none');


%%
disp('Running Advection')
x_sparse = linspace(2.26,2.36,199);
y_sparse = linspace(0.33,0.47,200);
z_sparse = 2.55;

[x_gr,y_gr,z_gr]=ndgrid(x_sparse,y_sparse,z_sparse);
[x_mesh,y_mesh,z_mesh]=meshgrid(x_sparse,y_sparse,z_sparse);


tfin=10^(-3);
sVec=linspace(0,tfin,10000);
[xt,yt,zt,time_note,TSE_Bar,TSE,TRA_Bar,TRA] = Advect_and_Calculate_3DSteady(sVec,x_gr(:),y_gr(:),z_gr(:),u_interp,v_interp,w_interp,NCores);

disp(' ')
disp('Plotting')

%%
col=linspace(0.25,0,256)';
gray_col=colormap(gray(256));
figure
ax(1)=subplot(2,2,1);
TRA_Bar=reshape(real(TRA_Bar),size(x_gr));
surf(x_gr,y_gr,TRA_Bar)
view(0,90)
shading interp
title('$\overline{\mathrm{TRA}}_0^{10^{-3}}$','Interpreter','latex')
axis tight
colormap(ax(1),[gray_col(:,1)+col, gray_col(:,2), gray_col(:,3)])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar

ax(2)=subplot(2,2,2);
TRA=reshape(real(TRA),size(x_gr));
surf(x_gr,y_gr,TRA)
view(0,90)
shading interp
title('$\mathrm{TRA}_0^{10^{-3}}$','Interpreter','latex')
axis tight
colormap(ax(2),[gray_col(:,1)+col, gray_col(:,2), gray_col(:,3)])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
hold on
colorbar

ax(3)=subplot(2,2,3);
TSE_Bar=reshape(real(TSE_Bar),size(x_gr));
surf(x_gr,y_gr,TSE_Bar)
view(0,90)
shading interp
title('$\overline{\mathrm{TSE}}_0^{10^{-3}}$','Interpreter','latex')
axis tight
colormap(ax(3),[gray_col(:,1), gray_col(:,2), gray_col(:,3)+col])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar

ax(4)=subplot(2,2,4);
TSE=reshape(real(TSE),size(x_gr));
surf(x_gr,y_gr,TSE)
view(0,90)
shading interp
title('$\mathrm{TSE}_0^{10^{-3}}$','Interpreter','latex')
axis tight
colormap(ax(4),[gray_col(:,1), gray_col(:,2), gray_col(:,3)+col])
daspect([1 1 1])
colorbar
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
drawnow