%% Steady 3D Demo
clear all
close all

%% Preparing Data and Path

% Determine where this script is.
folder = fileparts(which('Demo_Steady3D_Geometry.m'));
% Add that folder and all subfolders to the path.
addpath(genpath(folder));
load('MomentumBarrierField.mat')

NCores=maxNumCompThreads-2;
if isempty(gcp)
    parpool(max([NCores,1]));
end

%%%%%%%%%% Normed Approach for Geometry only %%%%%%%%%%%%%%%
Lap_vx_NORM=Lap_vx./sqrt(Lap_vx.^2+Lap_vy.^2+Lap_vz.^2);
Lap_vy_NORM=Lap_vy./sqrt(Lap_vx.^2+Lap_vy.^2+Lap_vz.^2);
Lap_vz_NORM=Lap_vz./sqrt(Lap_vx.^2+Lap_vy.^2+Lap_vz.^2);
Lap_vx=Lap_vx_NORM;
Lap_vy=Lap_vy_NORM;
Lap_vz=Lap_vz_NORM;
clear Lap_vx_NORM Lap_vy_NORM Lap_vz_NORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Spline and cubic interpolants produce better visualizations at higher
%%%% computational cost
u_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vx,'linear','none');
v_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vy,'linear','none');
w_interp=griddedInterpolant({x_span,y_span,z_span},Lap_vz,'linear','none');


%%
disp('Running Advection')
x_sparse = linspace(2.26,2.36,299);
y_sparse = linspace(0.33,0.47,300);
z_sparse = 2.55;

[x_gr,y_gr,z_gr]=ndgrid(x_sparse,y_sparse,z_sparse);
[x_mesh,y_mesh,z_mesh]=meshgrid(x_sparse,y_sparse,z_sparse);


tfin=0.75;
sVec=linspace(0,tfin,15000);
[xt,yt,zt,time_note,TSE_Bar,~,TRA_Bar] = Advect_and_Calculate_3DSteady(sVec,x_gr(:),y_gr(:),z_gr(:),u_interp,v_interp,w_interp,NCores);

disp(' ')
disp('Plotting')

%%
col=linspace(0.25,0,256)';
gray_col=colormap(gray(256));
figure
ax(1)=subplot(1,2,1);
TRA_Bar=reshape(real(TRA_Bar),size(x_gr));
surf(x_gr,y_gr,TRA_Bar)
view(0,90)
shading interp
title('$\overline{\mathrm{NTRA}}_0^{0.75}(\mathbf{x}_0)$','Interpreter','latex')
axis tight
colormap(ax(1),[gray_col(:,1)+col, gray_col(:,2), gray_col(:,3)])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar
caxis([20 60])

ax(2)=subplot(1,2,2);
TSE_Bar=reshape(real(TSE_Bar),size(x_gr));
surf(x_gr,y_gr,TSE_Bar)
view(0,90)
shading interp
title('$\overline{\mathrm{NTSE}}_0^{0.75}(\mathbf{x}_0)$','Interpreter','latex')
axis tight
colormap(ax(2),[gray_col(:,1), gray_col(:,2), gray_col(:,3)+col])
daspect([1 1 1])
set(gca,'fontsize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
caxis([5 25])
colorbar