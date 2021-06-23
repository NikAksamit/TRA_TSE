% Input arguments:
% tspan                        : Advection timesteps (e.g. tspan=linspace(t_initial,t_final,100)
% xx,yy,zz                        : Initial positions for advection
% U_Interp,V_Interp,W_Interp            : Velocity field interpolants with inputs (tspan,xx,yy)
% NCores                       : Number of Cores for parpool

% Output arguments:
%   xt,yt,zt    : xx-component, yy-component, zz-component of trajectory final position
%   time_note    : tspan time if a trajectory left interpolant domain
%   TSE_Bar_EPS,TRA_Bar_EPS,TSE_EPS    : Single trajectory metrics from
%   section IV in [1] (Theorem 3)

%--------------------------------------------------------------------------
% Author: Nikolas Aksamit  naksamit@ethz.ch
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References:
% ﻿[1] Haller, G., Aksamit, N. O., & Bartos, A. P. E. (2021). Quasi-Objective Coherent Structure Diagnostics from Single Trajectories.
%  Chaos, 31, 043131-1–17. https://doi.org/10.1063/5.0044151
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate extended phase space (EPS) TRA and TSE diagnostics for unsteady
%%% 3D flows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function     [xt,yt,zt,time_note,TSE_Bar_EPS,TRA_Bar_EPS,TSE_EPS] = Advect_and_Calculate_3DUnsteady(tspan,xx,yy,zz,U_Interp,V_Interp,W_Interp,NCores)

%%% Calculate v_0 as spatio-temporal mean of the flow in your domain of
%%% interest
v_0=mean(sqrt(U_Interp.Values.^2+V_Interp.Values.^2));
%%%% v_0=0; %%%% Equivalent to steady version of TRA/TSE

if isempty(gcp)
    parpool(NCores);
end

% Shared variables
Np=numel(xx);
cpu_num = min(NCores,Np);
id = ceil( linspace(0,Np,cpu_num+1) );

% Shared variables
xt = zeros(2,length(xx));
yt = zeros(2,length(yy));
zt = zeros(2,length(zz));

xt(1,:) = xx;
yt(1,:) = yy;
zt(1,:) = zz;

spmd(cpu_num)
    tic
    Range = id(labindex)+1:id(labindex+1);
    x_spmd=xt(1,Range);
    y_spmd=yt(1,Range);
    z_spmd=zt(1,Range);
    
    time_note=zeros(size(x_spmd));
    TSE_EPS_spmd=zeros(size(x_spmd));
    TRA_EPS_spmd=zeros(size(x_spmd));
    
    for s=1:numel(tspan)-1
        %%%% This only records the last point, not entire path of particles
        ds = tspan(s+1) - tspan(s);
        
        [UK1,VK1,WK1] = r_prime(tspan(s),x_spmd(1,:),y_spmd(1,:),z_spmd(1,:),U_Interp,V_Interp,W_Interp);
        
        xx = x_spmd(1,:) + 0.5 * ds * UK1;
        yy = y_spmd(1,:) + 0.5 * ds * VK1;
        zz = z_spmd(1,:) + 0.5 * ds * WK1;
        [UK2,VK2,WK2] = r_prime(tspan(s),xx,yy,zz,U_Interp,V_Interp,W_Interp);
        
        
        xx = x_spmd(1,:) + 0.5 * ds * UK2;
        yy = y_spmd(1,:) + 0.5 * ds * VK2;
        zz = z_spmd(1,:) + 0.5 * ds * WK2;
        [UK3,VK3,WK3] = r_prime(tspan(s),xx,yy,zz,U_Interp,V_Interp,W_Interp);
        
        
        xx = x_spmd(1,:) + ds * UK3;
        yy = y_spmd(1,:) + ds * VK3;
        zz = z_spmd(1,:) + ds * WK3;
        [UK4,VK4,WK4] = r_prime(tspan(s),xx,yy,zz,U_Interp,V_Interp,W_Interp);
        
        %increment in trajectories (RK4 displacement)
        deltax(2,:) = ds / 6 * (UK1 + 2 * UK2 + 2 * UK3 + UK4);
        deltay(2,:) = ds / 6 * (VK1 + 2 * VK2 + 2 * VK3 + VK4);
        deltaz(2,:) = ds / 6 * (WK1 + 2 * WK2 + 2 * WK3 + WK4);
        
        %update particle positions
        x_spmd(2,:) = x_spmd(1,:) + deltax(2,:);
        y_spmd(2,:) = y_spmd(1,:) + deltay(2,:);
        z_spmd(2,:) = z_spmd(1,:) + deltaz(2,:);

        %%%% If particle leaves domain, keep record of last position.
        %%%% Record this time of leaving domain in time_note.
        x_spmd(1,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=x_spmd(2,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
        y_spmd(1,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=y_spmd(2,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
        z_spmd(1,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=z_spmd(2,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
        time_note(1,(isnan(x_spmd(2,:)) & time_note(1,:)==0) | (isnan(y_spmd(2,:)) & time_note(1,:)==0) | (isnan(z_spmd(2,:)) & time_note(1,:)==0))=s;
        
        if s==1
            smooth_vx=deltax(2,:)'/ds;
            smooth_vy=deltay(2,:)'/ds;
            smooth_vz=deltaz(2,:)'/ds;
            
            V1=[smooth_vx,smooth_vy,smooth_vz];
            Speed_sqrd(1,:)=sum(V1.^2,2);
            
            deltax(1,:) = deltax(2,:);
            deltay(1,:) = deltay(2,:);
            deltaz(1,:) = deltaz(2,:);
            
            %%% Save this first veloicty measurement for TRA measurement at end
            V1_spmd=V1;
        end
        
        if s>1
            smooth_vx=deltax(2,:)'/ds;
            smooth_vy=deltay(2,:)'/ds;
            smooth_vz=deltaz(2,:)'/ds;
            
            V2=[smooth_vx,smooth_vy,smooth_vz];
            Speed_sqrd(2,:)=sum(V2.^2,2);
            
            
            TSE_Inst_EPS=abs(log(sqrt(Speed_sqrd(2,:)+v_0^2)./sqrt(Speed_sqrd(1,:)+v_0^2)));
            TSE_EPS_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)))=TSE_EPS_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)))+TSE_Inst_EPS(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)));
            
            TRA_Inst_EPS=acos((sum(V2.*V1,2)+v_0^2)./(sqrt(sum(V2.^2,2)+v_0^2).*sqrt(sum(V1.^2,2)+v_0^2)))';
            TRA_EPS_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)))=TRA_EPS_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)))+TRA_Inst_EPS(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)) & ~isnan(z_spmd(2,:)));
            
            deltax(1,:) = deltax(2,:);
            deltay(1,:) = deltay(2,:);
            deltaz(1,:) = deltaz(2,:);
            Speed_sqrd(1,:)=Speed_sqrd(2,:);
            V1=V2;
            V2_spmd=V2;
        end
        
    end
    
    toc
end

V1 = cat(1,V1_spmd{:});
V2 = cat(1,V2_spmd{:});
TSE_EPS=1/(tspan(end)-tspan(1))*log(sqrt((sum(V2.*V2,2)+v_0.^2)./(sum(V1.*V1,2).^2+v_0.^2)));

xt = cat(2,x_spmd{:});
yt = cat(2,y_spmd{:});
time_note = cat(2,time_note{:});

TSE_Bar_EPS=cat(2,TSE_EPS_spmd{:})/(tspan(end)-tspan(1));
TRA_Bar_EPS=cat(2,TRA_EPS_spmd{:});

clear x_spmd y_spmd

% disp('Advection Complete')
end


function [uk,vk,wk] = r_prime(s,xx,yy,zz,U_Interp,V_Interp,W_Interp)

Bly = 1+0*xx;
bounds.x=U_Interp.GridVectors{1};
bounds.y=U_Interp.GridVectors{2};
bounds.z=U_Interp.GridVectors{3};

Bly(yy>max(bounds.y)) = 0;      %Freeze when outside domain of interpolant
Bly(yy<min(bounds.y)) = 0;
Bly(xx>max(bounds.x)) = 0;
Bly(xx<min(bounds.x)) = 0;
Bly(zz>max(bounds.z)) = 0;
Bly(zz<min(bounds.z)) = 0;

uk=U_Interp(double(s)*ones(size(xx)),double(xx),double(yy),double(zz));
uk=uk.*Bly;

vk=V_Interp(double(s)*ones(size(xx)),double(xx),double(yy),double(zz));
vk=vk.*Bly;

wk=W_Interp(double(s)*ones(size(xx)),double(xx),double(yy),double(zz));
wk=wk.*Bly;

end
