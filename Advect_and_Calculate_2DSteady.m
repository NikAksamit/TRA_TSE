% Input arguments:
    % tspan                        : Advection timesteps (e.g. tspan=linspace(t_initial,t_final,100)
    % xx,yy                        : Initial positions for advection
    % U_Interp,V_Interp            : Velocity field interpolants with inputs (tspan,xx,yy)
    % NCores                       : Number of Cores for parpool

    % Output arguments:
    %   xt,yt    : xx-component, yy-component of trajectory final position
    %   time_note    : tspan time if a trajectory left interpolant domain
    %   TSE_Bar,TSE,TRA_Bar,TRA    : Single trajectory metrics from
    %   section II and III in [1]

    %--------------------------------------------------------------------------
% Author: Nikolas Aksamit  naksamit@ethz.ch
%--------------------------------------------------------------------------

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References:
% ﻿[1] Haller, G., Aksamit, N. O., & Bartos, A. P. E. (2021). Quasi-Objective Coherent Structure Diagnostics from Single Trajectories. 
%  Chaos, 31, 043131-1–17. https://doi.org/10.1063/5.0044151
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate TRA and TSE diagnostics for steady 2D flows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function     [xt,yt,time_note,TSE_Bar,TSE,TRA_Bar,TRA] = Advect_and_Calculate_2DSteady(tspan,xx,yy,U_Interp,V_Interp,NCores)

%%%%% Start parallel pool with NCores. Make sure the size of pool
%%%%% >= number requested by spmd
if isempty(gcp)
    parpool(NCores);
end


Np=numel(xx);
cpu_num = min(NCores,Np);
id = ceil( linspace(0,Np,cpu_num+1) );

% Shared variables
xt = zeros(2,length(xx));
yt = zeros(2,length(yy));

xt(1,:) = xx;
yt(1,:) = yy;

spmd(cpu_num)
    tic
    
    %%%% Split up initial conditions for separate cores
    Range = id(labindex)+1:id(labindex+1);
    x_spmd=xt(1,Range);
    y_spmd=yt(1,Range);
    
    time_note=zeros(size(x_spmd));
    TSE_spmd=zeros(size(x_spmd));
    TRA_spmd=zeros(size(x_spmd));
    TSE_spmd_NM=zeros(size(x_spmd));
    
    for s=1:numel(tspan)-1
        
        %%%% This only records the last point, not entire path of particles
        ds = tspan(s+1) - tspan(s);
        
        [UK1,VK1] = r_prime(x_spmd(1,:),y_spmd(1,:),U_Interp,V_Interp);
        
        xx = x_spmd(1,:) + 0.5 * ds * UK1;
        yy = y_spmd(1,:) + 0.5 * ds * VK1;
        [UK2,VK2] = r_prime(xx,yy,U_Interp,V_Interp);
        
        
        xx = x_spmd(1,:) + 0.5 * ds * UK2;
        yy = y_spmd(1,:) + 0.5 * ds * VK2;
        [UK3,VK3] = r_prime(xx,yy,U_Interp,V_Interp);
        
        
        xx = x_spmd(1,:) + ds * UK3;
        yy = y_spmd(1,:) + ds * VK3;
        [UK4,VK4] = r_prime(xx,yy,U_Interp,V_Interp);
        
        %increment in trajectories (RK4 displacement)
        deltax(2,:) = ds / 6 * (UK1 + 2 * UK2 + 2 * UK3 + UK4);
        deltay(2,:) = ds / 6 * (VK1 + 2 * VK2 + 2 * VK3 + VK4);
        
        %update particle positions
        x_spmd(2,:) = x_spmd(1,:) + deltax(2,:);
        y_spmd(2,:) = y_spmd(1,:) + deltay(2,:);
        
        
        
        %%%% If particle leaves domain, keep record of last position.
        %%%% Record this time of leaving domain in time_note.
        x_spmd(1,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=x_spmd(2,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
        y_spmd(1,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=y_spmd(2,~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
        time_note(1,(isnan(x_spmd(2,:)) & time_note(1,:)==0) | (isnan(y_spmd(2,:)) & time_note(1,:)==0))=s;
        
        if s==1
            smooth_vx=deltax(2,:)'/ds;
            smooth_vy=deltay(2,:)'/ds;

            V1=[smooth_vx,smooth_vy];
            Speed_sqrd(1,:)=sum(V1.^2,2);
            
            deltax(1,:) = deltax(2,:);
            deltay(1,:) = deltay(2,:);
            
            %%% Save this first veloicty measurement for TRA measurement at end 
            V1_spmd=V1;
        end
        
        
        %%% After two timesteps, now begin calculating TSE and TRA 
        if s>1
            smooth_vx=deltax(2,:)'/ds;
            smooth_vy=deltay(2,:)'/ds;

            V2=[smooth_vx,smooth_vy];
            Speed_sqrd(2,:)=sum(V2.^2,2);
            

            %%%% Calculate instantaneous TSE and TRA values. Add to sum if particle did not leave domain.
            TSE_Inst=abs(log(sqrt(Speed_sqrd(2,:))./sqrt(Speed_sqrd(1,:))));
            TSE_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=TSE_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))+TSE_Inst(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
            TSE_Inst_NM=log(sqrt(Speed_sqrd(2,:))./sqrt(Speed_sqrd(1,:)));
            TSE_spmd_NM(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=TSE_spmd_NM(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))+TSE_Inst_NM(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
            
            TRA_Inst=acos(sum(V2.*V1,2)./(sqrt(sum(V2.^2,2)).*sqrt(sum(V1.^2,2))))';
            TRA_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))=TRA_spmd(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)))+TRA_Inst(~isnan(x_spmd(2,:)) & ~isnan(y_spmd(2,:)));
            
            deltax(1,:) = deltax(2,:);
            deltay(1,:) = deltay(2,:);
            Speed_sqrd(1,:)=Speed_sqrd(2,:);
            V1=V2;
            V2_spmd=V2;
        end
    end
    
    
    
    TSE_spmd=TSE_spmd(1,:);
    TRA_spmd=TRA_spmd(1,:);
    TSE_spmd_NM=TSE_spmd_NM(1,:);
    toc
end
%
V1 = cat(1,V1_spmd{:});
V2 = cat(1,V2_spmd{:});
V1=V1./(sqrt(sum(V1.*V1,2)));
V2=V2./(sqrt(sum(V2.*V2,2)));
xt = cat(2,x_spmd{:});
yt = cat(2,y_spmd{:});
time_note = cat(2,time_note{:});
TSE=cat(2,TSE_spmd_NM{:})/(tspan(end)-tspan(1));
TSE_Bar=cat(2,TSE_spmd{:})/(tspan(end)-tspan(1));
TRA_Bar=cat(2,TRA_spmd{:});

TRA=real(acos(V1(:,1).*V2(:,1)+V1(:,2).*V2(:,2)));

clear x_spmd y_spmd

disp('Advection and Single Trajectory Metric Calculations Complete')
end


function [uk,vk] = r_prime(xx,yy,U_Interp,V_Interp)


Bly = 1+0*xx;
bounds.x=U_Interp.GridVectors{1};
bounds.y=U_Interp.GridVectors{2};

Bly(yy>max(bounds.y)) = 0;      %Freeze when outside domain of interpolant
Bly(yy<min(bounds.y)) = 0;
Bly(xx>max(bounds.x)) = 0;
Bly(xx<min(bounds.x)) = 0;

uk=U_Interp(double(xx),double(yy));
uk=uk.*Bly;

vk=V_Interp(double(xx),double(yy));
vk=vk.*Bly;
end