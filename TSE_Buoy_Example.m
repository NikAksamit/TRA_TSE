%%%%%%%%%%%%%%% Input Data %%%%%%%%%%%%%%%%
%%% Longitude: vector of lon positions %%%%
%%% Latitude: vector of lat positions %%%%%
%%% time: vector of GPS sampling times %%%%
%%%         - units of days            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% User Variables %%%%%%%%%%%%%%%%
Integration_Time = 3;                     % Number of days for summing TSE and \bar{TSE}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Generate buoy speed times series from lat/lon %%%%
ds=abs(time(2)-time(1));        % Assuming uniform timestep. Units of days.
window = floor(Integration_Time/ds);
Distances=gc_dist_length(Longitdue(1:end-2),Latitude(1:end-2),Longitude_(3:end),Latitude(3:end));   %km
Speed_sqrd=(Distances.^2)/(ds).^2;    %km/day

%%%% Calculate instantaneous TSE measurements %%%%
TSE_Instant=log((Speed_sqrd(2:end))./(Speed_sqrd(1:end-1)));
TSE_Instant(abs(TSE_Instant)==Inf)=NaN;

%%%% Sum up TSE and \bar{TSE} %%%%
TSE_Bar=1/(time(window+1)-time(1)) * movsum(abs(TSE_Instant),[0 window],'Endpoints','fill');
TSE=1/(time(window+1)-time(1)) * log((Speed_sqrd(window+1:end))./(Speed_sqrd(1:end-window)));
TSE_Bar(end+1:numel(TSE_Instant))=NaN;
TSE(end+1:numel(TSE_Instant))=NaN;

%%%% Great circle distance assuming spherical earth %%%%
function [d]=gc_dist_length(long_initial,lat_initial,long_final,lat_final)
r=6371;     %km
d=r*2*asin(sqrt(sind((lat_initial-lat_final)/2).^2 + cosd(lat_initial).*cosd(lat_final).*sind((long_initial-long_final)/2).^2) );
end