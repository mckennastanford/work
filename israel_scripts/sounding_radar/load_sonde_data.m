%% load sonde data.

function Sonde = load_sonde_data(Flist_struct) 

siteName = 'awr';
config_sites;

Path = Flist_struct.path;
fileName = Flist_struct.name;
% Path = '/bear/s2/data/AWARE/SONDE/';
% fileName = 'awrsondewnpnS1.b1.20160103.055900.cdf';
% ncdisp([Path fileName])

% dimentional orientation vectors.
ttmp = ncread([Path fileName], 'time_offset'); %in seconds (since 2015-11-30) since midnight.
bttmp = ncread([Path fileName], 'base_time'); %Base time in Epoch (1970-1-1) in seconds).
Sonde.time = datenum(1970,1,1,0,0,double(ttmp) + double(bttmp)); %measurement time in datenum format.
clear ttmp bttmp

% Sonde measurements.
Sonde.dewpoint_temp = double(ncread([Path fileName], 'dp')); %in Celcius.
Sonde.drybulb_temp = double(ncread([Path fileName], 'tdry')); %in Celcius.
Sonde.RH = double(ncread([Path fileName], 'rh')); % in '%'.
Sonde.pressure = double(ncread([Path fileName], 'pres')); %hPa
Sonde.wind_direction = double(ncread([Path fileName], 'deg')); %wind direction in degrees.
Sonde.v_wind = double(ncread([Path fileName], 'v_wind')); %in m/s
Sonde.u_wind = double(ncread([Path fileName], 'u_wind')); %in m/s
Sonde.wind_speed = double(ncread([Path fileName], 'wspd')); %in m/s
% Sonde.ascent_rate = double(ncread([Path fileName], 'asc')); %m/s

% radiosonde location.
Sonde.lat = double(ncread([Path fileName], 'lat')); % degrees_N
Sonde.lon = double(ncread([Path fileName], 'lon')); % degrees_E
Sonde.alt = double(ncread([Path fileName], 'alt')); % altitude in m above MSL.

%'nan'ing bad data.
Sonde.dewpoint_temp(Sonde.dewpoint_temp == -9999) = nan;        Sonde.drybulb_temp(Sonde.drybulb_temp == -9999) = nan;
Sonde.RH(Sonde.RH == -9999) = nan;                              Sonde.pressure(Sonde.pressure == -9999) = nan;
Sonde.wind_direction(Sonde.wind_speed == -9999) = nan;          Sonde.wind_direction(Sonde.wind_speed == -9999) = nan;
Sonde.v_wind(Sonde.v_wind == -9999) = nan;                      Sonde.u_wind(Sonde.u_wind == -9999) = nan;
% Sonde.ascent_rate(Sonde.ascent_rate == -9999) = nan;

% quality check (qc) sub structure (value of 0 (zero) means that the data passed all the tests. for other qc information per parameter, look at the netcdf files).
%     Sonde.qc.pressure = ncread([Path fileName], 'qc_pres');
%     Sonde.qc.drybulb_temp = ncread([Path fileName], 'qc_tdry');
% % Sonde.qc.dewpoint_temp = ncread([Path fileName], 'qc_dp');
% % Sonde.qc.wind_speed = ncread([Path fileName], 'qc_wspd');
% % Sonde.qc.wind_direction = ncread([Path fileName], 'qc_deg');
%     Sonde.qc.RH = ncread([Path fileName], 'qc_rh');
%     Sonde.qc.u_wind = ncread([Path fileName], 'qc_u_wind');
%     Sonde.qc.v_wind = ncread([Path fileName], 'qc_v_wind');
% % Sonde.qc.ascent_rate = ncread([Path fileName], 'qc_asc');