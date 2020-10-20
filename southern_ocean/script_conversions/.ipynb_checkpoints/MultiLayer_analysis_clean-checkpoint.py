clear all
%% Load long-term data.
Site = 0; % 0 - nsa, 1 - awr.
if Site == 0
    STR = 'nsaC1';
    Site_name = 'North Slope, Alaska';
    Site_loc = [71.32 -156.61];
    Search_string = ['nsa*b1*.cdf'];
    Sonde_path = '/bear/s1/data/nsa/hsrl.wisconsin/sonde/'; % ARM file directory
    Site_alt = 7; % for Barrow.
elseif Site == 1
    STR = 'awrM1';
    Site_name = 'McMurdo Station, Antarctica';
    Site_loc = [-77.85 166.65];
    Search_string = ['awr*b1*.cdf'];
    Sonde_path = '/bear/s2/data/AWARE/SONDE/'; % ARM file directory
    Site_alt = 70; % for McM.
end
Res_2_use = 15; % in m
Max_cloud_height = 10e3; % highest altitude to analyze in this study.
[Flist_struct] = give_me_files_and_subfolders(Search_string, Sonde_path );
Alt = (0:Res_2_use:Max_cloud_height)'+ Site_alt;
Num_layers = length(Alt);
Num_profiles = length(Flist_struct);

Sonde.release_t = NaN(1, Num_profiles);
Sonde.vert_res_mean = NaN(1, Num_profiles);
Sonde.asc_mean = NaN(1, Num_profiles);
Sonde.time = NaN(Num_layers, Num_profiles);
Sonde.lat = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.lon = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.p = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.T = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.q = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.RH = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.RH_i = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.theta = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.theta_e = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.u = single(ones(Num_layers, Num_profiles).* -9999);
Sonde.v = single(ones(Num_layers, Num_profiles).* -9999);
for ii = Num_profiles: -1: 1
    if mod(Num_profiles - ii, 100) == 0;    disp(['Profile ', num2str(Num_profiles - ii + 1), ' out of ', num2str(Num_profiles), ' (', num2str((Num_profiles - ii)/Num_profiles*100, '%.2f'), '% done!)']);     end
    Sondetmp = load_sonde_data(Flist_struct(ii));
    % Try two types of field names for the ascent rate
    try        Sondetmp.asc_mean = ncread([Flist_struct(ii).path Flist_struct(ii).name], 'asc');        catch;      end
    try        Sondetmp.asc_mean = ncread([Flist_struct(ii).path Flist_struct(ii).name], 'wcmp');       catch;      end
    [~, Unique_loc] = unique(Sondetmp.alt); % unique sounding values.
    if length(Unique_loc) > 1
        Sonde.release_t(ii) = Sondetmp.time(1);
        Sonde.vert_res_mean(ii) = nanmean(diff(Sondetmp.alt(Unique_loc)));
        try        Sonde.asc_mean(ii) = nanmean(Sondetmp.asc_mean(Unique_loc));        catch;      warning('likely no asc data'); end
        Sonde.time(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.time(Unique_loc), Alt);
        Sonde.lat(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.lat(Unique_loc), Alt);
        Sonde.lon(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.lon(Unique_loc), Alt);
        Sonde.p(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.pressure(Unique_loc), Alt);
        Sonde.T(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.drybulb_temp(Unique_loc), Alt);
        Sonde.RH(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.RH(Unique_loc), Alt);
        Moretmp = calculate_theta_and_more(Sonde.T(:, ii), Sonde.p(:, ii), Sonde.RH(:, ii));
        Sonde.q(:, ii) = Moretmp.q;
        Sonde.RH_i(:, ii) = Moretmp.RH_i;
        Sonde.theta(:, ii) = Moretmp.Theta;
        Sonde.theta_e(:, ii) = Moretmp.Theta_e;
        Sonde.u(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.u_wind(Unique_loc), Alt);
        Sonde.v(:, ii) = interp1(Sondetmp.alt(Unique_loc), Sondetmp.v_wind(Unique_loc), Alt);
    end
    clear Sondetmp Moretmp
end
% no data or data only outside the examined range.
Valid_locs = find(~(isnan(Sonde.release_t) | sum(~isnan(Sonde.T)) == 0));
Sonde.time = Sonde.time(:, Valid_locs);
Sonde.lat = Sonde.lat(:, Valid_locs);               Sonde.lon = Sonde.lon(:, Valid_locs);
Sonde.p = Sonde.p(:, Valid_locs);                   Sonde.T = Sonde.T(:, Valid_locs);
Sonde.q = Sonde.q(:, Valid_locs);
Sonde.RH = Sonde.RH(:, Valid_locs);                 Sonde.RH_i = Sonde.RH_i(:, Valid_locs);
Sonde.theta = Sonde.theta(:, Valid_locs);           Sonde.theta_e = Sonde.theta_e(:, Valid_locs);
Sonde.u = Sonde.u(:, Valid_locs);                   Sonde.v = Sonde.v(:, Valid_locs);
Sonde.release_t = Sonde.release_t(Valid_locs);
Sonde.vert_res_mean = Sonde.vert_res_mean(Valid_locs);
Sonde.asc_mean = Sonde.asc_mean(Valid_locs);
% order the release time.
[~, Order] = sort(Sonde.release_t);
Sonde.time = Sonde.time(:, Order);
Sonde.lat = Sonde.lat(:, Order);               Sonde.lon = Sonde.lon(:, Order);
Sonde.p = Sonde.p(:, Order);                   Sonde.T = Sonde.T(:, Order);
Sonde.q = Sonde.q(:, Order);
Sonde.RH = Sonde.RH(:, Order);                 Sonde.RH_i = Sonde.RH_i(:, Order);
Sonde.theta = Sonde.theta(:, Order);           Sonde.theta_e = Sonde.theta_e(:, Order);
Sonde.u = Sonde.u(:, Order);                   Sonde.v = Sonde.v(:, Order);
Sonde.release_t = Sonde.release_t(Order);
Sonde.vert_res_mean = Sonde.vert_res_mean(Order);
Sonde.asc_mean = Sonde.asc_mean(Order);
Num_profiles = length(Order);

%--------------------------------------------------------------------------------------------
%% Export raw fields to a NetCDF file
Release_time = (Sonde.release_t - datenum(1970,1,1)).*24.*3600;
Time = (Sonde.time - datenum(1970,1,1)).*24.*3600;
Time(Time < 0) = -9999;     Release_time(Release_time < 0) = -9999;
Sonde.time(isnan(Sonde.time)) = -9999;
Sonde.lat(isnan(Sonde.lat)) = -9999;               Sonde.lon(isnan(Sonde.lon)) = -9999;
Sonde.p(isnan(Sonde.p)) = -9999;                   Sonde.T(isnan(Sonde.T)) = -9999;
Sonde.q(isnan(Sonde.q)) = -9999;
Sonde.RH(isnan(Sonde.RH)) = -9999;                 Sonde.RH_i(isnan(Sonde.RH_i)) = -9999;
Sonde.theta(isnan(Sonde.theta)) = -9999;           Sonde.theta_e(isnan(Sonde.theta_e)) = -9999;
Sonde.u(isnan(Sonde.u)) = -9999;                   Sonde.v(isnan(Sonde.v)) = -9999;
Sonde.release_t(isnan(Sonde.release_t)) = -9999;
Sonde.vert_res_mean(isnan(Sonde.vert_res_mean)) = -9999;
Sonde.asc_mean(isnan(Sonde.asc_mean)) = -9999;

save_path = ['~/'];
nc_Filename = [STR, '_full_sounding.nc'];
nccreate([save_path, nc_Filename], 'site_altitude', 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'site_altitude', Site_alt)
ncwriteatt([save_path, nc_Filename], 'site_altitude', 'Units', 'meters')
ncwriteatt([save_path, nc_Filename], 'site_altitude', 'Description', 'Site altitude above MSL')
nccreate([save_path, nc_Filename], 'site_latitude', 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'site_latitude', Site_loc(1))
ncwriteatt([save_path, nc_Filename], 'site_latitude', 'Units', 'degrees N')
ncwriteatt([save_path, nc_Filename], 'site_latitude', 'Description', 'Site latitude')
nccreate([save_path, nc_Filename], 'site_longitude', 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'site_longitude', Site_loc(2))
ncwriteatt([save_path, nc_Filename], 'site_longitude', 'Units', 'degrees E')
ncwriteatt([save_path, nc_Filename], 'site_longitude', 'Description', 'Site longitude')
nccreate([save_path, nc_Filename], 'height', 'dimensions', {'Num_layers', Num_layers}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'height', Alt - Site_alt)
ncwriteatt([save_path, nc_Filename], 'height', 'Units', 'meters')
ncwriteatt([save_path, nc_Filename], 'height', 'description', 'height AGL (equivalent to range in cloud mask files)')
nccreate([save_path, nc_Filename], 'release_time', 'dimensions', {'Num_profiles', Num_profiles}, 'datatype', 'double', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'release_time', Release_time)
ncwriteatt([save_path, nc_Filename], 'release_time', 'Units', 'Seconds')
ncwriteatt([save_path, nc_Filename], 'release_time', 'Description', ['Radiosonde release time (seconds since 1970-01-01)'])
nccreate([save_path, nc_Filename], 'vert_res_mean', 'dimensions', {'Num_profiles', Num_profiles}, 'datatype', 'double', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'vert_res_mean', Sonde.vert_res_mean)
ncwriteatt([save_path, nc_Filename], 'vert_res_mean', 'Units', 'meters')
ncwriteatt([save_path, nc_Filename], 'vert_res_mean', 'Description', ['Mean vertical resolution (a function of the ascent rate and the temporal resolution)'])
nccreate([save_path, nc_Filename], 'asc_mean', 'dimensions', {'Num_profiles', Num_profiles}, 'datatype', 'double', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'asc_mean', Sonde.asc_mean)
ncwriteatt([save_path, nc_Filename], 'asc_mean', 'Units', 'm/s')
ncwriteatt([save_path, nc_Filename], 'asc_mean', 'Description', ['Mean ascent rate'])
nccreate([save_path, nc_Filename], 'latitude', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'latitude', Sonde.lat)
ncwriteatt([save_path, nc_Filename], 'latitude', 'Units', 'Degrees N')
ncwriteatt([save_path, nc_Filename], 'latitude', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'latitude', 'Description', ['Radiosonde latitude'])
nccreate([save_path, nc_Filename], 'longitude', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'longitude', Sonde.lon)
ncwriteatt([save_path, nc_Filename], 'longitude', 'Units', 'Degrees E')
ncwriteatt([save_path, nc_Filename], 'longitude', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'longitude', 'Description', ['Radiosonde longitude'])
nccreate([save_path, nc_Filename], 'time', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'double', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'time', Time)
ncwriteatt([save_path, nc_Filename], 'time', 'Units', 'Seconds')
ncwriteatt([save_path, nc_Filename], 'time', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'time', 'Description', ['Radiosonde time (seconds since 1970-01-01)'])
nccreate([save_path, nc_Filename], 'p', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'p', Sonde.p)
ncwriteatt([save_path, nc_Filename], 'p', 'Units', 'hPa')
ncwriteatt([save_path, nc_Filename], 'p', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'p', 'Description', ['Pressure'])
nccreate([save_path, nc_Filename], 't', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 't', Sonde.T)
ncwriteatt([save_path, nc_Filename], 't', 'Units', 'Celsius')
ncwriteatt([save_path, nc_Filename], 't', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 't', 'Description', ['Temperature'])
nccreate([save_path, nc_Filename], 'q', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'q', Sonde.q)
ncwriteatt([save_path, nc_Filename], 'q', 'Units', 'g/kg')
ncwriteatt([save_path, nc_Filename], 'q', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'q', 'Description', ['Specific humidity'])
nccreate([save_path, nc_Filename], 'rh', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'rh', Sonde.RH)
ncwriteatt([save_path, nc_Filename], 'rh', 'Units', '%')
ncwriteatt([save_path, nc_Filename], 'rh', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'rh', 'Description', ['Relative humidity'])
nccreate([save_path, nc_Filename], 'rh_i', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'rh_i', Sonde.RH_i)
ncwriteatt([save_path, nc_Filename], 'rh_i', 'Units', '%')
ncwriteatt([save_path, nc_Filename], 'rh_i', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'rh_i', 'Description', ['Relative humidity with respect to ice; using Murphy and Koop (2005)'])
nccreate([save_path, nc_Filename], 'theta', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'theta', Sonde.theta)
ncwriteatt([save_path, nc_Filename], 'theta', 'Units', 'K')
ncwriteatt([save_path, nc_Filename], 'theta', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'theta', 'Description', ['Potential temperature'])
nccreate([save_path, nc_Filename], 'theta_e', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'theta_e', Sonde.theta_e)
ncwriteatt([save_path, nc_Filename], 'theta_e', 'Units', 'K')
ncwriteatt([save_path, nc_Filename], 'theta_e', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'theta_e', 'Description', ['Equivalent potential temperature'])
nccreate([save_path, nc_Filename], 'u', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'u', Sonde.u)
ncwriteatt([save_path, nc_Filename], 'u', 'Units', 'm/s')
ncwriteatt([save_path, nc_Filename], 'u', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'u', 'Description', ['u component of wind'])
nccreate([save_path, nc_Filename], 'v', 'dimensions', {'Num_layers', Num_layers, 'Num_profiles', Num_profiles}, 'datatype', 'single', 'format', 'netcdf4', 'DeflateLevel', 9)
ncwrite([save_path, nc_Filename], 'v', Sonde.v)
ncwriteatt([save_path, nc_Filename], 'v', 'Units', 'm/s')
ncwriteatt([save_path, nc_Filename], 'v', 'Missing_value', '-9999')
ncwriteatt([save_path, nc_Filename], 'v', 'Description', ['v component of wind'])
ncwriteatt([save_path, nc_Filename], '/', 'Site', Site_name)
ncwriteatt([save_path, nc_Filename], '/', 'Generated by', 'Israel Silber')
ncwriteatt([save_path, nc_Filename], '/', 'Date generated', datestr(now))

%--------------------------------------------------------------------------------------------
%% Generate radar-sounding database.
Site = 1; % 0 - nsa, 1 - awr
if Site == 0
    kazr_analysis = true; % true - analysis of kazr AND sounding, false - only sounding analysis.
    STR = 'nsaC1';
    Site_name = 'North Slope, Alaska';
    KAZR_path = {'/bear/s1/data/nsa/kazr/moment_tmp/md/', '/bear/s1/data/nsa/kazr/moments/'};  % for each mode.
    radar_fieldname = {'signal_to_noise_ratio_copol', 'reflectivity_copol', 'mean_doppler_velocity_copol', 'spectral_width_copol', 'reflectivity_xpol'};
    SNR_fieldname = 'signal_to_noise_ratio_copol';
elseif Site == 1
    kazr_analysis = true; % true - analysis of kazr AND sounding, false - only sounding analysis.
    STR = 'awrM1';
    Site_name = 'McMurdo Station, Antarctica';
    KAZR_path = {'/bear/s2/data/AWARE/kazr/moments_tmp/md/', '/bear/s2/data/AWARE/kazr/'};  % for each mode.
    radar_fieldname = {'signal_to_noise_ratio_copol', 'reflectivity_copol', 'mean_doppler_velocity_copol', 'spectral_width_copol', 'reflectivity_copol'}; % NO XPOL IN AWR DATA SO JUST GRABBING COPOL TO AVOID ERRORS.
end
if kazr_analysis;       nc_Filename = [STR, '_full_cld_analysis_radar_sounding.nc'];
else                    nc_Filename = [STR, '_full_cld_analysis_sounding.nc'];        end
Sounding_nc_filename = ['~/', STR,'_full_sounding.nc'];
save_path = ['~/'];
t_sounding = datenum(1970,1,1,0,0, double(ncread(Sounding_nc_filename, 'release_time')));
rng_sounding = double(ncread(Sounding_nc_filename, 'height'));
sounding_mean_asc = double(ncread(Sounding_nc_filename, 'asc_mean'));
sounding_mean_vert_res = double(ncread(Sounding_nc_filename, 'vert_res_mean'));
Connectivity = 8; % ONLY relevant for 2D matrix;
RH_uncertainty = 5; % radiosonde RH measurement uncertainty
RH_sat_threshold = 100 - RH_uncertainty; % RH value to determine saturation.
Num_layers = length(rng_sounding);
l_cloud_thickness_threshold = 30; % sounding liquid cloud threshold in m
Window_size = [0 15]; % first component is for time prior to sounding release, second component for time following sounding release.
if length(Window_size) == 1;        Window_size = [Window_size      Window_size];       end
md_bin_removal = 0; % number of the first bins to remove from the md data (due to potential leakage and low chances for adding real information at low level) - 12 at McMurdo, 0 for Barrow (data seems good and starts higher).
if Site == 1;       md_bin_removal = 12;        end
KAZR_snr_cloud_thresh = -16; % threshold for cloud determination.
C_fraction_thresh = 1.00; % window fraction radar threshold
Max_allowed_layers = 20;
if kazr_analysis
    for ii = 1: size(KAZR_path, 1)
        if ii == 1;
            if Site ~= 1
                Flist_ge = give_me_files_and_subfolders([STR(1:3), 'kazrcorge*nco*.nc'], KAZR_path{ii, 2});
                Flist_md = give_me_files_and_subfolders([STR(1:3), 'kazrcormd*nco*.nc'], KAZR_path{ii, 1});
            else
                Flist_ge = give_me_files_and_subfolders([STR(1:3), 'kazrge*.cdf'], KAZR_path{ii, 2});
                Flist_md = give_me_files_and_subfolders([STR(1:3), 'kazrmd*.cdf'], KAZR_path{ii, 1});
            end
        else
            if Site ~= 1
                Flist_ge = [Flist_ge; give_me_files_and_subfolders([STR(1:3), 'kazrcorge*nco*.nc'], KAZR_path{ii, 2})];
                Flist_md = [Flist_md; give_me_files_and_subfolders([STR(1:3), 'kazrcormd*nco*.nc'], KAZR_path{ii, 1})];
            else
                Flist_ge = [Flist_ge; give_me_files_and_subfolders([STR(1:3), 'kazrge*.cdf'], KAZR_path{ii, 2})];
                Flist_md = [Flist_md; give_me_files_and_subfolders([STR(1:3), 'kazrmd*.cdf'], KAZR_path{ii, 1})];
            end
            clear tmp_ge tmp_md
        end
    end
    %     return
    for ii = length(Flist_ge):-1:1
        if Site ~= 1
            ge_t(ii) = datenum(str2double(Flist_ge(ii).name(22:25)),str2double(Flist_ge(ii).name(26:27)),str2double(Flist_ge(ii).name(28:29)),str2double(Flist_ge(ii).name(31:32)),str2double(Flist_ge(ii).name(33:34)),str2double(Flist_ge(ii).name(35:36)));
        else
            ge_t(ii) = datenum(str2double(Flist_ge(ii).name(16:19)),str2double(Flist_ge(ii).name(20:21)),str2double(Flist_ge(ii).name(22:23)),str2double(Flist_ge(ii).name(25:26)),str2double(Flist_ge(ii).name(27:28)),str2double(Flist_ge(ii).name(29:30)));
        end
    end
    for ii = length(Flist_md):-1:1
        if Site ~= 1
            md_t(ii) = datenum(str2double(Flist_md(ii).name(22:25)),str2double(Flist_md(ii).name(26:27)),str2double(Flist_md(ii).name(28:29)),str2double(Flist_md(ii).name(31:32)),str2double(Flist_md(ii).name(33:34)),str2double(Flist_md(ii).name(35:36)));
        else
            md_t(ii) = datenum(str2double(Flist_md(ii).name(16:19)),str2double(Flist_md(ii).name(20:21)),str2double(Flist_md(ii).name(22:23)),str2double(Flist_md(ii).name(25:26)),str2double(Flist_md(ii).name(27:28)),str2double(Flist_md(ii).name(29:30)));
        end
    end
    Extra_loc = zeros(size(t_sounding));
    Extra_loc(round(t_sounding) - t_sounding < datenum(0,0,0,0,Window_size(2),0) & round(t_sounding) - t_sounding > 0) = 1; % end of day (load additional day)
    Extra_loc(t_sounding - round(t_sounding) < datenum(0,0,0,0,Window_size(1),0) & t_sounding - round(t_sounding) > 0) = 2; % begining of day (load additional day)
end
% Create output netCDF file and fields.
Control = {RH_uncertainty RH_sat_threshold l_cloud_thickness_threshold Window_size KAZR_snr_cloud_thresh C_fraction_thresh Max_allowed_layers Site_name};
if kazr_analysis
    MultiLayer_analysis_save_radar_nc(save_path, nc_Filename, Sounding_nc_filename, [], 1, [], Num_layers, Control, [], 1); % giving 1 to kazr because of this field is empty then the kazr data and analysis are not saved.
else
    MultiLayer_analysis_save_radar_nc(save_path, nc_Filename, Sounding_nc_filename, [], [], [], Num_layers, Control, [], 1); % Don't add KAZR fields
end
current_ind = 0; % index for writing the nc fields.
%
% Begin analysis
for ii =  1:length(t_sounding)
    
    % skip the period with bad KAZR data.
    if kazr_analysis
        if Site == 0; if (t_sounding(ii) > datenum(2012,12,16,5,0,0) && t_sounding(ii) < datenum(2013,3,17,18,0,0)) || (t_sounding(ii) > datenum(2017,9,25,23,0,0) && t_sounding(ii) < datenum(2017,9,26,0,0,0));    ...
                    disp('Bad kazr period (16-Dec-2012 05:30:00 to 17-Mar-2013 17:30:00. Skipping');  continue;       end;    end; % for nsa
        
        if Extra_loc(ii) == 0
            loc_ge = find(t_sounding(ii) - ge_t < 1 & t_sounding(ii) - ge_t >= 0);        loc_md = find(t_sounding(ii) - md_t < 1 & t_sounding(ii) - md_t >= 0);
        elseif Extra_loc(ii) == 1
            loc_ge = find(t_sounding(ii) - ge_t < 1 & t_sounding(ii) - ge_t >= -1);       loc_md = find(t_sounding(ii) - md_t < 1 & t_sounding(ii) - md_t >= -1);
        elseif Extra_loc(ii) == 2
            loc_ge = find(t_sounding(ii) - ge_t < 2 & t_sounding(ii) - ge_t >= 0);        loc_md = find(t_sounding(ii) - md_t < 2 & t_sounding(ii) - md_t >= 0);
        end
        if isempty(loc_ge) && isempty(loc_md);      continue;   else       disp(['Now processing ', datestr(t_sounding(ii))]);       end % checking for data.
        
        % load data (DEMANDING BOTH GE AND MD DATA).
        time_vec_md = []; % init time vector.
        time_vec_ge = []; % init time vector.
        SNR_data_md = [];
        SNR_data_ge = [];
        m0_data_ge = [];
        m0_data_md = [];
        m1_data_ge = [];
        m1_data_md = [];
        m2_data_ge = [];
        m2_data_md = [];
        m0_xpol_data_ge = [];
        m0_xpol_data_md = [];
        
        % GE processing
        if ~isempty(loc_ge)
            for jj = 1: length(loc_ge)
                time_tmp = datenum(0,0,0,0,0, ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], 'time')) + floor(ge_t(loc_ge(jj)));
%                 t_locs = find(abs(time_tmp - t_sounding(ii)) < datenum(0,0,0,0,Window_size,0));
                t_locs = find(time_tmp - t_sounding(ii) < datenum(0,0,0,0,Window_size(2),0)  &  time_tmp - t_sounding(ii) > -datenum(0,0,0,0,Window_size(1),0)); % ADDED 5/10/20: New assymetric window treatment.
                if isempty(t_locs);      continue;       end % checking for relevant data.
                alt_vec_ge = ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], 'range');
                r_locs = find(alt_vec_ge < rng_sounding(end));
                SNR_data_ge(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], radar_fieldname{1}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                m0_data_ge(1: length(r_locs), end + 1: end + length(t_locs)) = 10.^(ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], radar_fieldname{2}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)])./ 10); % assuming consecutive relevant indices for the count part.
                m1_data_ge(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], radar_fieldname{3}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                m2_data_ge(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], radar_fieldname{4}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                m0_xpol_data_ge(1: length(r_locs), end + 1: end + length(t_locs)) = 10.^(ncread([Flist_ge(loc_ge(jj)).path, Flist_ge(loc_ge(jj)).name], radar_fieldname{5}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)])./ 10); % assuming consecutive relevant indices for the count part.
                time_vec_ge(end + 1: end + length(t_locs)) = time_tmp(t_locs);
            end
            if size(m0_data_ge,2) == 1 || isempty(m0_data_ge) || all(isnan(m0_data_ge(:)));     continue;       end
            alt_vec_ge = double(alt_vec_ge(r_locs));
            m0_data_ge(m0_data_ge == -9999 | m0_data_ge == 0) = nan; % adding an option for 0 because of the 10 power calculation above.
            m1_data_ge(m1_data_ge == -9999) = nan;
            m2_data_ge(m2_data_ge == -9999) = nan;
            m0_xpol_data_ge(m0_xpol_data_ge == -9999 | m0_xpol_data_ge == 0) = nan; % adding an option for 0 because of the 10 power calculation above.
            [~,I] = unique(time_vec_ge); % remove redundancy.
            if length(I) < length(time_vec_ge)
                m0_data_ge = m0_data_ge(I);     m1_data_ge = m1_data_ge(I);     m2_data_ge = m2_data_ge(I);     m0_xpol_data_ge = m0_xpol_data_ge(I);     time_vec_ge = time_vec_ge(I);     SNR_data_ge = SNR_data_ge(I);     
            end
            
            % low-level GE artifact at Barrow on the first 2 range gates ( 3-range gates amplified region, but demanding at least 2 gates for cloud and 3rd is rather weakly amplified anyway...).
            Bad_data_highest_bin = 2;
            SNR_data_ge(1:Bad_data_highest_bin, :) = nan;
            m0_data_ge(1:Bad_data_highest_bin, :) = nan;
            m1_data_ge(1:Bad_data_highest_bin, :) = nan;
            m2_data_ge(1:Bad_data_highest_bin, :) = nan;
            m0_xpol_data_ge(1:Bad_data_highest_bin, :) = nan;
            Loc_no_kazr_data = find(rng_sounding < alt_vec_ge(Bad_data_highest_bin + 1), 1, 'last');
            
            % GE cloud
            ge_cloud = SNR_data_ge > KAZR_snr_cloud_thresh;
            Cloud_summary_padded = reshape([ge_cloud; zeros(size(ge_cloud))], size(ge_cloud, 1), 2* size(ge_cloud, 2));
            Cluster_mat = bwconncomp(Cloud_summary_padded, Connectivity);            Cluster_mat = Cluster_mat.PixelIdxList;
            Cluster_size = cell2mat(cellfun(@(x) length(x), Cluster_mat, 'UniformOutput', false)); % Finding size of each cluster.
            Cluster_ind = cell2mat(Cluster_mat(Cluster_size == 1)'); % Cropping only sporadic counts (i.e., thin non-persistent counts, i.e., single pixel noise).
            Cloud_summary_padded(Cluster_ind) = 0; % removing "noise detctions", i.e., min cloud thickness is 60 m.
            ge_cloud = Cloud_summary_padded(:, 1:2: end);
            clear Cloud_summary_padded Cluster_mat Cluster_size Cluster_ind jj
        end
        
        % MD processing
        if ~isempty(loc_md)
            for jj = 1: length(loc_md)
                time_tmp = datenum(0,0,0,0,0, ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], 'time')) + floor(md_t(loc_md(jj)));
%                 t_locs = find(abs(time_tmp - t_sounding(ii)) < datenum(0,0,0,0,Window_size,0));
                t_locs = find(time_tmp - t_sounding(ii) < datenum(0,0,0,0,Window_size(2),0)  &  time_tmp - t_sounding(ii) > -datenum(0,0,0,0,Window_size(1),0)); % ADDED 5/10/20: New assymetric window treatment.
                if isempty(t_locs);      continue;       end % checking for relevant data.
                alt_vec_md = ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], 'range');
                r_locs = find(alt_vec_md < rng_sounding(end));
                SNR_data_md(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], radar_fieldname{1}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                m0_data_md(1: length(r_locs), end + 1: end + length(t_locs)) = 10.^(ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], radar_fieldname{2}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)])./ 10); % assuming consecutive relevant indices for the count part.
                m1_data_md(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], radar_fieldname{3}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                m2_data_md(1: length(r_locs), end + 1: end + length(t_locs)) = ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], radar_fieldname{4}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)]); % assuming consecutive relevant indices for the count part.
                time_vec_md(end + 1: end + length(t_locs)) = time_tmp(t_locs);
                m0_xpol_data_md(1: length(r_locs), end + 1: end + length(t_locs)) = 10.^(ncread([Flist_md(loc_md(jj)).path, Flist_md(loc_md(jj)).name], radar_fieldname{5}, ...
                    [r_locs(1) t_locs(1)], [length(r_locs) length(t_locs)])./ 10); % assuming consecutive relevant indices for the count part.
            end
            if size(m0_data_md,2) == 1 || isempty(m0_data_md) || all(isnan(m0_data_md(:)));     continue;       end
            alt_vec_md = double(alt_vec_md(r_locs));
            m0_data_md(m0_data_md == -9999 | m0_data_md == 0) = nan; % adding an option for 0 because of the 10 power calculation above.
            m1_data_md(m1_data_md == -9999) = nan;
            m2_data_md(m2_data_md == -9999) = nan;
            m0_xpol_data_md(m0_xpol_data_md == -9999 | m0_xpol_data_md == 0) = nan; % adding an option for 0 because of the 10 power calculation above.
            [~,I] = unique(time_vec_md); % remove redundancy.
            if length(I) < length(time_vec_md)
                m0_data_md = m0_data_md(I);     m1_data_md = m1_data_md(I);     m2_data_md = m2_data_md(I);     m0_xpol_data_md = m0_xpol_data_md(I);     time_vec_md = time_vec_md(I);     SNR_data_md = SNR_data_md(I);     
            end
            
            % MD cloud   (and leakage and residuals removal).
            run Big_analysis_KAZR_HSRL_remove_md_leakage.m; % running the leakage removal process.
            if md_bin_removal > 0 % if there is something to remove.
                SNR_data_md(1:md_bin_removal, :) = nan;
                m0_data_md(1:md_bin_removal, :) = nan;
                m1_data_md(1:md_bin_removal, :) = nan;
                m2_data_md(1:md_bin_removal, :) = nan;
                m0_xpol_data_md(1:md_bin_removal, :) = nan;
            end
            %         md_eff_remove_alt = alt_vec_md(1) + md_alt_thresh; % conditional removal up to 'md_alt_thresh' above first md range gate.
            
            % finding highest KAZR cloud top (can generally use the cloud top relevant to the highest probed HSRL layer, but saving computation time needed for loops and gridding).
            md_cloud = SNR_data_md > KAZR_snr_cloud_thresh;
            %         md_cloud(sum(md_cloud, 2)./ size(md_cloud, 2) < cc_md_thresh & alt_vec_md < md_eff_remove_alt, :) = 0; % (For MCMURDO): removing low level low cloud fraction as these are the prime md leakage suspects for the full data window.
            Cloud_summary_padded = reshape([md_cloud; zeros(size(md_cloud))], size(md_cloud, 1), 2* size(md_cloud, 2));
            Cluster_mat = bwconncomp(Cloud_summary_padded, Connectivity);            Cluster_mat = Cluster_mat.PixelIdxList;
            Cluster_size = cell2mat(cellfun(@(x) length(x), Cluster_mat, 'UniformOutput', false)); % Finding size of each cluster.
            Cluster_ind = cell2mat(Cluster_mat(Cluster_size == 1)'); % Cropping only sporadic counts (i.e., thin non-persistent counts, i.e., single pixel noise).
            Cloud_summary_padded(Cluster_ind) = 0; % removing "noise detctions", i.e., min cloud thickness is 60 m.
            md_cloud = Cloud_summary_padded(:, 1:2: end);
            clear Cloud_summary_padded Cluster_mat Cluster_size Cluster_ind jj
        end
        
        if size(m0_data_ge,2) == 1 || isempty(m0_data_ge) || size(m0_data_md,2) == 1 || isempty(m0_data_md) || all(isnan(m0_data_ge(:))) || all(isnan(m0_data_md(:)));  continue;     end; % continue in case of bad data.
        current_ind = current_ind + 1;
        kazr.lowest_valid_kazr_ind = Loc_no_kazr_data + 1;
        % grid MD onto GE
        [x1,y1] = ndgrid(time_vec_md, alt_vec_md);
        [x2,y2] = ndgrid(time_vec_ge, alt_vec_ge);
        F1 = griddedInterpolant(x1, y1, md_cloud', 'linear', 'none'); % creates a gridded interpolant class which is orders of magnitude faster than griddata.
        cloud_mask_unite = any(cat(3, ge_cloud, F1(x2, y2)'), 3);  % ADDED 2/5/19 - GRID HSRL ONTO KAZR GE GRID
        kazr.frac = interp1(alt_vec_ge, nansum(cloud_mask_unite, 2)./ length(time_vec_ge), rng_sounding);       kazr.frac(1:Loc_no_kazr_data) = -9999;
        clear md_cloud ge_cloud F1
        F1 = griddedInterpolant(x1, y1, m0_data_md', 'linear', 'none'); % creates a gridded interpolant class which is orders of magnitude faster than griddata.
        m0_unite = nanmean(cat(3, m0_data_ge, F1(x2, y2)'), 3);  % ADDED 2/5/19 - GRID HSRL ONTO KAZR GE GRID
        m0_unite(~cloud_mask_unite) = nan;
        kazr.m0_mean = interp1(alt_vec_ge, 10.* log10(abs(nanmean(m0_unite, 2))), rng_sounding);       kazr.m0_mean(1:Loc_no_kazr_data) = -9999;
        kazr.m0_sd = interp1(alt_vec_ge, nanstd(10.* log10(abs(m0_unite)), 0, 2), rng_sounding);       kazr.m0_sd(1:Loc_no_kazr_data) = -9999; % SD of m0 in log scale.
        clear m0_data_md m0_data_ge m0_unite F1
        F1 = griddedInterpolant(x1, y1, m1_data_md', 'linear', 'none'); % creates a gridded interpolant class which is orders of magnitude faster than griddata.
        m1_unite = nanmean(cat(3, m1_data_ge, F1(x2, y2)'), 3);  % ADDED 2/5/19 - GRID HSRL ONTO KAZR GE GRID
        m1_unite(~cloud_mask_unite) = nan;
        kazr.m1_mean = interp1(alt_vec_ge, nanmean(m1_unite, 2), rng_sounding);       kazr.m1_mean(1:Loc_no_kazr_data) = -9999;
        kazr.m1_sd = interp1(alt_vec_ge, nanstd(m1_unite, 0, 2), rng_sounding);       kazr.m1_sd(1:Loc_no_kazr_data) = -9999;
        clear m1_data_md m1_data_ge m1_unite F1
        F1 = griddedInterpolant(x1, y1, m2_data_md', 'linear', 'none'); % creates a gridded interpolant class which is orders of magnitude faster than griddata.
        m2_unite = nanmean(cat(3, m2_data_ge, F1(x2, y2)'), 3);   % ADDED 2/5/19 - GRID HSRL ONTO KAZR GE GRID
        m2_unite(~cloud_mask_unite) = nan;
        kazr.m2_mean = interp1(alt_vec_ge, nanmean(m2_unite, 2), rng_sounding);       kazr.m2_mean(1:Loc_no_kazr_data) = -9999;
        kazr.m2_sd = interp1(alt_vec_ge, nanstd(m2_unite, 0, 2), rng_sounding);       kazr.m2_sd(1:Loc_no_kazr_data) = -9999;
        clear m2_data_md m2_data_ge m2_unite F1
        F1 = griddedInterpolant(x1, y1, m0_xpol_data_md', 'linear', 'none'); % creates a gridded interpolant class which is orders of magnitude faster than griddata.
        m0_xpol_unite = nanmean(cat(3, m0_xpol_data_ge, F1(x2, y2)'), 3);  % ADDED 2/5/19 - GRID HSRL ONTO KAZR GE GRID
        m0_xpol_unite(~cloud_mask_unite) = nan;
        kazr.m0_xpol_mean = interp1(alt_vec_ge, 10.* log10(abs(nanmean(m0_xpol_unite, 2))), rng_sounding);       kazr.m0_xpol_mean(1:Loc_no_kazr_data) = -9999;
        kazr.m0_xpol_sd = interp1(alt_vec_ge, nanstd(10.* log10(abs(m0_xpol_unite)), 0, 2), rng_sounding);       kazr.m0_xpol_sd(1:Loc_no_kazr_data) = -9999; % SD of m0 in log scale.
        clear m0_data_md m0_data_ge m0_unite F1 x1 x2 y1 y2 cloud_mask_unite
    else
        current_ind = current_ind + 1; % propagating writing counter.
    end
    
    % load sounding data.
    sonde.p = ncread(Sounding_nc_filename, 'p', [1 ii], [Num_layers, 1]);
    sonde.t = ncread(Sounding_nc_filename, 't', [1 ii], [Num_layers, 1]);
    sonde.q = ncread(Sounding_nc_filename, 'q', [1 ii], [Num_layers, 1]);
    sonde.rh = ncread(Sounding_nc_filename, 'rh', [1 ii], [Num_layers, 1]);
    sonde.rhi = ncread(Sounding_nc_filename, 'rh_i', [1 ii], [Num_layers, 1]);
    sonde.th = ncread(Sounding_nc_filename, 'theta', [1 ii], [Num_layers, 1]);
    sonde.the = ncread(Sounding_nc_filename, 'theta_e', [1 ii], [Num_layers, 1]);
    sonde.u = ncread(Sounding_nc_filename, 'u', [1 ii], [Num_layers, 1]);
    sonde.v = ncread(Sounding_nc_filename, 'v', [1 ii], [Num_layers, 1]);
    sonde.lat = ncread(Sounding_nc_filename, 'latitude', [1 ii], [Num_layers, 1]);
    sonde.lon = ncread(Sounding_nc_filename, 'longitude', [1 ii], [Num_layers, 1]);
    sonde.release_t = t_sounding(ii);
    sonde.asc_mean = sounding_mean_asc(ii);
    sonde.vert_res_mean = sounding_mean_vert_res(ii);
    
    % radar analysis
    if kazr_analysis
        Cluster_mat = bwconncomp(kazr.frac(:, end) >= C_fraction_thresh, Connectivity);
        Cluster_mat = Cluster_mat.PixelIdxList; % Getting the location of each of the cluster's elements.
        Cloud_size = cell2mat(cellfun(@(x) length(x), Cluster_mat, 'UniformOutput', false)); % Finding size of each cluster.
    else
        Cloud_size = []; % setting an empty array to make the 'if' below skip to the only-liquid conditions.
    end
    
    % sounding analysis
    sonde_res = diff(rng_sounding(1:2));
    Cluster_mat2 = bwconncomp(sonde.rh >= RH_sat_threshold, Connectivity);
    Cluster_mat2 = Cluster_mat2.PixelIdxList; % Getting the location of each of the cluster's elements.
    Cloud_size2 = cell2mat(cellfun(@(x) length(x), Cluster_mat2, 'UniformOutput', false)); % Finding size of each cluster.
    Cluster_mat2 = Cluster_mat2(Cloud_size2 >= l_cloud_thickness_threshold/sonde_res);
    Cloud_size2 = Cloud_size2(Cloud_size2 >= l_cloud_thickness_threshold/sonde_res);
    
    % summary analysis
    if kazr_analysis
        summary.c_num = int16(length(Cloud_size)); % number of cloud layers in profile
    end
    summary.l_num = int16(length(Cloud_size2)); % number of liquid cloud layers in profile
    summary.sat_rhi_mask = int16(sonde.rhi >= RH_sat_threshold);        summary.sat_rhi_mask(sonde.rhi == -9999) = -9999; % mask for saturated RH_i (same threshold as for liquid).
    if kazr_analysis
        summary.c_mask = int16(zeros(Num_layers, 1));                       summary.c_mask(1:Loc_no_kazr_data) = -9999;
        summary.l_mask = int16(zeros(Num_layers, 1));                       summary.l_mask(sonde.rhi == -9999) = -9999; % only based on sounding
    end
    summary.l_mask_unassigned = int16(zeros(Num_layers, 1));            summary.l_mask_unassigned(sonde.rhi == -9999) = -9999; % only based on sounding
    if kazr_analysis
        summary.c_h_bel_top = int16(zeros(Num_layers, 1));                  summary.c_h_bel_top(1:Loc_no_kazr_data) = -9999; % depth below top of cloud - could also be useful for liquid.
        summary.l_num_ic = int16(ones(Max_allowed_layers,1)) .* -9999;
        summary.l_num_ic_max = int16(zeros(1)); % max number of liquid layers per cloud.
        summary.l_ind = int16(ones(Max_allowed_layers,1)) .* -9999; % index for liquid layer in cloud layer.
        summary.c_thick = ones(Max_allowed_layers,1) .* -9999;
        summary.c_thick_max = 0;
        summary.c_base = ones(Max_allowed_layers,1) .* -9999;
        summary.c_top = ones(Max_allowed_layers,1) .* -9999;
        summary.l_is_ct = int16(ones(Max_allowed_layers,1)) .* -9999; % highest cloud bin contains liquid ( == 2 if more than one index below the lowest valid KAZR grid cell to avoid biases).
        summary.l_is_highest_l_ic = int16(ones(Max_allowed_layers,1)) .* -9999; % highest liquid layer in a cloud (will be one for pure liquid clouds (no KAZR signal). ( == 2 if more than one index below the lowest valid KAZR grid cell to avoid biases).
    end
    summary.l_h_bel_top = int16(zeros(Num_layers, 1));                  summary.l_h_bel_top(sonde.rhi == -9999) = -9999; % only based on sounding.
    summary.l_thick = ones(Max_allowed_layers,1) .* -9999;
    summary.l_thick_max = 0;
    summary.l_base = ones(Max_allowed_layers,1) .* -9999;
    summary.l_top = ones(Max_allowed_layers,1) .* -9999;
    if ~isempty(Cloud_size)
        match_l = zeros(1, length(Cloud_size2));
        for cc = 1: length(Cloud_size) % filling cloud mask.
            summary.c_mask(Cluster_mat{cc}) = cc;
            tot_num_layers = length(Cloud_size);
            if ~isempty(Cloud_size2)
                for ll = 1: length(Cloud_size2) % filling liquid mask (NOTE: matching numbers to the cloud array).
                    if match_l(ll) == 0;    M = ismember(Cluster_mat2{ll}, Cluster_mat{cc});        else       M = 0;       end % don't associate a cloud twice.
                    if any(M); % part of a detected cloud
                        summary.c_mask(Cluster_mat2{ll}) = cc;       summary.l_mask(Cluster_mat2{ll}) = cc;      match_l(ll) = cc; % match the cloud index and make sure full liquid overlap with the cloud mask.
                    end
                end
                if cc == length(Cloud_size)
                    if ~all(match_l > 0) % assign remaining liquid layers in the last step.
                        start_ind = max([max(summary.c_mask)  max(summary.l_mask)]) + 1;
                        for ll = find(match_l == 0);
                            summary.c_num = summary.c_num + 1;
                            match_l(ll) = start_ind;       start_ind = start_ind + 1;
                        end
                    end
                    % Complete the cloud properties now that the full mask
                    % is final.
                    tot_num_layers = max([length(Cloud_size) max(match_l)]);
                    summary.l_num_ic(1: tot_num_layers) = histcounts(match_l, (0: tot_num_layers) + 0.5);        summary.l_num_ic_max = max(summary.l_num_ic(1: tot_num_layers));
                    summary.l_is_highest_l_ic(1: length(Cloud_size2)) = 0; % first define all layers as non-highest
                    for llc = 1: tot_num_layers;    if summary.l_num_ic(llc) > 0;    llc_t = find(match_l == llc, 1, 'last');    summary.l_is_highest_l_ic(llc_t) = 1;       end;       end;        clear llc llc_t % now find highest per layer.
                    for ll = 1: length(Cloud_size2);
                        summary.l_mask_unassigned(Cluster_mat2{ll}) = ll;
                        summary.l_h_bel_top(Cluster_mat2{ll}) = (Cluster_mat2{ll}(end) - Cluster_mat2{ll}).* sonde_res;
                        summary.l_thick(ll) = Cloud_size2(ll) * sonde_res;
                        summary.l_base(ll) = rng_sounding(Cluster_mat2{ll}(1));          summary.l_top(ll) = rng_sounding(Cluster_mat2{ll}(end));
                        if match_l(ll) <= length(Cloud_size) % detected with KAZR
                            if Cluster_mat{match_l(ll)}(end) <= Cluster_mat2{ll}(end);       summary.l_is_ct(ll) = 1;
                            else                                                             summary.l_is_ct(ll) = 0;                              end
                        elseif Cluster_mat2{ll}(end) <= Loc_no_kazr_data;     % mark data not connected (from below) to the lowest valid KAZR data grid cell.
                            summary.l_is_ct(ll) = -1;                                   summary.l_is_highest_l_ic(ll) = -1;
                            summary.l_mask(Cluster_mat2{ll}, end) = -match_l(ll);       summary.c_mask(Cluster_mat2{ll}, end) = -match_l(ll);
                            summary.c_h_bel_top(Cluster_mat2{ll}) = -1;
                            summary.c_thick(match_l(ll)) = -1;                                        summary.c_base(match_l(ll)) = -1;                                    summary.c_top(match_l(ll)) = -1;
                            match_l(ll) = -match_l(ll);
                        else
                            summary.l_is_ct(ll) = 1;                                   summary.l_is_highest_l_ic(ll) = 1;
                            summary.l_mask(Cluster_mat2{ll}, end) = match_l(ll);       summary.c_mask(Cluster_mat2{ll}, end) = match_l(ll);
                            summary.c_h_bel_top(Cluster_mat2{ll}) = (Cluster_mat2{ll}(end) - Cluster_mat2{ll}).* sonde_res;
                            summary.c_thick(match_l(ll)) = Cloud_size2(ll) * sonde_res;
                            summary.c_base(match_l(ll)) = rng_sounding(Cluster_mat2{ll}(1));          summary.c_top(match_l(ll)) = rng_sounding(Cluster_mat2{ll}(end));
                        end
                        
                    end
                    summary.l_thick_max = max(summary.l_thick(1: length(Cloud_size2)));
                    summary.l_ind(1: length(Cloud_size2)) = match_l;
                end
            else
                summary.l_num_ic(1: length(Cloud_size)) = 0;
            end
        end
        for cc = 1: length(Cloud_size);    % rerunning for cloud thickness, as thickness was added due to liquid not detected by radar.
            c_indices = find(summary.c_mask == cc);
            summary.c_thick(cc) = length(c_indices) * sonde_res;
            summary.c_h_bel_top(c_indices) = (c_indices(end) - c_indices).* sonde_res;
            summary.c_base(cc) = rng_sounding(c_indices(1));          summary.c_top(cc) = rng_sounding(c_indices(end));
        end
        summary.c_thick_max = max(summary.c_thick(1: tot_num_layers));
    elseif ~isempty(Cloud_size2) % only a liquid cloud was detected or skip kazr analysis.
        for ll = 1: length(Cloud_size2);
            summary.l_mask_unassigned(Cluster_mat2{ll}) = ll;
            summary.l_h_bel_top(Cluster_mat2{ll}) = (Cluster_mat2{ll}(end) - Cluster_mat2{ll}).* sonde_res;
            summary.l_thick(ll) = Cloud_size2(ll) * sonde_res;
            summary.l_base(ll) = rng_sounding(Cluster_mat2{ll}(1));          summary.l_top(ll) = rng_sounding(Cluster_mat2{ll}(end));
            if kazr_analysis
                if Cluster_mat2{ll}(end) <= Loc_no_kazr_data;
                    summary.l_is_ct(ll) = -1;                                   summary.l_is_highest_l_ic(ll) = -1;
                    summary.l_mask(Cluster_mat2{ll}, end) = -ll;                summary.c_mask(Cluster_mat2{ll}, end) = -ll;
                    summary.c_h_bel_top(Cluster_mat2{ll}) = -1;
                    summary.c_thick(ll) = -1;                                   summary.c_base(ll) = -1;                                    summary.c_top(ll) = -1;
                    summary.l_ind(ll) = -ll;
                else
                    summary.l_is_ct(ll) = 1;                                   summary.l_is_highest_l_ic(ll) = 1;
                    summary.l_mask(Cluster_mat2{ll}, end) = ll;                summary.c_mask(Cluster_mat2{ll}, end) = ll;
                    summary.c_h_bel_top(Cluster_mat2{ll}) = (Cluster_mat2{ll}(end) - Cluster_mat2{ll}).* sonde_res;
                    summary.c_thick(ll) = Cloud_size2(ll) * sonde_res;
                    summary.c_base(ll) = rng_sounding(Cluster_mat2{ll}(1));          summary.c_top(ll) = rng_sounding(Cluster_mat2{ll}(end));
                    summary.l_ind(ll) = ll;
                end
            end
        end
        summary.l_thick_max = max(summary.l_thick(1: length(Cloud_size2)));
        if kazr_analysis
            summary.c_thick_max = max(summary.c_thick(1: length(Cloud_size2)));
            summary.l_num_ic(1: length(Cloud_size2)) = 1;                                   summary.l_num_ic_max = 1; % only liquid layers so a single layer per cloud.
        end
    end
    if kazr_analysis
        MultiLayer_analysis_save_radar_nc(save_path, nc_Filename, Sounding_nc_filename, sonde, kazr, summary, Num_layers, Control, current_ind, 0)
    else
        MultiLayer_analysis_save_radar_nc(save_path, nc_Filename, Sounding_nc_filename, sonde, [], summary, Num_layers, Control, current_ind, 0)
    end
end
clear ii start_ind cc ll match_l Cluster_mat Cluster_mat Cloud_size Cluster_mat2 Cluster_mat2 Cloud_size2
