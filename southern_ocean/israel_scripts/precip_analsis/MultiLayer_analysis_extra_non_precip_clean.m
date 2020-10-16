%% 8/31/2020
% Precipitation analysis and plotting of the sounding-KAZR dataset.
%---------------------------------------------------------------------------------------------------------------------------------
%%
clear all
%---------------------------------------------------------------------------------------------------------------------------------
%% set site (1 - McMurdo, 0 - NSA)
Site = 0;
%---------------------------------------------------------------------------------------------------------------------------------
%% Load data
if Site == 0;        nc_filename = ['/bear/s1/data/nsa/nsaC1_full_cld_analysis_radar_sounding_v3.nc']; % NSA
elseif Site == 1;    nc_filename = ['/bear/s1/data/nsa/awrM1_full_cld_analysis_radar_sounding_v3.nc']; % McMurdo
end
cl = double(ncread(nc_filename, 'liq_c_mask'));
ct = double(ncread(nc_filename, 'tot_c_mask'));
l_base = double(ncread(nc_filename, 'l_base'));
cf = double(ncread(nc_filename, 'win_c_fraction'));
cl_un = double(ncread(nc_filename, 'liq_c_mask_unassign'));     cl_un(cl_un <= 0) = 0;
h = double(ncread(nc_filename, 'height'));
t_sounding = double(ncread(nc_filename, 'sonde_time_datenum'));
m0 = double(ncread(nc_filename, 'win_c_m0_mean'));
m1 = double(ncread(nc_filename, 'win_c_m1_mean'));      m1(m1 == -9999 | ct == 0) = nan;
T = double(ncread(nc_filename, 't'));                   T(T == -9999) = nan;
Min_radar_height = h(ncread(nc_filename, 'lowest_valid_kazr_ind')); % get the vector
lwp = double(ncread(nc_filename, 'lwp_mean'));

%---------------------------------------------------------------------------------------------------------------------------------
%% Ananlysis of mixed layers.

% set control parameters
c_frac_thresh = 0.5; % cloud occurrence fraction threshold in 15-min window.
Min_radar_height_val = 300; % min radar range to inspect.
max_height = 4300; % max altitude to inspect (e..g., 4300 m to consider the altitude reached by a radiosonde 15 minutes after its release.
Radar_precip_depth_thresh = 60; % in m.
Min_radar_height(:) = 300; % simply use a constant value - show clutter free echoes.

disp(['Min radar height = ', num2str(Min_radar_height_val), ' m'])

if Site == 0; T(:, [6476]) = nan; end % cleaning bad sounding profile in the nsa dataset in v.2.6476 is with unreasonable LWP.
ct = cf > c_frac_thresh | cl_un > 0; % reproduce the cloud mask array in case of c_frac_thresh != 0.5
m0(m0 == -9999 | ct == 0) = nan; % remove m0 data if KAZR cloud fraction is below threshold


Radar_precip_depth_bin_thresh = Radar_precip_depth_thresh/ diff(h(1:2));

run ModelE_simulator_init_instrument.m % init instrument info and constants.
if Site == 0
    Z_min = Instrument_info.specs.nsa.Z_min_1km(2) + 20.* log10(h)-20.* log10(1e3);         Z_min(1) = nan;
elseif Site == 1
    Z_min = Instrument_info.specs.awr.Z_min_1km(2) + 20.* log10(h)-20.* log10(1e3);         Z_min(1) = nan;
end
Min_loc = find(h == Min_radar_height(1)) - 1;
cl(1:Min_loc, :) = 0;           cl_un(1:Min_loc + Radar_precip_depth_bin_thresh, :) = 0;
ct(1:Min_loc, :) = 0;                                           m0(1:Min_loc, :) =nan;
for ii = 1: size(cl_un, 2) % require subfreezing cloud temperatures.
    Cluster_mat = bwconncomp(cl_un(:,ii) > 0, 8);
    Cluster_mat = Cluster_mat.PixelIdxList; % Getting the location of each of the cluster's elements.
    Counter = 1;
    cl_un(:,ii) = 0;
    for jj = 1: length(Cluster_mat)
        Cond_a = all(T(Cluster_mat{jj}, ii) < 0) &&  all(T(Cluster_mat{jj}, ii) >= -40);
        if Cond_a
            cl_un(Cluster_mat{jj}, ii) = Counter;
            Counter = Counter + 1;
        end
    end
end
m0_bel_lcbh = [];
T_bel_lcbh = [];
T_cth = [];
lwp_per_case = [];
h_base = []; %
indices_high = []; % t indices.
l_per_prof = [];
h_depth = [];
IWflx_bel = [];
Free_locs =zeros(5, 0);
h_base_np = [];
for ii = 1: size(cl_un, 2)
    if max(cl_un(:, ii)) > 0
        layer_cc = max(cl_un(:, ii));
        for jj = 1: layer_cc
            Highest_l_loc = find(cl_un(:, ii) == jj);
            if h(Highest_l_loc(1)) >= max_height;       continue;       end % skip if the layer is above the height limit.
            Cond_a = all(T(Highest_l_loc, ii) < 0) &&  all(T(Highest_l_loc, ii) >= -40);
            Cond = 10.*log10(nansum(10.^(m0(Highest_l_loc(1)-Radar_precip_depth_bin_thresh: Highest_l_loc(1)-1, ii)./10))./Radar_precip_depth_bin_thresh) >= 10*log10(nanmean(10.^(Z_min(Highest_l_loc(1)-Radar_precip_depth_bin_thresh: Highest_l_loc(1)-1)./10))) & Cond_a;
            
            if ~Cond || Highest_l_loc(1) == 1;
                if Cond_a;
                    h_base_np(end+1) = h(Highest_l_loc(1));
                    D = datevec(t_sounding(ii));                    
                    Free_locs(1:5, end + 1) = [T(Highest_l_loc(end), ii), D(2)  sum(l_base(:, ii) > 0)   h(Highest_l_loc(end))-h(Highest_l_loc(1))     lwp(ii)];
                end;  continue;
            end
            
            m0_bel_lcbh(end + 1) =  10.*log10(nansum(10.^(m0(Highest_l_loc(1)-Radar_precip_depth_bin_thresh: Highest_l_loc(1)-1, ii)./10))./Radar_precip_depth_bin_thresh);
            T_bel_lcbh(end + 1) = T(Highest_l_loc(1), ii);
            T_cth(end + 1) = T(Highest_l_loc(end), ii);
            IWflx_bel(end + 1) = 10.^(2.42e-4.* (10.*log10(10.^(m0_bel_lcbh(end)./10))) .* T_bel_lcbh(end) + 6.99e-2.* (10.*log10(10.^(m0_bel_lcbh(end)./10))) - 1.86e-2.* T_bel_lcbh(end) - 1.63); % IWC based on Hogan et al., 2006, for 35 GHz.
            IWflx_bel(end) = IWflx_bel(end).*nanmin(m1(Highest_l_loc(1)-Radar_precip_depth_bin_thresh: Highest_l_loc(1)-1, ii)).*-3.6; % taking the most downward (negative) m1 in the examined range gates) and converting from g/(m^2*s) to mm/h.
            lwp_per_case(end + 1) = lwp(ii);
            h_base(end+1) = h(Highest_l_loc(1));
            h_depth(end + 1) = (Highest_l_loc(end) - Highest_l_loc(1) + 1) * diff(h(1:2));
            indices_high(end + 1) = ii; % index in the full dataset.
            l_per_prof(end + 1) = sum(l_base(:, ii) > 0); % number of layers per profile
        end
    end
end

%---------------------------------------------------------------------------------------------------------------------------------
%% Plot precipitation rate PDF
Fontsize = 14;
IWflx_bel(IWflx_bel <= 0) = nan; % removing updrafts (no IW flux).
R_bins = -6:0.5:1;
Cond = 1: length(IWflx_bel); % no constraints.
Fig = figure('position', [1400 0 800 700]); hold on;
Hist_type = 'pdf';
R2 = IWflx_bel.*0.1; % Rain rate IWC from Hogan et al. (2006) times m1 (after I divided by 1e3 to convert g/m^2 to mm, and converting from s to h, and changing sign).
a = histcounts(log10(R2(Cond)), R_bins, 'normalization', Hist_type); plot(10.^(R_bins(1:end-1) + diff(R_bins(1:2))/2), a, ':', 'linewidth', 4, 'color', rgb('black')); set(gca, 'xscale', 'log'); xlabel('SR [mm/h]'); ylabel('Occurrence frequency'); grid on; box on;% title('Ice Precipitation Rates');
R2 = IWflx_bel.*3; % Rain rate IWC from Hogan et al. (2006) times m1 (after I divided by 1e3 to convert g/m^2 to mm, and converting from s to h, and changing sign).
a = histcounts(log10(R2(Cond)), R_bins, 'normalization', Hist_type); plot(10.^(R_bins(1:end-1) + diff(R_bins(1:2))/2), a, ':', 'linewidth', 4, 'color', rgb('black')); set(gca, 'xscale', 'log'); xlabel('SR [mm/h]'); ylabel('Occurrence frequency'); grid on; box on; %title('Ice Precipitation Rates');
R = IWflx_bel; % Rain rate IWC from Hogan et al. (2006) times m1 (after I divided by 1e3 to convert g/m^2 to mm, and converting from s to h, and changing sign).
a = histcounts(log10(R(Cond)), R_bins, 'normalization', Hist_type); plot(10.^(R_bins(1:end-1) + diff(R_bins(1:2))/2), a, '--', 'linewidth', 4, 'color', rgb('black')); set(gca, 'xscale', 'log'); xlabel('SR [mm/h]'); ylabel('Occurrence frequency'); grid on; box on;% title('Ice Precipitation Rates');
R2 = IWflx_bel; % Rain rate IWC from Hogan et al. (2006) times m1 (after I divided by 1e3 to convert g/m^2 to mm, and converting from s to h, and changing sign).
hold on; a = histcounts(log10(R2(Cond)), R_bins, 'normalization', 'pdf'); fill(10.^([R_bins(1:end-1) R_bins(1)] + diff(R_bins(1:2))/2), [a 0], rgb('skyblue'), 'facealpha', 0.3, 'linewidth', 4); set(gca, 'xscale', 'log'); xlabel('R [mm/h]'); ylabel('PDF'); grid on; box on; %title('Ice Precipitation Rates'); % removing PDF units for paper revisions.

set(gca, 'fontsize', Fontsize, 'fontweight', 'bold'); xlabel('R [mm/h]'); grid on; box on; set(gca, 'xtick', 10.^(-6:1), 'xminorgrid', 'on', 'xminortick', 'on', 'xlim', 10.^([-6, 1]), 'ytick', 0:0.05:0.45, 'yticklabel', {'0'; num2str((0.05:0.05:0.45)', '%.2f')}, 'ylim', [0 0.40+eps]);

%---------------------------------------------------------------------------------------------------------------------------------
%% Histograms
Dtmp = datevec(t_sounding(indices_high));       Dtmp = Dtmp(:,2)';
H_bins = [0:250:8e3]./1e3; figure; HH = axes('position', [0.1380    0.1520    0.7400    0.7600]); hold on; grid; box on; histogram([h_base_np h_base]./1e3, H_bins, 'normalization', 'probability', 'facecolor', rgb('RosyBrown')); hold on;  histogram(h_base_np./1e3, H_bins, 'normalization', 'probability', 'facecolor', rgb('seagreen')); xlabel('Cloud top altitude [km]'); ylabel('Occurrence frequency'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold')
H_bins = 0.5:12.5; figure; HH = axes('position', [0.1380    0.1520    0.7400    0.7600]); hold on; grid; box on; histogram([Free_locs(2,:) Dtmp], H_bins, 'normalization', 'probability', 'facecolor', rgb('RosyBrown')); hold on;  histogram(Free_locs(2,:), H_bins, 'normalization', 'probability', 'facecolor', rgb('seagreen')); legend('All supercooled', 'Non-precipitating', 'location', 'best'); legend('boxoff'); xlabel('Month'); ylabel('Occurrence frequency'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold', 'xtick', 1:12, 'xticklabel', datestr(datenum(1,1:12,1), 'mmm')); axis tight;
H_bins = -40:2:0; figure; HH = axes('position', [0.1380    0.1520    0.7400    0.7600]); hold on; grid; box on; histogram([Free_locs(1,Free_locs(3,:) >= 1) T_cth(l_per_prof >= 1)], H_bins, 'normalization', 'probability', 'facecolor', rgb('RosyBrown')); hold on;  histogram(Free_locs(1,Free_locs(3,:) >= 1), H_bins, 'normalization', 'probability', 'facecolor', rgb('seagreen'));  xlabel('Cloud top temperature [^oC]'); ylabel('Occurrence frequency'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold'); ylim([0 0.15]); set(HH, 'ytick', 0:0.03:0.15)
a = histcounts(T_cth(l_per_prof >= 1), H_bins, 'normalization', 'count') ./ histcounts([Free_locs(1,Free_locs(3,:) >= 1) T_cth(l_per_prof >= 1)], H_bins, 'normalization', 'count'); H = axes('position', get(gca, 'position'), 'color', 'none', 'yaxislocation', 'right');  hold on; plot(H_bins(1:end-1) + diff(H_bins)./2, a, 'k', 'linewidth', 3); ylabel('Precipitating fraction'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold', 'xlim', HH.XLim, 'ylim', [0 1], 'xticklabel', []);
H_bins = 0:50:1.5e3; figure; HH = axes('position', [0.1380    0.1520    0.7400    0.7600]); hold on; grid; box on; histogram([Free_locs(4,Free_locs(3,:) >= 1) h_depth(l_per_prof >= 1)], H_bins, 'normalization', 'probability', 'facecolor', rgb('RosyBrown')); hold on;  histogram(Free_locs(4,Free_locs(3,:) >= 1), H_bins, 'normalization', 'probability', 'facecolor', rgb('seagreen'));  xlabel('Cloud depth [m]'); ylabel('Occurrence frequency'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold')
H_bins = 0:25:300; figure; HH = axes('position', [0.1380    0.1520    0.7400    0.7600]); hold on; grid; box on; histogram([Free_locs(5,Free_locs(3,:) == 1) lwp_per_case(l_per_prof == 1)], H_bins, 'normalization', 'probability', 'facecolor', rgb('RosyBrown')); hold on;  histogram(Free_locs(5,Free_locs(3,:) == 1), H_bins, 'normalization', 'probability', 'facecolor', rgb('seagreen')); xlabel('LWP [g/m^2]'); ylabel('Occurrence frequency'); set(gca, 'fontsize', Fontsize, 'fontweight', 'bold')
return
%---------------------------------------------------------------------------------------------------------------------------------
%% Precipitation detection sensitivity emulation (supercooled cloud)
Min_radar_height = 300;
if Site == 0; STR_SITE = 'nsa';
elseif Site == 1; STR_SITE = 'awr'; end
cl = double(ncread(nc_filename, 'liq_c_mask'));
ct = double(ncread(nc_filename, 'tot_c_mask'));
cf = double(ncread(nc_filename, 'win_c_fraction'));
cl_un = double(ncread(nc_filename, 'liq_c_mask_unassign'));
h = double(ncread(nc_filename, 'height'));
l_base = double(ncread(nc_filename, 'l_base'));
m0 = double(ncread(nc_filename, 'win_c_m0_mean'));
t = double(ncread(nc_filename, 't'));       t(t == -9999) = nan;
if Site == 0; t(:, [7026]) = nan;       end
Radar_precip_depth_thresh_init = 60; % in m.
Radar_precip_depth_thresh = 60; % in m.
Radar_precip_depth_bin_thresh = Radar_precip_depth_thresh/ diff(h(1:2));
Min_loc = find(h == Min_radar_height) - 1;
cl(1:Min_loc + Radar_precip_depth_bin_thresh, :) = 0;           cl_un(1:Min_loc + Radar_precip_depth_bin_thresh, :) = 0;
ct(1:Min_loc, :) = 0;                                           m0(1:Min_loc, :) =nan;
cf(1:Min_loc, :) = 0;
Max_z_loop = 66;
Max_depth_loop = 37;
Ice_above_liq_thresh = 5; % grid cell threshold. Ice must be below this number of grid cells, in order to consider the liquid as cloud top

% Allocate analysis summary arrays.
Counter_lowest = zeros(Max_depth_loop, Max_z_loop, size(cl, 2));
Counter_unshielded_highest = zeros(Max_depth_loop, Max_z_loop, size(cl, 2));
Counter_highest = zeros(Max_depth_loop, Max_z_loop, size(cl, 2));
Counter_per_layer = zeros(Max_depth_loop, Max_z_loop, size(cl, 2));
Counter_per_layer_kazr = zeros(Max_depth_loop, size(cl, 2));
Counter_lowest_kazr = zeros(Max_depth_loop, size(cl, 2));
Counter_highest_kazr = zeros(Max_depth_loop, size(cl, 2));
Counter_unshielded_highest_kazr = zeros(Max_depth_loop, size(cl, 2));
Counter_profiles = zeros(Max_depth_loop, size(cl, 2));
Counter_unshielded_profiles = zeros(Max_depth_loop, size(cl, 2));

cl_un(cl_un <= 0) = nan;
Counter_40minus = 0;
max_height = 4300; % max altitude to inspect (e..g., 4300 m to consider the altitude reached by a radiosonde 15 minutes after its release.
ct = cf > 0.5 | cl_un > 0; % 5/9/20  - changed the required hourly fraction to account for LCBH drift (over the NSA - 250 m mean, 300 m 3rd quartile. - no need to rerun the full Multi-layer processing
m0(m0 == -9999 | ct == 0) = nan; % remove m0 data if KAZR cloud fraction is below 0.25
run ModelE_simulator_init_instrument.m % init instrument info and constants.

Z_min = Instrument_info.specs.(STR_SITE).Z_min_1km(2) + 20.* log10(h)-20.* log10(1e3);         Z_min(1) = nan; % minimum detectable radar echo (based on a separated in-depth analysis)
too_high = zeros(1, size(cl_un, 2));
disp(['min radar height = ', num2str(Min_radar_height), ' m'])
for ii = 1: size(cl_un, 2) % require subfreezing cloud temperatures.
    Cluster_mat = bwconncomp(cl_un(:,ii) > 0, 8);
    Cluster_mat = Cluster_mat.PixelIdxList; % Getting the location of each of the cluster's elements.
    Counter = 1;
    cl_un(:, ii) = 0;
    for jj = 1: length(Cluster_mat)
        if all(t(Cluster_mat{jj}, ii) < 0) && ~any(t(Cluster_mat{jj}, ii) < -40)
            cl_un(Cluster_mat{jj}, ii) = Counter;
            Counter = Counter + 1;
        else
            Counter_40minus = Counter_40minus + 1;
        end
    end
    too_high(ii) = sum(l_base(:, ii) >= max_height);
end
Counter_profiles(1, :) = max(cl_un) - min(cl_un) + 1 - too_high;     Counter_profiles(1, isnan(Counter_profiles(1, :))) = 0; % number of layers above threshold height
Cloudsat_m0_thresh_init = -50; % detection threshold not considering attenuation in dBZ. see https://cloudsat.atmos.colostate.edu/instrument

for dd = 1:Max_depth_loop
    
    % remove data below liquid for detection
    disp(['Loop #', num2str(dd)])
    Radar_precip_depth_thresh = Radar_precip_depth_thresh_init + (dd - 1)*15; % in m.
    Radar_precip_depth_bin_thresh = Radar_precip_depth_thresh/ diff(h(1:2));
    cl_un(1:Min_loc + Radar_precip_depth_bin_thresh, :) = 0;
    cl_un(cl_un == 0) = nan;        Counter_profiles(dd, :) = nanmax(cl_un) - nanmin(cl_un) + 1 - too_high;     Counter_profiles(dd, isnan(Counter_profiles(dd, :)) | nanmax(cl_un) == 0) = 0; % number of layers above threshold height (NOTE: ASSUMING THAT Min_loc + Radar_precip_depth_bin_thresh < max_height
    
    for zz = 1: Max_z_loop
        Cloudsat_m0_thresh = Cloudsat_m0_thresh_init + zz - 1;
        for ii = 1: size(cl_un,2)
            if Counter_profiles(dd, ii) > 0
                Cluster_mat2 = bwconncomp(cl_un(:,ii) > 0, 8);
                Cluster_mat2 = Cluster_mat2.PixelIdxList; % Getting the location of each of the cluster's elements.
                for ll = 1: length(Cluster_mat2)
                    if zz == 1 % kazr data
                        if ll == 1; Counter_unshielded_profiles(dd, ii) = find(ct(:,ii) > 0, 1, 'last') - Cluster_mat2{length(Cluster_mat2)}(end) < Ice_above_liq_thresh & h(Cluster_mat2{length(Cluster_mat2)}(1)) < max_height;       end % UPDATED 5/14/20 - ADDED MAX HEIGHT.
                        Cond = 10.*log10(nansum(10.^(m0(Cluster_mat2{ll}(1)-Radar_precip_depth_bin_thresh: Cluster_mat2{ll}(1)-1, ii)./10))./Radar_precip_depth_bin_thresh) >= 10*log10(nanmean(10.^(Z_min(Cluster_mat2{ll}(1)-Radar_precip_depth_bin_thresh: Cluster_mat2{ll}(1)-1)./10))) & h(Cluster_mat2{ll}(1)) < max_height; % KAZR detectiblity% UPDATED 5/14/20 - ADDED MAX HEIGHT.
                        if Cond
                            Counter_per_layer_kazr(dd, ii) = Counter_per_layer_kazr(dd, ii) + 1;
                            if ll == length(Cluster_mat2); Counter_highest_kazr(dd, ii) = 1;   if Counter_unshielded_profiles(dd, ii) == 1;  Counter_unshielded_highest_kazr(dd, ii) = 1;   end;   end% highest non-shielded layer.
                        end
                        if 10 * log10(nansum(10.^(m0(Min_loc+1: Min_loc+Radar_precip_depth_bin_thresh, ii)./10))./Radar_precip_depth_bin_thresh) >= 10*log10(nanmean(10.^(Z_min(Cluster_mat2{ll}(1)-Radar_precip_depth_bin_thresh: Cluster_mat2{ll}(1)-1)./10)));
                            Counter_lowest_kazr(dd, ii) = 1;
                        end; % checking for cloud in the lowest bin
                    end
                    Cond = 10.*log10(nansum(10.^(m0(Cluster_mat2{ll}(1)-Radar_precip_depth_bin_thresh: Cluster_mat2{ll}(1)-1, ii)./10))./Radar_precip_depth_bin_thresh) >= Cloudsat_m0_thresh & h(Cluster_mat2{ll}(1)) < max_height;
                    if Cond;
                        Counter_per_layer(dd, zz, ii) = Counter_per_layer(dd, zz, ii) + 1;
                        if ll == length(Cluster_mat2); Counter_highest(dd, zz, ii) = 1;   if Counter_unshielded_profiles(dd, ii) == 1;  Counter_unshielded_highest(dd, zz, ii) = 1;   end;    end% highest non-shielded layer.
                    end
                    if 10 * log10(nansum(10.^(m0(Min_loc+1: Min_loc+Radar_precip_depth_bin_thresh, ii)./10))./Radar_precip_depth_bin_thresh) >= Cloudsat_m0_thresh
                        Counter_lowest(dd, zz, ii) = 1;
                    end; % checking for cloud in the lowest bin
                end
            end
        end
    end
end
if Site ~= 1
    save(['~/', STR_SITE, 'C1_precip_supercooled_no40_new_analysis_radar_sounding_v3_0.50_4.3km_min_h_', num2str(Min_radar_height, '%04.0f'), '.mat'], 'Counter_per_layer_kazr', 'Counter_lowest_kazr', 'Counter_profiles', 'Counter_per_layer', 'Counter_lowest', 'Max_z_loop', 'Max_depth_loop', ...
        'Counter_40minus', 'Min_radar_height', 'Radar_precip_depth_thresh_init', 'Cloudsat_m0_thresh_init', 'Counter_highest_kazr', 'Counter_highest', 'Counter_unshielded_profiles', 'Counter_unshielded_highest_kazr', 'Counter_unshielded_highest');
else
    save(['~/', STR_SITE, 'M1_precip_supercooled_no40_new_analysis_radar_sounding_v3_0.50_4.3km_min_h_', num2str(Min_radar_height, '%04.0f'), '.mat'], 'Counter_per_layer_kazr', 'Counter_lowest_kazr', 'Counter_profiles', 'Counter_per_layer', 'Counter_lowest', 'Max_z_loop', 'Max_depth_loop', ...
        'Counter_40minus', 'Min_radar_height', 'Radar_precip_depth_thresh_init', 'Cloudsat_m0_thresh_init', 'Counter_highest_kazr', 'Counter_highest', 'Counter_unshielded_profiles', 'Counter_unshielded_highest_kazr', 'Counter_unshielded_highest');
end

%---------------------------------------------------------------------------------------------------------------------------------
%% Plot precipitation sensitivity analysis results
Min_radar_height_2use = 300;
    load('jet_modified.mat');       load('HSV_modified2.mat');
    if Site == 1;
        load(['~/awrM1_precip_supercooled_no40_new_analysis_radar_sounding_v3_0.50_4.3km_min_h_', num2str(Min_radar_height_2use, '%04.0f'),'.mat'])
    elseif Site == 0;
        load(['~/nsaC1_precip_supercooled_no40_new_analysis_radar_sounding_v3_0.50_4.3km_min_h_', num2str(Min_radar_height_2use, '%04.0f'),'.mat'])
    end
    % load('/bear/s1/data/nsa/nsaC1_precip_warm_new_analysis_radar_sounding.mat')
    Fontsize = 26;
    Plot_data = {sum(Counter_per_layer > 0,3)./  repmat(sum(Counter_profiles > 0, 2), 1,Max_z_loop).*100, ...
        sum(Counter_per_layer,3) ./  repmat(sum(Counter_profiles, 2), 1,Max_z_loop).*100, ...
        sum(Counter_lowest,3) ./  repmat(sum(Counter_profiles > 0, 2), 1,Max_z_loop).*100, ...
        sum(Counter_highest,3) ./  repmat(sum(Counter_profiles > 0, 2), 1,Max_z_loop).*100, ...
        sum(Counter_unshielded_highest,3) ./  repmat(sum(Counter_unshielded_profiles > 0, 2), 1,Max_z_loop).*100};
    Title = {'Precipitation from liquid clouds per profile',...
        'Precipitation from liquid clouds per detected liquid cloud',...
        ['Surface precipitation (', num2str(Min_radar_height_2use), ' m AGL) from liquid clouds'],...
        'Precipitation from the topmost layer', ...
        'Precipitation from the topmost unshielded layer'};
    
for ii = [2 5 3]
    Fig = figure('position', [1400 0 1100 1000], 'visible', 'on');     hold on;  axis square;        box  on;        grid on;     shading flat %     grid minor;
    scatter(60, -50, 300, 'r', 'filled', 'markeredgecolor', 'k')
    scatter(480, -29, 300, 'vg', 'filled', 'markeredgecolor', 'k')
    scatter(240, -29, 300, '^g', 'filled', 'markeredgecolor', 'k')
    scatter(240, -15, 300, 'sb', 'filled', 'markeredgecolor', 'k')
    scatter(240, -5, 300, 'sm', 'filled', 'markeredgecolor', 'k')
    scatter(250, 12, 300, 'dc', 'filled', 'markeredgecolor', 'k')
    Plot_data{ii}(1,1)
    pcolor(Radar_precip_depth_thresh_init + (0:Max_depth_loop - 1).*15, Cloudsat_m0_thresh_init + (0:Max_z_loop-1), Plot_data{ii}');
    shading flat
    grid on;
    caxis([0 100]);     xlabel('Vertical resolution [m]');        ylabel('Reflectivity threshold [dBZ]');
    axis tight;
    set(gca, 'fontsize', Fontsize, 'fontweight', 'bold', 'position', [0.1300    0.1100    0.7200    0.8150]);
    colormap(jet_modified)
    if ii == 2 && Site == 0; legend('KAZR', 'CPR', 'CPR (oversampling)', '2C-PC & 2C-SP (possible)', '2C-PC (certain)', 'KaPR (GPM)'); end
    [h1, h2] = contour(Radar_precip_depth_thresh_init + (0:Max_depth_loop - 1).*15, Cloudsat_m0_thresh_init + (0:Max_z_loop-1), Plot_data{ii}', 0:10:100, 'k', 'showtext', 'on', 'linewidth', 2);
    clabel(h1, h2, 'fontweight', 'bold', 'fontsize', 20)
    scatter(60, -50, 300, 'r', 'filled', 'markeredgecolor', 'k')
    scatter(480, -29, 300, 'vg', 'filled', 'markeredgecolor', 'k')
    scatter(240, -29, 300, '^g', 'filled', 'markeredgecolor', 'k')
    scatter(240, -15, 300, 'sb', 'filled', 'markeredgecolor', 'k')
    scatter(240, -5, 300, 'sm', 'filled', 'markeredgecolor', 'k')
    scatter(250, 12, 300, 'dc', 'filled', 'markeredgecolor', 'k')
end