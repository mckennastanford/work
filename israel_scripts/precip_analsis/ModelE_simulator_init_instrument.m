%-------------------------------------------------------------------------
% Israel Silber
% Last update: 09/04/2019
% Last update: 04/19/2020 (added (M)WACR and CPR, updated Ze_min for KAZR).
% Last update: 06/09/2020 (added CEIL and RL.
%-------------------------------------------------------------------------
% Init instrument (radar/lidar) information (see all the comments for more
% information).
%-------------------------------------------------------------------------
% NOTE: either frequency or wavelength of the instruments must be provided.
%-------------------------------------------------------------------------

Instrument_info.c = 299792458; % speed of light in m.
R_d = 287.058; % gas constant for dry air [J/(kg*K)].

% init instrument data.
Instrument_info.intrument_str = {'HSRL', 'KAZR', '1064nm', 'WACR', 'CPR', 'RL', 'CEIL', 'XSACR'}; % WACR is either for WACR or MWACR.
Instrument_info.is_radar = [false true false true true false false true]; % true - radar, false - lidar.
Instrument_info.freq = [nan 34.860 nan 95.04 94.05 nan nan 9.710]; % in GHz (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
Instrument_info.lambda = [0.532 nan 1.064 nan nan 0.355 0.910 nan]; % in microns
Instrument_info.ext_OD = [4 nan 4 nan nan 4 4 nan]; % OD value where we get full extinction of the lidar signal (should be ~3-5 by default).
Instrument_info.K_w = [nan 0.88 nan 0.84 0.75 nan nan 0.93]; % index of refraction for water used for Ze calculation |K_w|^2 - see the ARM KAZR handbook (Widener et al.,2012), WACR handbook (Widener and Johnson, 2006), and Stephens et al., 2008
Instrument_info.eps_liq = [complex(1.337273,1.7570744e-9)^2 complex(5.489262,2.8267679)^2 complex(1.320416,1.2588968e-6)^2 complex(3.468221,2.1423486)^2 complex(3.468221,2.1423486)^2 complex(1.357247,2.4198595e-9)^2 complex(1.323434,5.6988883e-7)^2 complex(8.112180,1.7811075)^2]; % Complex dielectric constant for liquid water from refractive indices taken from Segelstein 1981.
Instrument_info.specs.sgp.pt = [nan 2000 nan 1513 1820 nan nan nan]; % transmitting power in watts (see prosensing webpage at: https://www.prosensing.com/crb-product/ka-band-zenith-radar-kazr/, Mead, 2010, WACR Calibration Study", and Table 1 in Stephens et al., 2008 for CPR)
Instrument_info.specs.nsa.pt = [nan 2000 nan nan 1820 nan nan nan]; % transmitting power in watts (see prosensing webpage at: https://www.prosensing.com/crb-product/ka-band-zenith-radar-kazr/, Widener and Johnson, 2006 (for MWACR), and Table 1 in Stephens et al., 2008 for CPR)
Instrument_info.specs.awr.pt = [nan 2000 nan nan 1820 nan nan nan]; % transmitting power in watts (see prosensing webpage at: https://www.prosensing.com/crb-product/ka-band-zenith-radar-kazr/
Instrument_info.specs.ena.pt = [nan 2200 nan nan 1820 nan nan nan]; % transmitting power in watts (see: https://www.atmos-meas-tech.net/12/4931/2019/amt-12-4931-2019.pdf/ and Stephens et al., 2008 for CPR)
% Instrument_info.specs.sgp.pt = [nan 150 nan nan nan]; % transmitting power in watts (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
% Instrument_info.specs.nsa.pt = [nan 150 nan nan nan]; % transmitting power in watts (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
% Instrument_info.specs.awr.pt = [nan 150 nan nan nan]; % transmitting power in watts (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
Instrument_info.specs.sgp.theta = [nan 0.19 nan 0.19 0.12 nan nan 1.40]; % 3 dB beamwidth in degrees (see \muN), Table 17-2 in Kollias et al. (2016) https://journals.ametsoc.org/doi/full/10.1175/AMSMONOGRAPHS-D-15-0037.1
Instrument_info.specs.nsa.theta = [nan 0.31 nan nan 0.12 nan nan 1.40]; % 3 dB beamwidth in degrees (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
Instrument_info.specs.awr.theta = [nan 0.31 nan 0.38 0.12 nan nan 1.40]; % 3 dB beamwidth in degrees (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY and Kneifel et al., JGR, 2015 for MWACR)
Instrument_info.specs.ena.theta = [nan 0.31 nan nan 0.12 nan nan nan]; % 3 dB beamwidth in degrees (see: https://www.atmos-meas-tech.net/12/4931/2019/amt-12-4931-2019.pdf/)
Instrument_info.specs.sgp.gain = [nan 10^(57.48/10) nan 10^(39.4/10) 10^(63.1/10) nan nan 10^(42.0/10)]; % gain (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
Instrument_info.specs.nsa.gain = [nan 10^(53.37/10) nan nan 10^(63.1/10) nan nan 10^(42.0/10)]; % gain (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY)
Instrument_info.specs.awr.gain = [nan 10^(52.73/10) nan 10^(37.8/10) 10^(63.1/10) nan nan 10^(42.0/10)]; % gain (see Bharadwaj et al., KA-BAND ARM ZENITH PROFILING RADAR NETWORK FOR CLIMATE STUDY) - NO INFORMATION YET FOR MWACR
Instrument_info.specs.ena.gain = [nan nan nan nan 10^(63.1/10) nan nan nan]; % gain - NO INFORMATION YET FOR MWACR AND KAZR2
Instrument_info.specs.sgp.Z_min_1km = [nan -51.5 nan nan -30 nan nan nan]; % Z_min in dBZ at 1 km based on 8-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold), Table 1 in Stephens et al., 2008 for CPR - NOT CHARACTERIZED FOR SGP WACR.
% Instrument_info.specs.nsa.Z_min_1km = [nan -49.0 nan nan nan]; % Z_min in dBZ at 1 km based on 8-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold).
% Instrument_info.specs.awr.Z_min_1km = [nan -47.2 nan nan nan]; % Z_min in dBZ at 1 km based on 1-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold).
Instrument_info.specs.nsa.Z_min_1km = [nan -48.5 nan -46.0 -30 nan nan nan]; % Z_min in dBZ at 1 km based on 8-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold), Table 1 in Stephens et al., 2008 for CPR.
Instrument_info.specs.awr.Z_min_1km = [nan -45.5 nan -40.0 -30 nan nan nan]; % Z_min in dBZ at 1 km based on 1-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold), Table 1 in Stephens et al., 2008 for CPR - NOT CHARACTERIZED FOR AWR MWACR.
Instrument_info.specs.ena.Z_min_1km = [nan -56.5 nan nan -30 nan nan nan]; % Z_min in dBZ at 1 km based on 1-year analysis of data (minimum detected hydrometeor signal using -16 dB SNR threshold), Table 1 in Stephens et al., 2008 for CPR.
Instrument_info.specs.nsa.lr = [nan 10^(4/10) nan nan nan nan nan nan]; % attenuation based on the general attributes in the spectra files
Instrument_info.specs.nsa.pr_noise_ge = [nan 10^(-68.5/10) nan nan nan nan nan nan]; % minimum detectable signal in mW based on the general attributes in the spectra files
Instrument_info.specs.nsa.pr_noise_md = [nan 10^(-72.3/10) nan nan nan nan nan nan]; % minimum detectable signal in mW based on the general attributes in the spectra files
Instrument_info.specs.sgp.tau_ge = [nan nan nan nan 0.3 nan nan nan nan]; % pulse width in mus based on Tab;e 17-1 in Kollias et al. (2016) https://journals.ametsoc.org/doi/full/10.1175/AMSMONOGRAPHS-D-15-0037.1
Instrument_info.specs.nsa.tau_ge = [nan 0.3 nan nan nan 3.3 nan nan nan]; % pulse width in mus based on the general attributes in the spectra files and Table 1 in Stephens et al., 2008 for CPR
Instrument_info.specs.nsa.tau_md = [nan 4.0 nan nan nan 3.3 nan nan nan]; % pulse width in mus based on the general attributes in the spectra files and Table 1 in Stephens et al., 2008 for CPR
Instrument_info.specs.ena.tau_ge = [nan 0.2 nan nan nan 3.3 nan nan nan]; % pulse width in mus based on the general attributes in the spectra files and Table 1 in Stephens et al., 2008 for CPR
Instrument_info.specs.ena.tau_md = [nan 4.0 nan nan nan 3.3 nan nan nan]; % pulse width in mus based on the general attributes in the spectra files and Table 1 in Stephens et al., 2008 for CPR

% calculate lambda/f based on f/lambda
for instinst = 1: length(Instrument_info.freq)
    if isnan(Instrument_info.freq(instinst))
        Instrument_info.freq(instinst) = Instrument_info.c / (Instrument_info.lambda(instinst) * 1e-6) * 1e-9;
    end
    if isnan(Instrument_info.lambda(instinst))
        Instrument_info.lambda(instinst) = Instrument_info.c / (Instrument_info.freq(instinst) * 1e9) * 1e6;
    end
end
clear instinst