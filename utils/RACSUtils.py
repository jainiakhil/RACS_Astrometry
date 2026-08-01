from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob
import os
import gc

def CleanRACS(df_racs, racs_raw):
# ## Cleaning the RACS catalogue as per our requirements:
#     remove sources with neighbouring sources < 30 arcsec
#     remove non-point sources (int/peak < 1.5)
#     snr > 8

    racs_final = []

    for i in range(len(df_racs)):
        racs_coord = SkyCoord(racs_raw[i]['ra'], racs_raw[i]['dec'], unit=u.degree)
        
        try:
            _, sep, _ = racs_coord.match_to_catalog_sky(racs_coord, nthneighbor=2)
            compactness = racs_raw[i]['total_flux_gaussian'] / racs_raw[i]['peak_flux']
            snr = racs_raw[i]['peak_flux'] / racs_raw[i]['e_peak_flux']
            
            racs = racs_raw[i][(sep > 30*u.arcsec) & (compactness < 1.5) & (snr > 6)]
        
        except:
            racs = racs_raw[i].copy()
        
        racs_final.append(racs)
        
    return racs_final 


def CleanRACS2(df_racs, racs_raw):
# ## Cleaning the RACS catalogue as per our requirements:
#     remove sources with neighbouring sources < 5 arcsec
#     remove non-point sources (int/peak < 1.5)
#     snr > 6

    racs_final = []

    for i in range(len(df_racs)):
        racs_coord = SkyCoord(racs_raw[i]['ra'], racs_raw[i]['dec'], unit=u.degree)
        
        try:
            _, sep, _ = racs_coord.match_to_catalog_sky(racs_coord, nthneighbor=2)
            compactness = racs_raw[i]['total_flux_gaussian'] / racs_raw[i]['peak_flux']
            snr = racs_raw[i]['peak_flux'] / racs_raw[i]['e_peak_flux']
            
            racs = racs_raw[i][(sep > 5*u.arcsec) & (compactness < 1.5) & (snr > 6)]
        
        except:
            racs = racs_raw[i].copy()
        
        racs_final.append(racs)
        
    return racs_final 


def CleanRACSLow3(df_racs, racs_raw):
# ## Cleaning the RACS catalogue as per our requirements:
#     remove sources with neighbouring sources < 30 arcsec
#     remove non-point sources (int/peak < 1.5)
#     snr > 6

    racs_final = []

    for i in range(len(df_racs)):
        racs_coord = SkyCoord(racs_raw[i]['col_ra_deg_cont'], racs_raw[i]['col_dec_deg_cont'], unit=u.degree)
        
        try:
            _, sep, _ = racs_coord.match_to_catalog_sky(racs_coord, nthneighbor=2)
            compactness = racs_raw[i]['col_flux_int'] / racs_raw[i]['col_flux_peak']
            snr = racs_raw[i]['col_flux_peak'] / racs_raw[i]['col_flux_peak_error']
            
            racs = racs_raw[i][(sep > 30*u.arcsec) & (compactness < 1.5) & (snr > 6)]
        
        except:
            racs = racs_raw[i].copy()
        
        racs_final.append(racs)
        
    return racs_final 


def CleanRACSMid(df_racs, racs_raw):
# ## Cleaning the RACS catalogue as per our requirements:
#     remove sources with neighbouring sources < 5 arcsec
#     remove non-point sources (int/peak < 1.5)
#     snr > 6

    racs_final = []

    for i in range(len(df_racs)):
        racs_coord = SkyCoord(racs_raw[i]['ra'], racs_raw[i]['dec'], unit=u.degree)
        
        try:
            _, sep, _ = racs_coord.match_to_catalog_sky(racs_coord, nthneighbor=2)
            compactness = racs_raw[i]['int_flux'] / racs_raw[i]['peak_flux']
            snr = racs_raw[i]['peak_flux'] / racs_raw[i]['err_peak_flux']
            
            racs = racs_raw[i][(sep > 5*u.arcsec) & (compactness < 1.5) & (snr > 6)]
        
        except:
            racs = racs_raw[i].copy()
        
        racs_final.append(racs)
        
    return racs_final 


def CleanRACSHigh(df_racs, racs_raw):
# ## Cleaning the RACS catalogue as per our requirements:
#     remove sources with neighbouring sources < 30 arcsec
#     remove non-point sources (int/peak < 1.5)
#     snr > 6

    racs_final = []

    for i in range(len(df_racs)):
        racs_coord = SkyCoord(racs_raw[i]['col_ra_deg_cont'], racs_raw[i]['col_dec_deg_cont'], unit=u.degree)
        
        try:
            _, sep, _ = racs_coord.match_to_catalog_sky(racs_coord, nthneighbor=2)
            compactness = racs_raw[i]['col_flux_int'] / racs_raw[i]['col_flux_peak']
            snr = racs_raw[i]['col_flux_peak'] / racs_raw[i]['col_flux_peak_err']
            
            racs = racs_raw[i][(sep > 30*u.arcsec) & (compactness < 1.5) & (snr > 6)]
        
        except:
            racs = racs_raw[i].copy()
        
        racs_final.append(racs)
        
    return racs_final 


def CleanWISE(df_racs, wise):
## Cleaning the WISE catalogue and calculate RA and DEC errors as per our requirements    
    
    wise_final = []
    
    for i in range(len(df_racs)):
        wise_coord = SkyCoord(wise[i]['RAJ2000'], wise[i]['DEJ2000'], unit=u.degree)
        
        try:
            _, sep, _ = wise_coord.match_to_catalog_sky(wise_coord, nthneighbor=2)
            ra_err = ((wise[i]['eeMaj']**2) * np.sin(np.deg2rad(wise[i]['eePA']))**2) + ((wise[i]['eeMin']**2) * np.cos(np.deg2rad(wise[i]['eePA']))**2)
            dec_err = ((wise[i]['eeMaj']**2) * np.cos(np.deg2rad(wise[i]['eePA']))**2) + ((wise[i]['eeMin']**2) * np.sin(np.deg2rad(wise[i]['eePA']))**2)
            
            wise1 = wise[i].copy()
            wise1['RA_ERR'] = ra_err
            wise1['DEC_ERR'] = dec_err
            wise1 = wise1[(sep > 5*u.arcsec)]
        
        except:
            wise1 = wise[i].copy()
        
        wise_final.append(wise1)
    
    return wise_final


def CleanFIRST(df_racs, first):
# ## Cleaning the FIRST catalogue and calculate RA and DEC errors as per our requirements:

    first_final = []

    for i in range(len(df_racs)):
        try:
            snr = (first[i]['Fpeak'] - 0.25) / first[i]['Rms']
            maj_err = first[i]['fMaj'] * (snr**-1 + 20**-1) / 1.645
            min_err = first[i]['fMin'] * (snr**-1 + 20**-1) / 1.645
            
            # ra_err = ((maj_err**2) * np.sin(np.deg2rad(first3d[k][i]['fPA']))**2) + ((min_err**2) * np.cos(np.deg2rad(first3d[k][i]['fPA']))**2)
            # dec_err = ((maj_err**2) * np.cos(np.deg2rad(first3d[k][i]['fPA']))**2) + ((min_err**2) * np.sin(np.deg2rad(first3d[k][i]['fPA']))**2)
            
            first1 = first[i].copy()
            first1['RA_ERR'] = maj_err
            first1['DEC_ERR'] = min_err
        
        except:
            first1 = first[i].copy()
        
        # first1 = first3d[k][i][(compactness < 2) & (snr > 6)]
        
        first_final.append(first1)
    
    return first_final


def compare_survey(target_coord, reference_coord, 
                target_ra_err=None, target_dec_err=None, reference_ra_err=None, reference_dec_err=None, 
                survey='', radius=5*u.arcsec, floor=0, ax=None, beam_num=None):
    
    idx, sep, _ = target_coord.match_to_catalog_sky(reference_coord)

    ind1 = sep < radius
    
    target_match_coord = target_coord[ind1]
    reference_match_coord = reference_coord[idx][ind1]
    
    tot_sources = np.unique(idx[ind1]).shape[0]
    
    dra, ddec = target_match_coord.spherical_offsets_to(reference_match_coord)
    dra, ddec = dra.arcsec, ddec.arcsec
        
    if target_ra_err is not None and target_dec_err is not None and reference_ra_err is not None and reference_dec_err is not None:
        xerr = np.sqrt(target_ra_err[ind1]**2 + reference_ra_err[idx][ind1]**2)
        yerr = np.sqrt(target_dec_err[ind1]**2 + reference_dec_err[idx][ind1]**2)
        
        # Floor condition: If the error value is less than floor value, set it to floor value
        # xerr[xerr < floor] = floor
        # yerr[yerr < floor] = floor
        xerr = np.sqrt(xerr**2 + floor**2)
        yerr = np.sqrt(yerr**2 + floor**2)
        
        ax.errorbar(dra, ddec, xerr=xerr, yerr=yerr, fmt='o', alpha=0.4, c='black')
        dra_mean = np.average(dra, weights=1/xerr**2)
        ddec_mean = np.average(ddec, weights=1/yerr**2)
        
        dra_uncert = np.sqrt(1 / np.sum(1/xerr**2))
        ddec_uncert = np.sqrt(1 / np.sum(1/yerr**2))

    else:
        ax.scatter(dra, ddec, marker='o', alpha=0.4, c='black')
        dra_mean, ddec_mean = dra.mean(), ddec.mean()
        dra_rms, ddec_rms = dra.std(), ddec.std()
        
    ax.axvline(dra_mean, color='red', ls='--', label=r'$\bar\Delta RA={:.2f}\pm{:.3f}$"'.format(dra_mean, dra_uncert))
    ax.axhline(ddec_mean, color='red', ls='--', label=r'$\bar\Delta DEC={:.2f}\pm{:.3f}$"'.format(ddec_mean, ddec_uncert))

    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$DEC (arcsec)')
    ax.set_title(f"{survey} - Beam {beam_num}")
    ax.text(0.01, 0.99, f"Total sources compared in field: {tot_sources}", bbox=dict(facecolor='white', alpha=0.3), 
            ha='left', va='top', transform=ax.transAxes)

    ax.legend()

    return dra_mean, ddec_mean, dra_uncert, ddec_uncert, tot_sources


def compare_survey_noplot(target_coord, reference_coord, 
                        target_ra_err=None, target_dec_err=None, reference_ra_err=None, reference_dec_err=None, 
                        radius=5*u.arcsec, floor=0):
    
    idx, sep, _ = target_coord.match_to_catalog_sky(reference_coord)

    ind1 = sep < radius
    
    target_match_coord = target_coord[ind1]
    reference_match_coord = reference_coord[idx][ind1]
    
    tot_sources = np.unique(idx[ind1]).shape[0]
    
    dra, ddec = target_match_coord.spherical_offsets_to(reference_match_coord)
    dra, ddec = dra.arcsec, ddec.arcsec

    if target_ra_err is not None and target_dec_err is not None and reference_ra_err is not None and reference_dec_err is not None:
        xerr = np.sqrt(target_ra_err[ind1]**2 + reference_ra_err[idx][ind1]**2)
        yerr = np.sqrt(target_dec_err[ind1]**2 + reference_dec_err[idx][ind1]**2)
        
        # Floor condition: If the error value is less than floor value, set it to floor value
        # xerr[xerr < floor] = floor
        # yerr[yerr < floor] = floor
        xerr = np.sqrt(xerr**2 + floor**2)
        yerr = np.sqrt(yerr**2 + floor**2)
            
        dra_mean = np.average(dra, weights=1/xerr**2)
        ddec_mean = np.average(ddec, weights=1/yerr**2)
        
        dra_uncert = np.sqrt(1 / np.sum(1/xerr**2))
        ddec_uncert = np.sqrt(1 / np.sum(1/yerr**2))

    else:
        dra_mean, ddec_mean = dra.mean(), ddec.mean()
        dra_rms, ddec_rms = dra.std(), ddec.std()
    

    return dra_mean, ddec_mean, dra_uncert, ddec_uncert, tot_sources


def objective_function(coefficients, offset_obs_vals, offset_obs_uncert):
    offset_beam, offset_scan = coefficients[0:36], coefficients[36:]
    
    offset_modelled = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_modelled
    offset_residual = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_residual
    
    for i in range(np.size(offset_scan)):
        for j in range(np.size(offset_beam)):
            offset_modelled[i, j] = offset_beam[j] + offset_scan[i]
            
            if offset_obs_vals[i, j] == 0:
                offset_modelled[i, j] = 0
    
    offset_residual = offset_obs_vals - offset_modelled
    
    return np.sum((offset_residual/offset_obs_uncert)**2)  # minimize the chi-squared value


def generate_offset_model(coefficients, offset_obs_vals):
    offset_beam, offset_scan = coefficients[0:36], coefficients[36:]
    
    offset_modelled = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_modelled
    offset_residual = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_residual
    
    for i in range(np.size(offset_scan)):
        for j in range(np.size(offset_beam)):
            offset_modelled[i, j] = offset_beam[j] + offset_scan[i]
            
            # if offset_obs_vals[i, j] == 0:
            #     offset_modelled[i, j] = 0
                
    offset_residual = offset_obs_vals - offset_modelled
    
    return offset_beam, offset_scan, offset_modelled, offset_residual


def objective_function_bad_regions(coefficients, offset_beam, offset_obs_vals, offset_obs_uncert):
    # offset_beam, offset_scan = coefficients[0:36], coefficients[36:]
    offset_scan = coefficients
    
    offset_modelled = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_modelled
    offset_residual = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_residual
    
    for i in range(np.size(offset_scan)):
        for j in range(np.size(offset_beam)):
            offset_modelled[i, j] = offset_beam[j] + offset_scan[i]
            
            if offset_obs_vals[i, j] == 0:
                offset_modelled[i, j] = 0
    
    offset_residual = offset_obs_vals - offset_modelled
    
    return np.sum((offset_residual/offset_obs_uncert)**2)  # minimize the chi-squared value


def generate_offset_model_bad_regions(coefficients, offset_beam, offset_obs_vals):
    offset_scan = coefficients
    
    offset_modelled = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_modelled
    offset_residual = np.zeros((np.size(offset_scan), np.size(offset_beam)))  # Initialize offset_residual
    
    for i in range(np.size(offset_scan)):
        for j in range(np.size(offset_beam)):
            offset_modelled[i, j] = offset_beam[j] + offset_scan[i]
            
            # if offset_obs_vals[i, j] == 0:
            #     offset_modelled[i, j] = 0
                
    offset_residual = offset_obs_vals - offset_modelled
    
    return offset_beam, offset_scan, offset_modelled, offset_residual


def plot_beam_offset(directory, df_racs, 
                    offset_scan_ra, offset_scan_dec, offset_beam_ra, offset_beam_dec, dra_plot3d_cal, ddec_plot3d_cal, offset_residual_ra, offset_residual_dec, 
                    chi_squared_ra_observed, chi_squared_dec_observed, chi_squared_observed, chi_squared_ra_residual, chi_squared_dec_residual, chi_squared_residual, 
                    radius, cal_sbids, directory_in=''):
    
    # Plot the beam-wise mean offsets for each scan
    os.makedirs(os.path.join(directory, directory_in), exist_ok=True)

    ii = 0

    for k in range(len(df_racs)): 
        field_name = df_racs[k].iloc[0]['FIELD_NAME'][-7:]
        utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START']).replace(':', '-').replace(' ', '_').split('.')[0]
        sbid_val = df_racs[k].iloc[0]['SBID']
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']

        if str(cal_sbid_val) == cal_sbids:
            print(f"Plotting for Scan {k} having SBID: {sbid_val} and CAL_SBID: {cal_sbid_val}")
            
            offset_obs_beam_avg = np.sqrt(offset_scan_ra[ii]**2 + offset_scan_dec[ii]**2)
            offset_obs_scan_avg = np.sqrt(offset_beam_ra**2 + offset_beam_dec**2)
            offset_obs = np.sqrt(np.array(dra_plot3d_cal[ii])**2 + np.array(ddec_plot3d_cal[ii])**2)
            offset_res = np.sqrt(offset_residual_ra[ii]**2 + offset_residual_dec[ii]**2)
            
            plt.rcParams.update({'font.size': 10})
            
            fig, bx = plt.subplots(2, 2, figsize=(20, 18), dpi=300)
            
            norm = matplotlib.colors.Normalize(vmin=0, vmax=1.5)
            cm = matplotlib.cm.copper

            sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
            sm.set_array([])
            
            bx[0, 0].quiver(df_racs[k]['RA_DEG'].mean(), df_racs[k]['DEC_DEG'].mean(), offset_scan_ra[ii], offset_scan_dec[ii], 
                    scale=1, scale_units='xy', angles='xy', color=cm(norm(offset_obs_beam_avg)), alpha=0.8)
            # print(f"0,0: {offset_scan_m1_ra[ii]} {offset_scan_m1_dec[ii]}")
            # print(f"Magnitude: {np.sqrt(offset_scan_m1_ra[ii]**2 + offset_scan_m1_dec[ii]**2)} Angle: {np.arctan(offset_scan_m1_dec[ii]/offset_scan_m1_ra[ii])}")
            
            # add labels
            bx[0, 0].text(df_racs[k]['RA_DEG'].mean()-0.15, df_racs[k]['DEC_DEG'].mean()-0.15, f"{offset_obs_beam_avg:.2f}", fontsize=8)
            # for _, row in df_racs[k].iterrows():
            #     bx[0, 0].text(row['RA_DEG'], row['DEC_DEG'], int(row['BEAM_NUM']), fontsize=8)

            bx[0, 0].set_title(f"Beam-independent Scan Correction")
            bx[0, 0].set(xlabel="RA (deg)", ylabel="DEC (deg)", xlim=(min(df_racs[k]['RA_DEG'])-0.1, max(df_racs[k]['RA_DEG'])+0.1), 
                    ylim=(min(df_racs[k]['DEC_DEG'])-0.1, max(df_racs[k]['DEC_DEG'])+0.1))
            
            bx[0, 1].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], offset_beam_ra, offset_beam_dec, 
                    scale=1, scale_units='xy', angles='xy', color=cm(norm(offset_obs_scan_avg)), alpha=0.8)
            # print(f"0,1: {offset_beam_m1_ra} {offset_beam_m1_dec}")
            # print(f"Magnitude: {np.sqrt(offset_beam_m1_ra**2 + offset_beam_m1_dec**2)} Angle: {np.arctan(offset_beam_m1_dec/offset_beam_m1_ra)}")
            
            # add labels
            for zz, row in df_racs[k].iterrows():
                    bx[0, 1].text(row['RA_DEG'], row['DEC_DEG'], int(row['BEAM_NUM']), fontsize=8)
                    bx[0, 1].text(row['RA_DEG'], row['DEC_DEG']-0.15, f"{offset_obs_scan_avg[zz]:.2f}", fontsize=8)

            bx[0, 1].set_title(f"Scan-independent Beam Correction")
            bx[0, 1].set(xlabel="RA (deg)", ylabel="DEC (deg)")

            bx[1, 0].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], dra_plot3d_cal[ii], ddec_plot3d_cal[ii], 
                    scale=1, scale_units='xy', angles='xy', color=cm(norm(offset_obs)), alpha=0.8)
            bx[1, 0].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], dra_plot3d_cal[ii], np.zeros_like(dra_plot3d_cal[ii]),
                    scale=1, scale_units='xy', angles='xy', color='black', alpha=0.16)
            bx[1, 0].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], np.zeros_like(ddec_plot3d_cal[ii]), ddec_plot3d_cal[ii],
                    scale=1, scale_units='xy', angles='xy', color='black', alpha=0.16)
            clb = plt.colorbar(sm, ax=bx, format='%.2f')  # Specify the ax argument to steal space from bx
            # print(f"1,0: {np.asarray(dra_plot3d_190424[ii])} {np.asarray(ddec_plot3d_190424[ii])}")
            # print(f"Magnitude: {np.sqrt(np.asarray(dra_plot3d_190424[ii])**2 + np.asarray(ddec_plot3d_190424[ii])**2)} Angle: {np.arctan(np.asarray(ddec_plot3d_190424[ii])/np.asarray(dra_plot3d_190424[ii]))}")
            
            # add labels
            for zz, row in df_racs[k].iterrows():
                    bx[1, 0].text(row['RA_DEG'], row['DEC_DEG'], int(row['BEAM_NUM']), fontsize=8)
                    bx[1, 0].text(row['RA_DEG'], row['DEC_DEG']-0.15, f"{offset_obs[zz]:.2f}", fontsize=8)

            bx[1, 0].set_title(f"Mean Observed\n$\chi^2$: For RA: {np.sum(chi_squared_ra_observed[ii]):.2f}, For DEC: {np.sum(chi_squared_dec_observed[ii]):.2f}, Total: {chi_squared_observed[ii]:.2f}")
            bx[1, 0].set(xlabel="RA (deg)", ylabel="DEC (deg)")
            # bx[0].ylabel("DEC (deg)")
            clb.set_label("Offset Magnitude (in arcsec)")
            
            bx[1, 1].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], offset_residual_ra[ii], offset_residual_dec[ii], 
                    scale=1, scale_units='xy', angles='xy', color=cm(norm(offset_res)), alpha=0.8)
            bx[1, 1].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], offset_residual_ra[ii], np.zeros_like(offset_residual_ra[ii]),
                    scale=1, scale_units='xy', angles='xy', color='black', alpha=0.16)
            bx[1, 1].quiver(df_racs[k]['RA_DEG'], df_racs[k]['DEC_DEG'], np.zeros_like(offset_residual_dec[ii]), offset_residual_dec[ii],
                    scale=1, scale_units='xy', angles='xy', color='black', alpha=0.16)
            # clb = plt.colorbar(sm, ax=bx, format='%.2f')  # Specify the ax argument to steal space from bx
            # print(f"1,1: {offset_residual_m1_ra[ii]} {offset_residual_m1_dec[ii]}")
            # print(f"Magnitude: {np.sqrt(offset_residual_m1_ra[ii]**2 + offset_residual_m1_dec[ii]**2)} Angle: {np.arctan(offset_residual_m1_dec[ii]/offset_residual_m1_ra[ii])}")
            
            # add labels
            for zz, row in df_racs[k].iterrows():
                    bx[1, 1].text(row['RA_DEG'], row['DEC_DEG'], int(row['BEAM_NUM']), fontsize=8)
                    bx[1, 1].text(row['RA_DEG'], row['DEC_DEG']-0.15, f"{offset_res[zz]:.2f}", fontsize=8)
            
            bx[1, 1].set_title(f"Residual\n$\chi^2$: For RA: {np.sum(chi_squared_ra_residual[ii]):.2f}, For DEC: {np.sum(chi_squared_dec_residual[ii]):.2f}, Total: {chi_squared_residual[ii]:.2f}")
            bx[1, 1].set(xlabel="RA (deg)", ylabel="DEC (deg)")
            # bx[1].xlabel("RA (deg)")
            # bx[1].ylabel("DEC (deg)")
            clb.set_label("Offset Magnitude (in arcsec)")
            
            fig.suptitle(f"{ii} Model3 {utc_time[0:10]} Beam Offsets. SBID: {sbid_val} Field: {field_name}", fontsize=25)
            # fig.text(0.435, 0.465, f"$\Delta\chi^2$: \nFor RA: {(chi_squared_ra_observed[ii] - chi_squared_ra_residual[ii]):.2f}, \nFor DEC: {(chi_squared_dec_observed[ii] - chi_squared_dec_residual[ii]):.2f}, \nTotal: {(chi_squared_observed[ii] - chi_squared_residual[ii]):.2f}", ha='center', fontsize=12)

            fig.savefig(f"{directory}/{directory_in}/{utc_time}_RACSLow_vs_WISE_Offsets_SBID_{sbid_val}_Field_{field_name}_rad{radius}_model3_v5.png", dpi=300, bbox_inches='tight')
            fig.clf()
            # plt.close(fig)
            # plt.close('all')
            # gc.collect()
            
            ii += 1
                    
    return None


def CorrectRACS(df_racs, full_offset_modelled_ra, full_offset_modelled_dec, dra_uncert_plot3d, ddec_uncert_plot3d, cal_sbids, radius, 
                directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME'][-8:]
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            load_racs_query = os.path.join(directory_racs, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            
            racs = [Table(np.load(os.path.join(load_racs_query, f'Beam_{i}.npy'))) for i in range(len(os.listdir(load_racs_query))) 
                    if os.path.isfile(os.path.join(load_racs_query, f'Beam_{i}.npy'))]
            racs_final = CleanRACS(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_new = racs_final[i]['ra']*u.deg + full_offset_modelled_ra[k][i]*u.arcsec
                dec_new = racs_final[i]['dec']*u.deg + full_offset_modelled_dec[k][i]*u.arcsec
                e_ra_new = np.sqrt((racs_final[i]['e_ra']*u.arcsec)**2 + (dra_uncert_plot3d[k][i]*u.arcsec)**2) if dra_uncert_plot3d[k][i] != np.inf else racs_final[i]['e_ra']*u.arcsec
                e_dec_new = np.sqrt((racs_final[i]['e_dec']*u.arcsec)**2 + (ddec_uncert_plot3d[k][i]*u.arcsec)**2) if ddec_uncert_plot3d[k][i] != np.inf else racs_final[i]['e_dec']*u.arcsec
                
                racs_new = racs_final[i].copy()
                racs_new['ra_new'] = ra_new
                racs_new['dec_new'] = dec_new
                racs_new['e_ra_new'] = e_ra_new
                racs_new['e_dec_new'] = e_dec_new
                    
                racs_corrected.append(racs_new)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None


def CorrectRACS2(df_racs, full_offset_modelled_ra, full_offset_modelled_dec, dra_uncert_plot3d, ddec_uncert_plot3d, cal_sbids, radius, 
                directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected_fin = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME'][-8:]
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            load_racs_query = os.path.join(directory_racs, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            
            racs_final = [Table(np.load(os.path.join(load_racs_query, f'Beam_{i}.npy'))) for i in range(len(os.listdir(load_racs_query))) 
                    if os.path.isfile(os.path.join(load_racs_query, f'Beam_{i}.npy'))]
            
            for i in range(len(racs_final)):
                ra_fin = racs_final[i]['ra_new']*u.deg + full_offset_modelled_ra[k][i]*u.arcsec
                dec_fin = racs_final[i]['dec_new']*u.deg + full_offset_modelled_dec[k][i]*u.arcsec
                e_ra_fin = np.sqrt((racs_final[i]['e_ra_new']*u.arcsec)**2 + (dra_uncert_plot3d[k][i]*u.arcsec)**2) if dra_uncert_plot3d[k][i] != np.inf else racs_final[i]['e_ra_new']*u.arcsec
                e_dec_fin = np.sqrt((racs_final[i]['e_dec_new']*u.arcsec)**2 + (ddec_uncert_plot3d[k][i]*u.arcsec)**2) if ddec_uncert_plot3d[k][i] != np.inf else racs_final[i]['e_dec_new']*u.arcsec
                
                racs_fin = racs_final[i].copy()
                racs_fin['ra_fin'] = ra_fin
                racs_fin['dec_fin'] = dec_fin
                racs_fin['e_ra_fin'] = e_ra_fin
                racs_fin['e_dec_fin'] = e_dec_fin
                    
                racs_corrected_fin.append(racs_fin)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected_fin)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected_fin[i])
            
            racs_corrected_fin = []
    
    return None


def CorrectRACSLow3(df_racs, full_offset_modelled_ra, full_offset_modelled_dec, dra_uncert_plot3d, ddec_uncert_plot3d, cal_sbids, radius, 
                    directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME']
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            sbid_val = df_racs[k].iloc[0]['SBID']
            load_racs_query = os.path.join(directory_racs, f'selavy-image.i.{field_name}.SB{sbid_val}.cont.{field_name}.beam??.taylor.0.restored.components.xml')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSLow_Queries_{field_name[-7:]}_rad{radius}')
            
            racs = [Table.read(load_racs_query.replace('??', f'{i:02d}')) for i in range(len(glob.glob(load_racs_query)))]
            racs_final = CleanRACSLow3(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_new = racs_final[i]['col_ra_deg_cont'] + full_offset_modelled_ra[k][i]*u.arcsec
                dec_new = racs_final[i]['col_dec_deg_cont'] + full_offset_modelled_dec[k][i]*u.arcsec
                e_ra_new = np.sqrt((racs_final[i]['col_ra_err'].to(u.arcsec))**2 + (dra_uncert_plot3d[k][i]*u.arcsec)**2) if dra_uncert_plot3d[k][i] != np.inf else racs_final[i]['col_ra_err']
                e_dec_new = np.sqrt((racs_final[i]['col_dec_err'].to(u.arcsec))**2 + (ddec_uncert_plot3d[k][i]*u.arcsec)**2) if ddec_uncert_plot3d[k][i] != np.inf else racs_final[i]['col_dec_err']
                
                racs_new = racs_final[i].copy()
                racs_new['col_ra_deg_new'] = ra_new
                racs_new['col_dec_deg_new'] = dec_new
                racs_new['col_ra_err_new'] = e_ra_new
                racs_new['col_dec_err_new'] = e_dec_new
                    
                racs_corrected.append(racs_new)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None


def CorrectRACSLow3_2(df_racs, full_offset_modelled_ra, full_offset_modelled_dec, dra_uncert_plot3d, ddec_uncert_plot3d, cal_sbids, radius, 
                    directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected_fin = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME'][-7:]
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            load_racs_query = os.path.join(directory_racs, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSLow_Queries_{field_name}_rad{radius}')
            
            racs_final = [Table(np.load(os.path.join(load_racs_query, f'Beam_{i}.npy'))) for i in range(len(os.listdir(load_racs_query))) 
                        if os.path.isfile(os.path.join(load_racs_query, f'Beam_{i}.npy'))]
            
            for i in range(len(racs_final)):
                ra_fin = racs_final[i]['col_ra_deg_new']*u.deg + full_offset_modelled_ra[k][i]*u.arcsec
                dec_fin = racs_final[i]['col_dec_deg_new']*u.deg + full_offset_modelled_dec[k][i]*u.arcsec
                e_ra_fin = np.sqrt((racs_final[i]['col_ra_err_new']*u.arcsec)**2 + (dra_uncert_plot3d[k][i]*u.arcsec)**2) if dra_uncert_plot3d[k][i] != np.inf else racs_final[i]['col_ra_err_new']*u.arcsec
                e_dec_fin = np.sqrt((racs_final[i]['col_dec_err_new']*u.arcsec)**2 + (ddec_uncert_plot3d[k][i]*u.arcsec)**2) if ddec_uncert_plot3d[k][i] != np.inf else racs_final[i]['col_dec_err_new']*u.arcsec
                
                racs_fin = racs_final[i].copy()
                racs_fin['col_ra_deg_fin'] = ra_fin
                racs_fin['col_dec_deg_fin'] = dec_fin
                racs_fin['col_ra_err_fin'] = e_ra_fin
                racs_fin['col_dec_err_fin'] = e_dec_fin
                    
                racs_corrected_fin.append(racs_fin)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected_fin)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected_fin[i])
            
            racs_corrected_fin = []
    
    return None


def CorrectRACSMid1(df_racs, cal_sbids, radius, directory_racs1, directory_racs2, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME']
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            sbid_val = df_racs[k].iloc[0]['SBID']
            load_racs_query = sorted(glob.glob(os.path.join(directory_racs1, f'image.i.{field_name}.SB{sbid_val}.cont.{field_name}.beam??.taylor.0.restored.conv_comp.vot')) + 
                        glob.glob(os.path.join(directory_racs2, f'image.i.{field_name}.SB{sbid_val}.cont.{field_name}.beam??.taylor.0.restored.conv_comp.vot')))
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSMid_Queries_SB{sbid_val}_{field_name[-7:]}_rad{radius}')
            
            racs = [Table.read(load_racs_query[i]) if os.path.isfile(load_racs_query[i]) 
                    else Table(names=['ra', 'err_ra', 'dec', 'err_dec', 'peak_flux', 'err_peak_flux', 'int_flux', 'err_int_flux', 'a', 'err_a', 'b', 'err_b', 'pa', 'err_pa']) 
                    for i in range(len(load_racs_query))]
            racs_final = CleanRACSMid(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_new = racs_final[i]['ra']*u.deg + df_racs[k].iloc[i]['S0B_RA_Offset_Modelled']*u.arcsec
                dec_new = racs_final[i]['dec']*u.deg + df_racs[k].iloc[i]['S0B_DEC_Offset_Modelled']*u.arcsec
                e_ra_new = np.sqrt((racs_final[i]['err_ra']*u.arcsec)**2 + (df_racs[k].iloc[i]['Delta RA Uncertainity']*u.arcsec)**2) if df_racs[k].iloc[i]['Delta RA Uncertainity'] != np.inf else racs_final[i]['err_ra']
                e_dec_new = np.sqrt((racs_final[i]['err_dec']*u.arcsec)**2 + (df_racs[k].iloc[i]['Delta DEC Uncertainity']*u.arcsec)**2) if df_racs[k].iloc[i]['Delta DEC Uncertainity'] != np.inf else racs_final[i]['err_dec']
                
                racs_new = racs_final[i].copy()
                racs_new['ra_new'] = ra_new
                racs_new['dec_new'] = dec_new
                racs_new['err_ra_new'] = e_ra_new
                racs_new['err_dec_new'] = e_dec_new
                
                racs_corrected.append(racs_new)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None


def CorrectRACSMid1_2(df_racs, cal_sbids, radius, directory_racs_corr, directory_racs_fin):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME']
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            sbid_val = df_racs[k].iloc[0]['SBID']
            load_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSMid_Queries_SB{sbid_val}_{field_name[-7:]}_rad{radius}')
            save_racs_query = os.path.join(directory_racs_fin, f'{utc_time}_RACSMid_Queries_SB{sbid_val}_{field_name[-7:]}_rad{radius}')
            
            racs_final = [Table(np.load(os.path.join(load_racs_query, f'Beam_{i}.npy'))) for i in range(len(os.listdir(load_racs_query))) 
                    if os.path.isfile(os.path.join(load_racs_query, f'Beam_{i}.npy'))]
            # racs_final = CleanRACSMid(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_fin = racs_final[i]['ra_new']*u.deg + df_racs[k].iloc[i]['S2_RA_Offset_Modelled']*u.arcsec
                dec_fin = racs_final[i]['dec_new']*u.deg + df_racs[k].iloc[i]['S2_DEC_Offset_Modelled']*u.arcsec
                e_ra_fin = np.sqrt((racs_final[i]['err_ra_new']*u.arcsec)**2 + (df_racs[k].iloc[i]['S2_Delta_RA_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S2_Delta_RA_Uncertainty'] != np.inf else racs_final[i]['err_ra_new']
                e_dec_fin = np.sqrt((racs_final[i]['err_dec_new']*u.arcsec)**2 + (df_racs[k].iloc[i]['S2_Delta_DEC_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S2_Delta_DEC_Uncertainty'] != np.inf else racs_final[i]['err_dec_new']
                
                racs_fin = racs_final[i].copy()
                racs_fin['ra_fin'] = ra_fin
                racs_fin['dec_fin'] = dec_fin
                racs_fin['err_ra_fin'] = e_ra_fin
                racs_fin['err_dec_fin'] = e_dec_fin
                
                racs_corrected.append(racs_fin)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None


def CorrectRACSHigh1(df_racs, cal_sbids, radius, directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME']
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            sbid_val = df_racs[k].iloc[0]['SBID']
            load_racs_query = os.path.join(directory_racs, f'selavy-image.i.{field_name}.SB{sbid_val}.cont.{field_name}.beam??.taylor.0.restored.components.xml')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSHigh_Queries_{field_name[-7:]}_rad{radius}')
            
            racs_final = [Table.read(load_racs_query.replace('??', f'{i:02d}')) for i in range(len(glob.glob(load_racs_query)))]
            # racs_final = CleanRACSLow3(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_new = racs_final[i]['col_ra_deg_cont'] + df_racs[k].iloc[i]['S0B_RA_Offset_Modelled']*u.arcsec
                dec_new = racs_final[i]['col_dec_deg_cont'] + df_racs[k].iloc[i]['S0B_DEC_Offset_Modelled']*u.arcsec
                e_ra_new = np.sqrt((racs_final[i]['col_ra_err'].to(u.arcsec))**2 + (df_racs[k].iloc[i]['S0B_Delta_RA_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S0B_Delta_RA_Uncertainty'] != np.inf else racs_final[i]['col_ra_err']
                e_dec_new = np.sqrt((racs_final[i]['col_dec_err'].to(u.arcsec))**2 + (df_racs[k].iloc[i]['S0B_Delta_DEC_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S0B_Delta_DEC_Uncertainty'] != np.inf else racs_final[i]['col_dec_err']
                
                racs_new = racs_final[i].copy()
                racs_new['col_ra_deg_new'] = ra_new
                racs_new['col_dec_deg_new'] = dec_new
                racs_new['col_ra_err_new'] = e_ra_new
                racs_new['col_dec_err_new'] = e_dec_new
                    
                racs_corrected.append(racs_new)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None


def CorrectRACSHigh1_2(df_racs, cal_sbids, radius, directory_racs, directory_racs_corr):
    
    # Generating the corrected RACS catalogue by subtracting the modelled offsets from the original RACS catalogue positions
    racs_corrected = []

    for k in range(len(df_racs)): 
        cal_sbid_val = df_racs[k].iloc[0]['CAL_SBID']
        
        if str(cal_sbid_val) == cal_sbids:
            field_name = df_racs[k].iloc[0]['FIELD_NAME']
            utc_time = str(df_racs[k].iloc[0]['UTC_SCAN_START'])[0:10]
            sbid_val = df_racs[k].iloc[0]['SBID']
            load_racs_query = os.path.join(directory_racs, f'{utc_time}_RACSHigh_Queries_{field_name[-7:]}_rad{radius}')
            save_racs_query = os.path.join(directory_racs_corr, f'{utc_time}_RACSHigh_Queries_SB{sbid_val}_{field_name[-7:]}_rad{radius}')
            
            racs_final = [Table(np.load(os.path.join(load_racs_query, f'Beam_{i}.npy'))) for i in range(len(os.listdir(load_racs_query))) 
                    if os.path.isfile(os.path.join(load_racs_query, f'Beam_{i}.npy'))]
            # racs_final = CleanRACSLow3(df_racs[k], racs)
            
            for i in range(len(racs_final)):
                ra_fin = racs_final[i]['col_ra_deg_new']*u.deg + df_racs[k].iloc[i]['S2_RA_Offset_Modelled']*u.arcsec
                dec_fin = racs_final[i]['col_dec_deg_new']*u.deg + df_racs[k].iloc[i]['S2_DEC_Offset_Modelled']*u.arcsec
                e_ra_fin = np.sqrt((racs_final[i]['col_ra_err_new']*u.arcsec)**2 + (df_racs[k].iloc[i]['S2_Delta_RA_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S2_Delta_RA_Uncertainty'] != np.inf else racs_final[i]['col_ra_err_new']*u.arcsec
                e_dec_fin = np.sqrt((racs_final[i]['col_dec_err_new']*u.arcsec)**2 + (df_racs[k].iloc[i]['S2_Delta_DEC_Uncertainty']*u.arcsec)**2) if df_racs[k].iloc[i]['S2_Delta_DEC_Uncertainty'] != np.inf else racs_final[i]['col_dec_err_new']*u.arcsec
                
                racs_fin = racs_final[i].copy()
                racs_fin['col_ra_deg_fin'] = ra_fin
                racs_fin['col_dec_deg_fin'] = dec_fin
                racs_fin['col_ra_err_fin'] = e_ra_fin
                racs_fin['col_dec_err_fin'] = e_dec_fin
                    
                racs_corrected.append(racs_fin)
        
            # Saving the RACSLow Corrected Queries as NPY files
            os.makedirs(save_racs_query, exist_ok=True)
            for i in range(len(racs_corrected)):
                np.save(os.path.join(save_racs_query, f'Beam_{i}.npy'), racs_corrected[i])
            
            racs_corrected = []
    
    return None