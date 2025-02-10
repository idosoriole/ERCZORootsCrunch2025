import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from colour import Color
from matplotlib.ticker import FormatStrFormatter
import matplotlib.dates as mdates
from numpy.polynomial.polynomial import Polynomial
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import os
import time
import re
import random
import subprocess as S
import shutil 
from decimal import Decimal

def split_with_exception(case_name, exception_str):
    # Use a regular expression to split by '-' but ignore 'E-'
    parts = re.split(exception_str, case_name)
    return parts

def create_color_list(unique_values, palette):
# Create a color list
    cmap = plt.get_cmap(palette, len(unique_values))    # PiYG
    colors = []
    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        colors.append(mpl.colors.rgb2hex(rgba))
    return colors

def remove_outliers(df, col):
    Q1 = df[col].quantile(0.25)
    Q3 = df[col].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[col] >= lower_bound) & (df[col] <= upper_bound)]

def fix_Crunch_exp_over100(s):
    #Funtion to deal with the exportation problem crunchflow has when dealing with exponents greather than E+99
    if s == 'NaN':
        pass
    else:
        if 'E' not in s:        
            cont_m = 0
            cont_p = 0
            for i in range(len(s)):
                if s[i] == '-': 
                    cont_m+=1
                elif s[i] == '+': 
                    cont_p+=1
            if cont_m == 2:
                exp_pos = s.rfind('-')-1
            elif cont_p == 2:
                exp_pos = s.rfind('+')-1
                print('double +')
            elif (cont_m == 1) & (cont_p == 1) & (s[0] == '-'):    
                exp_pos = s.rfind('+')-1
            else:
                exp_pos = s.rfind('-')-1
            s = s[0:exp_pos+1] +'E'+s[exp_pos+1:len(s)]
        else:
            s = s
    return(s)

cols = ['Depth_m', 'year', 'Na_mM', 'K_mM', 'Ca_mM', 'Mg_mM', 'Si_mM']
def build_VMS_df(chem_df, start_year):
    df_field = chem_df.loc[:,cols].dropna()
    df_field = df_field.loc[df_field.year >= start_year]
    df_field = df_field.drop(columns = 'year')
    df_field.columns = ['Distance', 'Na+', 'K+', 'Ca++', 'Mg++', 'SiO2(aq)']
    df_field = df_field.groupby('Distance').agg(['mean', 'std']).reset_index()
    
    return df_field

def compile_sim_output(working_folder, case_name):
    BK_names = []
    cont = 0

    for filename in os.listdir(os.path.join(working_folder, case_name)):
        if 'timeseries' in filename: 
            data_f = list()
            with open(os.path.join(working_folder, case_name, filename)) as fp: 
                for cnt, line in enumerate(fp):
                    if cnt == 1:
                        cols = line.strip().split()[0:17] #17
                        cols[0] = 'time'
                    elif cnt > 1:
                        data = line.strip().split()[0:17] #17
                        data = list(map(fix_Crunch_exp_over100, data))
                        data = list(map(float, data))
                        data_f.append(data)
            BK_df = pd.DataFrame(data_f, columns = cols)
    return BK_df

def compile_Crunch_output(output_name, host_folder, case):
    #function to compile all the output files of a sensibility in a dictionnary containing the different parameters
    cont = 0

    for filename in os.listdir(os.path.join(host_folder, case)):
        if output_name in filename:       
            with open(os.path.join(host_folder, case, filename)) as fp:
                data_f = pd.DataFrame()
                for cnt, line in enumerate(fp):
                    if cnt == 0: #Recover time value
                        time = float("{}".format(line.strip()).split()[2])  
                    if cnt == 2: #recover column name
                        cols = "{}".format(line.strip()).split() 
                        cols.insert(0, 'time')
                    if cnt >= 3: #recover field value
                        aux = list()
                        data = "{}".format(line.strip())
                        aux.append(time)
                        for i in data.split():
                            i = fix_Crunch_exp_over100(i)
                            aux.append(float(i))

                        aux = pd.DataFrame([aux], columns = cols)
                        if cnt == 3: 
                            data_f = aux.copy()
                        else:
                            data_f = pd.concat([data_f, aux])

    data_f = data_f.sort_values(by=['time', 'Distance']).reset_index(drop = True)
    data_f.insert(0, 'Case', case)

    return data_f

def build_output_df(working_folder, case_name):
    BK_names = []
    cont = 0

    for filename in os.listdir(os.path.join(working_folder, case_name)):
        if 'timeseries' in filename: 
            data_f = list()
            with open(os.path.join(working_folder, case_name, filename)) as fp: 
                for cnt, line in enumerate(fp):
                    if cnt == 1:
                        cols = line.strip().split()[0:17] #17
                        cols[0] = 'time'
                    elif cnt > 1:
                        data = line.strip().split()[0:17] #17
                        data = list(map(fix_Crunch_exp_over100, data))
                        data = list(map(float, data))
                        data_f.append(data)
            BK_df = pd.DataFrame(data_f, columns = cols)

    dfs_output = []
    aux_output = []
    outvars = ['conc', 'totcon', 'saturation', 'exchange', 'gases'] 
    for outvar in outvars:
        aux = compile_Crunch_output(outvar, working_folder, case_name)
        aux_output.append(aux)

    aux_output = dict(zip(outvars, aux_output)) # create dict of outputs for the case ---> outputvariable : df_output
             
    return BK_df, aux_output

def plot_Figure3(dfs_output, cases_list, df_vms):

    spes_with_field_obs = ['Ca++', 'Mg++', 'Na+', 'K+', 'SiO2(aq)']

    # Create subplots
    fig, axs = plt.subplots(3, 3, figsize=(9, 9), sharex=False, sharey=True)
    legend_ax1 = fig.add_axes([0.62, 0.2, 0.35, 0.75], frame_on=False)  # Legend for cont == 6
    legend_ax2 = fig.add_axes([0.62, 0.1, 0.35, 0.25], frame_on=False)  # Legend for cont == 7
    legend_ax1.axis('off')  # Hide the external axis itself
    legend_ax2.axis('off')  # Hide the external axis itself

    colors = create_color_list(cases_list, 'tab20')[::-1]  # Reverse colors to match the cases

    for case_name, color_iter in zip(cases_list, colors):
        linestyle = ':' if 'No_Deep' in case_name else '-'
        legend_label = case_name

        # Loop over the species
        iter_cont = 0
        for cont, spe in enumerate(['Ca++', 'Mg++', 'Na+', 'K+', 'SiO2(aq)', 'Al+++', 'CO2(g)', 'O2(g)', 'pH']):
            row, col = divmod(cont, 3)  # Calculate row and column index for subplots
            ax = axs[row, col]

            if spe in ['CO2(g)', 'O2(g)']:  # Simulated Gases 
                df_t = dfs_output[case_name]['gases'].loc[dfs_output[case_name]['gases'].time == 50]
                ax.plot(df_t[spe] * 100 / 1.013, df_t['Distance'], label=legend_label, color=color_iter, linestyle=linestyle)
                ax.set_title(f'{spe} [%]')
                if spe == 'CO2(g)':
                    ax.text(x = 2, y = 9 , s = 'Field data \navailable in SI from \nTune et al., 2020', fontsize = 7)
                else:
                    ax.text(x = 10, y = 9 , s = 'Field data \navailable in SI from \nTune et al., 2020', fontsize = 7)

            elif spe == 'pH':  # Plot pH
                df_t = dfs_output[case_name]['conc'].loc[dfs_output[case_name]['conc'].time == 50]
                ax.plot(df_t[spe], df_t['Distance'], label=legend_label, color=color_iter, linestyle=linestyle)
                ax.set_title(spe)

            else:  # Other Species
                if iter_cont == 0 and spe in spes_with_field_obs:
                    ax.plot(df_vms[spe]['mean'], df_vms['Distance'], marker='o', linestyle='', color='grey')
                    ax.errorbar(df_vms[spe]['mean'], df_vms['Distance'], xerr=df_vms[spe]['std'], fmt='none', ecolor='black', capsize=5, alpha=0.1)
            
                df_t = dfs_output[case_name]['totcon'].loc[dfs_output[case_name]['totcon'].time == 50]
                ax.plot(df_t[spe] * 1000, df_t['Distance'], label=legend_label, color=color_iter, linestyle=linestyle)
                ax.set_title(f'{spe} [mM]')

            # Set common plot properties
            ax.set_ylim(18.9, -0.5)
            if spe == 'Na+':
                ax.set_xlim(1e-2, 1e0)
            elif spe == 'SiO2(aq)':
                ax.set_xlim(1e-2, 1e0)
            
            ax.set_xscale('log' if spe in ['Ca++', 'Mg++', 'Na+', 'K+', 'SiO2(aq)', 'Al+++'] else 'linear')
            ax.set_ylabel('Depth [m]' if col == 0 else '')

            # Handle legends for specific cases
            if cont == 6:
                handles, labels = ax.get_legend_handles_labels()
                legend_ax1.legend(handles, labels, loc=2, title=f'Simulation case:')
        iter_cont += 1
    
    plt.subplots_adjust(left=0.06, right=0.6, top=0.9, bottom=0.1, wspace=0.2, hspace=0.4)
    
    return fig

def compute_total_conc(df, index_above_watertable):
    total_conc = ((df['Ca++'].iloc[index_above_watertable] - df['Ca++'].iloc[0])*40 + \
    (df['Mg++'].iloc[index_above_watertable] - df['Mg++'].iloc[0])*24.3 +\
    (df['Na+'].iloc[index_above_watertable] - df['Na+'].iloc[0])*22.98 +\
    (df['K+'].iloc[index_above_watertable] - df['K+'].iloc[0])*39.09 +\
    (df['SiO2(aq)'].iloc[index_above_watertable] - df['SiO2(aq)'].iloc[0])*(28)) *1000
    return total_conc

def dfs_output_fun(working_folder, vars_to_assemble, case_name):
    print('running function dfs_output_fun...')
    cases = list()
    error_cases = list()
    cases = os.listdir(working_folder)

    dfs_output = []
    for case in cases:
        if os.path.isdir(os.path.join(working_folder, case)):
            # print('working on dfs_output %s ******************'%case)
            ##### Create dictionnary of dataframes for every selected output:
            aux_output = []
            outvars = vars_to_assemble # 'CaIsotopes','exchange', 'totexchange', 'DelGBiomass','gas', 'speciation', 'AqRate', 'pH'
            for outvar in outvars:
                # print(outvar)
                try:
                    # print(f'Launching compilation for {outvar}')
                    aux = compile_Crunch_output(outvar, working_folder, case)
                    aux_output.append(aux)
                except (ValueError):
                    if outvar=='pH':
                        pass
                    else:
                        print('*****ValueError in %s for %s'%(case, outvar))
                        pass
                except UnboundLocalError:
                    print('******UnboundLocalError in %s for %s'%(case, outvar))
                    pass
                except AssertionError:
                    print('******AssertionError in %s for %s'%(case, outvar))
                    pass
                except AttributeError: 
                    print('******AttributeError in %s for %s'%(case, outvar))
                    error_cases.append(case)
                    pass

        aux_output = dict(zip(outvars, aux_output)) # create dict of outputs for the case ---> outputvariable : df_output
        dfs_output.append(aux_output) #append to list containing output df's for all the cases

    dfs_output = dict(zip(cases, dfs_output)) #create dictionnary of case : df_outputs

    output_name = 'BK_'
    BK_dict = []
    BK_names = []
    cont = 0
    aux_BK_names = ['bottom', 'middle', 'top']
    # print('dfs_output_fun see case name is ', case_name)
    for filename in os.listdir(os.path.join(working_folder, case_name)):
        if 'BK_' in filename: 
            # BK_names.append(filename)
            print('working on %s'%filename, end = '\r')
            data_f = list()
            with open(os.path.join(working_folder, case_name, filename)) as fp: 
                for cnt, line in enumerate(fp):
                    if cnt == 1:
                        cols = line.strip().split()[0:17] #17
                        cols[0] = 'time'
                    elif cnt > 1:
                        data = line.strip().split()[0:17] #17
                        data = list(map(fix_Crunch_exp_over100, data))
                        data = list(map(float, data))
                        data_f.append(data)
            BK_df = pd.DataFrame(data_f, columns = cols)
            BK_dict.append(BK_df)
    BK_dict = dict(zip(aux_BK_names, BK_dict))
    return dfs_output, BK_dict

def compile_Crunch_output(output_name, host_folder, case):
    #function to compile all the output files of a sensibility in a dictionnary containing the different parameters

    files = list()
    SI_df = list()
    SI_files = list()
    cont = 0
    df_SI = list()

    for filename in os.listdir(os.path.join(host_folder, case)):
        if output_name in filename:   
            if ('totexchange' in filename) & (output_name == 'exchange'):
                next
            else: 
                output_num = int(filename[len(output_name):filename.find('.')]) #just to check if file is the last_one
                if output_name in ['pH', 'velocity']:
                    if (output_num == last_output_num) & (discard_output == 'yes'): #check condition for steady state reached
                        pass
                    else:
                        with open(os.path.join(host_folder, case, filename)) as fp: #with open(os.path.join(root, case, filename)) as fp:
                            data_f = list()
                            for cnt, line in enumerate(fp):
                                if cnt == 0:
                                    time = float("{}".format(line.strip()).split()[2])  
                                cols = ['time', output_name]
                                if cnt > 1:
                                    data = float("{}".format(line.strip()).split()[1])  
                                    data_f.append(time)
                                    data_f.append(data)
                elif output_name in ['speciation']:   
                    df = compile_single_speciation_file(host_folder, case, filename)
                else:        
                    with open(os.path.join(host_folder, case, filename)) as fp:
                        data_f = pd.DataFrame()
                        for cnt, line in enumerate(fp):
                            if cnt == 0: #Recover time value
                                time = float("{}".format(line.strip()).split()[2])  
                            if cnt == 2: #recover column name
                                cols = "{}".format(line.strip()).split() 
                                cols.insert(0, 'time')
                            if cnt >= 3: #recover field value
                                aux = list()
                                data = "{}".format(line.strip())
                                aux.append(time)
                                for i in data.split():
                                    i = fix_Crunch_exp_over100(i)
                                    aux.append(float(i))

                                aux = pd.DataFrame([aux], columns = cols)
                                if cnt == 3: 
                                    data_f = aux.copy()
                                else:
                                    data_f = pd.concat([data_f, aux])

            # To concatenate the different times
            if cont == 0:
                df_SI = data_f
            else: 
                df_SI = pd.concat([df_SI, data_f], ignore_index=True)
            cont += 1

    # print(f'\n***********\ncontents of df_SI for output_name={output_name} in case={case}: \n***********\n', len(df_SI))
    df_SI = df_SI.sort_values(by=['time', 'Distance']).reset_index(drop = True)
    df_SI.insert(0, 'Case', case)

    return df_SI

def transform_flow_to_flux(flow, watsat_level, bc_pco2_level, polynomial_coefficients):
    coeffs = polynomial_coefficients.get((watsat_level, bc_pco2_level))
    return sum(c * flow ** i for i, c in enumerate(coeffs)) if coeffs is not None else np.nan

def process_flux_data(fluxes_df, pco2_levels=[1.2, 2, 4], watsat_levels=[15, 23, 25, 30, 35, 40, 50], degree=2):
    """
    Code to transform ecosystem contribution estimates from Drainage space to F$_{CW}$ space by fitting a polynomial model between these two variables. 
    """    

    # Filter dataset based on selected BC_PCO2 and WatSat levels
    aux = fluxes_df.copy()
    aux = aux[aux.BC_PCO2.isin(pco2_levels) & aux.WatSat.isin(watsat_levels)]

    # Store polynomial coefficients
    polynomial_coefficients = {}

    # Fit polynomial regression for each (WatSat, BC_PCO2) combination
    for watsat_level in aux['WatSat'].unique():
        for bc_pco2_level in aux['BC_PCO2'].unique():
            df_subset = aux[(aux['WatSat'] == watsat_level) & (aux['BC_PCO2'] == bc_pco2_level)]
            if len(df_subset) > degree:  # Ensure enough data points
                p = Polynomial.fit(df_subset['Flow'], df_subset['Flux_CoxON'], degree)
                polynomial_coefficients[(watsat_level, bc_pco2_level)] = p.convert().coef

    # Apply transformation to dataset
    aux['Transformed_Flux_CoxON'] = aux.apply(
        lambda row: transform_flow_to_flux(row['Flow'], row['WatSat'], row['BC_PCO2'], polynomial_coefficients), axis=1
    )

    # Compute mean and 95% CI for each BC_PCO2 level
    mean_ci_curves = {}
    for bc_pco2_level in aux['BC_PCO2'].unique():
        df_subset = aux[aux['BC_PCO2'] == bc_pco2_level]
        grouped = df_subset.groupby('Flow')['Ecosystem_pct'].agg(['mean', 'sem']).reset_index()
        ci = stats.norm.interval(0.95, loc=grouped['mean'], scale=grouped['sem'])
        
        mean_ci_curves[bc_pco2_level] = {
            'Flow': grouped['Flow'],
            'Mean': grouped['mean'],
            'CI_Lower': ci[0],
            'CI_Upper': ci[1]
        }

    # Find best-matching WatSat for each BC_PCO2
    def curve_distance(y1, y2):
        return np.sqrt(np.sum((y1 - y2) ** 2))

    best_watsat_for_MUL = {}
    for bc_pco2_level, curves in mean_ci_curves.items():
        df_subset = aux[aux['BC_PCO2'] == bc_pco2_level]
        best_match = {'Mean': None, 'CI_Lower': None, 'CI_Upper': None}
        min_distances = {'Mean': float('inf'), 'CI_Lower': float('inf'), 'CI_Upper': float('inf')}
        
        for watsat_level in df_subset['WatSat'].unique():
            df_watsat = df_subset[df_subset['WatSat'] == watsat_level]
            grouped_watsat = df_watsat.groupby('Flow')['Ecosystem_pct'].mean().reindex(curves['Flow']).fillna(method='ffill').fillna(method='bfill')

            for key in ['Mean', 'CI_Lower', 'CI_Upper']:
                distance = curve_distance(curves[key], grouped_watsat)
                if distance < min_distances[key]:
                    min_distances[key] = distance
                    best_match[key] = watsat_level
        
        best_watsat_for_MUL[bc_pco2_level] = best_match

    # Transform M-U-L curves in Flux_CoxON space
    transformed_MUL = {}
    for bc_pco2_level, curves in mean_ci_curves.items():
        best_watsat = best_watsat_for_MUL[bc_pco2_level]
        
        transformed_curves = {
            'Flux_CoxON': np.array([transform_flow_to_flux(flow, best_watsat['Mean'], bc_pco2_level, polynomial_coefficients) for flow in curves['Flow']]),
            'Mean': np.array([transform_flow_to_flux(flow, best_watsat['Mean'], bc_pco2_level, polynomial_coefficients) for flow in curves['Flow']]),
            'CI_Lower': np.array([transform_flow_to_flux(flow, best_watsat['CI_Lower'], bc_pco2_level, polynomial_coefficients) for flow in curves['Flow']]),
            'CI_Upper': np.array([transform_flow_to_flux(flow, best_watsat['CI_Upper'], bc_pco2_level, polynomial_coefficients) for flow in curves['Flow']])
        }
        
        transformed_MUL[bc_pco2_level] = transformed_curves

    return mean_ci_curves, best_watsat_for_MUL, transformed_MUL, aux, polynomial_coefficients
