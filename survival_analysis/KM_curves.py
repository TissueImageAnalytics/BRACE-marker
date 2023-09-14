import pickle
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import math
import numpy as np
import matplotlib.pyplot as plt
import os

## script for generating Kaplan Meier (KM) curves for BRACE marker and other clinicopathological variables
## note that the Mean C-Indices might differ slightly from the reported results because of the 1000 bootstrap runs
## this script also generate .csv file which can then be used with the accompanied R script (/Forest_plots/forest_plot_R_script.R) to generate forest plots

def apply_preprocess(df_temp, subset, feats_list, time_col, event_col, censor_at, discov_mean_std_path):
    df_temp = df_temp.fillna(0)

    ## apply censoring. 
    df_temp.loc[df_temp[time_col] > censor_at, event_col] = 0
    df_temp.loc[df_temp[time_col] > censor_at, time_col] = censor_at

    ## convert all events other than 1 to 0
    df_temp.loc[df_temp[event_col] > 1, event_col] = 0

    ## remove chemotherapy cases
    df_temp.drop(df_temp[df_temp['Chemotherapy'] == 1].index, inplace=True)

     ##apply lymph node filtering i.e consider lymph node negative, or consider lymph node 0-3
    if subset == 'LN+':
        df_temp.drop(df_temp[df_temp['Endocrine Therapy'] != 1].index, inplace=True)

    if subset == 'LN-':
        df_temp.drop(df_temp[df_temp['Endocrine Therapy'] != 1].index, inplace=True)
        df_temp.drop(df_temp[df_temp['Lymph Node status'] == 1].index, inplace=True)

    if discov_mean_std_path != -1:
        train_mean_std_df = pd.read_csv(discov_mean_std_path)

        val_data_temp = df_temp[feats_list]
        val_data_temp = (val_data_temp - train_mean_std_df['mean'].values.tolist())/ train_mean_std_df['std'].values.tolist() ##standardise using training mean/std
        val_data_temp['time'] = df_temp[time_col]
        val_data_temp['event'] = df_temp[event_col]
    else:
        val_data_temp = df_temp[feats_list]
        val_data_temp['time'] = df_temp[time_col]
        val_data_temp['event'] = df_temp[event_col]

    return val_data_temp

##add new features based on some equation e.g Magee etc
def add_new_features(df):
    ## New Magee equation2: recurrence score = 18.8042 + (NPI)*2.34123 + (ER H-score)*(-0.03749) + (PR H-score)*(-0.03065) + (0 for negative, 1.82921 for equivocal and 11.51378 for HER2 positive)+tumour_size*0.04267.
    
    PR_Hscore = 0 ## PR Hscore not available in the clinical file
    HER2 = 0 ## all cases are HER2 negative

    df_new = pd.DataFrame()
    for ind, row in df.iterrows():
        
        magee_new2 = 18.8042 + row['NPI']*2.34123 + row['ER Hscore']*(-0.03749) + (PR_Hscore)*(-0.03065) + HER2 + row['Tumour_Size']*0.04267
        row['Magee_new2'] = magee_new2

        magee_new2_no_npi = 18.8042 + row['ER Hscore']*(-0.03749) + (PR_Hscore)*(-0.03065) + HER2 + row['Tumour_Size']*0.04267
        row['Magee_new2_no_npi'] = magee_new2_no_npi

        df_new = df_new.append(row, ignore_index=True)
    
    return df_new

def generate_KM_curves(path_to_model_a, path_to_model_b, features_df, cutoff_value, path_to_save_results, event_name, model_name, LN_status):
    ## load the Cox PH model that was fitted on the discovery set of Cohort-A
    with open(path_to_model_a, 'rb') as file:
        fitted_model = pickle.load(file)
    
    ## make predictions on the test set
    score_val = fitted_model.predict_partial_hazard(features_df)
    
    # Use the cut off value selected based on the discovery set to stratify subjects in the validation set into high and low risk groups
    upper = score_val >= cutoff_value
    T_upper_test = features_df['time'][upper]
    E_upper_test = features_df['event'][upper]
    lower = score_val < cutoff_value
    T_lower_test = features_df['time'][lower]
    E_lower_test = features_df['event'][lower]

    # Log-rank test: check if there is any significant difference between the groups being compared
    results = logrank_test(T_lower_test, T_upper_test, E_lower_test, E_upper_test)
    val_pvalue = round(results.p_value, 4)

    features_df[model_name] = score_val
    model_score_df = features_df[[model_name, 'time', 'event']]

    with open(path_to_model_b, 'rb') as file:
        fitted_model = pickle.load(file)

    val_cindex = fitted_model.score(features_df, scoring_method="concordance_index")
    
    ## compute the hazard ratios and confidence intervals for the test set
    cph = CoxPHFitter(baseline_estimation_method='breslow', l1_ratio=0.5, penalizer=0.001).fit(model_score_df, 'time', 'event')
    val_hzratio = round(cph.hazard_ratios_[model_name], 2)

    try:
        vhz_ci_low = round(math.exp(cph.confidence_intervals_.values[0][0]), 2)
    except OverflowError:
        vhz_ci_low = float('inf')
    
    try:
        vhz_ci_high = round(math.exp(cph.confidence_intervals_.values[0][1]), 2)
    except OverflowError:
        vhz_ci_high = float('inf')

    ## 1000 bootstraps on the test set to compute mean +- standard deviation for the C-Index
    c_indices_val = []

    rng = np.random.RandomState()

    for k in range(1000):
            EE = np.array(model_score_df["event"]).astype(int)

            index_val = list(rng.choice( np.nonzero(EE == 0)[0], size=len(EE) - np.sum(EE), replace=True)) + \
                        list(rng.choice( np.nonzero(EE == 1)[0], size=np.sum(EE), replace=True))
            bootval_set = model_score_df.iloc[index_val]
            bootvcindex = concordance_index(bootval_set["time"], -bootval_set[model_name], bootval_set["event"])
            c_indices_val.append(bootvcindex)

    cindex_std = round(np.asarray(c_indices_val).std(), 2)
    cindex_mean = round(np.asarray(c_indices_val).mean(), 2)

    ## create and save the KM curves
    # Initializing the KaplanMeierModel for each group
    km_upper = KaplanMeierFitter()
    km_lower = KaplanMeierFitter()
    font_size = 18
    fig_size = 10
    fig = plt.figure(figsize=(fig_size, fig_size-2)) ##adjust according to font size
    ax = fig.add_subplot(111)

    ax.set_xlabel('', fontsize=font_size)
    ax.set_ylabel('', fontsize=font_size)
            
    if val_pvalue < 0.0001:
        val_pvalue_txt = 'p < 0.0001'
    else:
        val_pvalue_txt = 'p = ' + str(np.round(val_pvalue, 4))

    from matplotlib.offsetbox import AnchoredText
    ax.add_artist(AnchoredText(val_pvalue_txt, loc=1, frameon=False, prop=dict(size=font_size)))
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)

    if event_name == 'Distant Metastasis':
        event_type = 'DMFS'
    elif event_name == 'Survival Status':
        event_type = 'BCSS'

    ax = km_upper.fit(T_upper_test, event_observed=E_upper_test, label=model_name+'-high').plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 5}, color='r', ci_show=False, xlabel='Months', ylabel= event_type + ' Probability')
    ax = km_lower.fit(T_lower_test, event_observed=E_lower_test, label=model_name+'-low').plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 5}, color='b', ci_show=False, xlabel='Months', ylabel= event_type + ' Probability')
    from lifelines.plotting import add_at_risk_counts
    add_at_risk_counts(km_upper, km_lower, ax=ax, fig=fig, fontsize=int(font_size*1) )
    ax.get_legend().remove()
    plt.subplots_adjust(bottom=0.4)
    plt.subplots_adjust(left=0.2)
    plt.savefig(path_to_save_results + model_name + '_' + event_type + '_' + LN_status + '.pdf', format='pdf', dpi=600)

    return val_pvalue, cindex_mean, cindex_std, val_hzratio, vhz_ci_low, vhz_ci_high, score_val

if __name__ == '__main__':
    print('Generating survival results...')
    
    ## in the provided clinical files the time and event columns are as follows
    event_of_interest = {'DMFS': {'time': 'TTDM/ month', 'event': 'Distant Metastasis'}, 'BCSS': {'time': 'Breast cancer specific survival/ month', 'event': 'Survival Status'}}
    censor_months = 120 #in months. e.g 10 years = 120 months. 180 for 15 years, 240 for 20 years. Use -1 if no censoring is required 
    
    for cohort in ['Cohort-B']: ## ['Cohort-A', 'Cohort-B]
        print('Cohort: ', cohort)
        results_df = pd.DataFrame()  ## Dataframe to save the results

        ## path to a csv file containing the features along with the time_to_event and event
        ## the provided clinical file is for the external cohort (ie. Cohort-B in the manuscript)
        
        feature_df = pd.read_csv('./features/' + cohort + '.csv')

        ## path where the KM plots should be saved
        results_main_path =  './results/'
        results_save_path = results_main_path + cohort + '/'
        
        if not os.path.isdir(results_main_path):
            os.mkdir(results_main_path)
        
        if not os.path.isdir(results_save_path):
            os.mkdir(results_save_path)

        for event_type in ['DMFS', 'BCSS']: ## for each event type i.e. Distant Metastatsis Free Survival and Breast Cancer Specific Survival
            print('Event type: ', event_type)
            for subset_name in ['LN-', 'LN+']: ## LN- means no postive lymph nodes, LN+ means 0 to 3 postives lymph nodes
                for model_name in ['Grade', 'NPI', 'BRACE']: ## three types of prognostic markers. the first two (Grade and NPI) are based on clinicopathological variables whereas the last (BRACE) is based on ML features
                    ## path to CoxPH model fitted on initial features. The score from this model_a is then used as a new variable to fit model_b below.
                    path_to_saved_model_a = './models/fitted_' + model_name + '_cox_model_a.pkl'

                    ## path to CoxPH model fitted on the score obtained from model_a above
                    path_to_saved_model_b = './models/fitted_' + model_name + '_cox_model_b.pkl'

                    ## path to the discovery statistics used to standardise the data
                    if model_name == 'BRACE':
                        discovery_mean = './models/' + model_name + '_discov_mean_std.csv'
                    else:
                        discovery_mean = -1

                    if model_name == 'BRACE':
                        ## list of ML features included in BRACE marker
                        feats_list = ['per_g1','per_g2','per_g3','tumor_pixels_per','inter_WSI_S_density_LL_cooccur','pcontrast','stromal_contrast','imitos_70p_score_thresh1_8_thresh2_20_mj']
                        cutoff_value = 1.28
                    elif model_name == 'NPI':
                        feats_list = ['NPI']
                        cutoff_value = 0.94
                    elif model_name == 'Grade':
                        feats_list = ['Grade']
                        cutoff_value = 2 #2
                    
                    time_column = event_of_interest[event_type]['time']
                    event_column = event_of_interest[event_type]['event']
                    
                    df = apply_preprocess(feature_df, subset_name, feats_list, time_column, event_column, censor_months, discovery_mean)

                    pvalue, cindex, cindex_std, hzratio, hz_ci_low, hz_ci_high, score = generate_KM_curves(path_to_saved_model_a, path_to_saved_model_b, df, cutoff_value, results_save_path, event_column, model_name, subset_name)
                    results_df = pd.concat([results_df, pd.DataFrame([{'Event': event_type, 'LN status':subset_name, 'n': df.shape[0], 'Model': model_name, 
                                                                    'Mean_C-Index': cindex, 'C-Index_Std': cindex_std, 'p-value': pvalue, 
                                                                    'HR': hzratio, 'CI-low': hz_ci_low, 'CI-high': hz_ci_high}])], ignore_index=True)
                    
                    ## save the features for Forest plots generation in R script (/Forest_plots/forest_plot_R_script.R)
                    if model_name == 'BRACE' and event_type == 'DMFS' and subset_name == 'LN+':
                        
                        df_new = df.copy()
                        df_new['BRACE'] = score
                        df_new.drop(feats_list, axis=1, inplace=True)
                        df_new.drop(['time', 'event'], axis=1, inplace=True)                              
                        
                        ## add the required columns
                        ## rename some of the column causing issues in R script
                        for c, cn in zip(['Invasive Tumour Size (cm)', 'Associated LCIS', 'Associated DCIS'], ['Tumour_Size', 'Associated_LCIS', 'Associated_DCIS']):
                            df_new[cn] = feature_df[c]

                        for c in ['Grade', 'NPI', 'Lymph Node status', 'Multifocality', 'Age at Diagnosis', 'ER Hscore', 'PR Hscore', 'LVI',
                                  'Survival Status', 'Breast cancer specific survival/ month', 'Distant Metastasis', 'TTDM/ month']:
                            df_new[c] = feature_df[c]
                        
                        ## change grade 1,2 to 1 and grade 3 to 2
                        df_new.loc[df_new['Grade'] < 3, 'Grade'] = 1
                        df_new.loc[df_new['Grade'] > 2, 'Grade'] = 2
                        
                        if cohort == 'Cohort-A':
                            df_new = add_new_features(df_new)

                        ## this csv file is used as an input to the R script under /Forest_plots/forest_plot_R_script.R to generate the forest plots
                        df_new.to_csv('./features/' +  cohort + '_for_forest_plots.csv')
        
        results_df.to_csv(results_save_path + 'results.csv')
        print(results_df)
        
    print('End')