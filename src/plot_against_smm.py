import argparse
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


# arguments: --estimated_rates_file, --output_filename, --method
def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot estimated rates from a file containing true rates and estimated rates.")
    
    # Add arguments
    parser.add_argument('--estimated_rates_file', type=str, required=True,
                        help='Path to the file containing estimated rates.')
    parser.add_argument('--output_filename', type=str, required=True, help='Output filename for the plot.')
    parser.add_argument('--method', type=str, required=True, 
                        help='Method used for estimation (can be either "linear" or "poly").')
    
    # Parse the arguments
    args = parser.parse_args()
    return args


def parse_data(estimated_rates_filename, method_name):
    """
    Read the estimated rates from the file and return them as a dataframe
    """
    # the column names should be:
    # ps	pd	d	seed	ksize	subst_rate_lin	del_rate_lin	ins_rate_lin	subst_rate_poly	del_rate_poly	ins_rate_poly	subst_rate_smm
    
    # parse using pandas    
    df = pd.read_csv(estimated_rates_filename, sep="\t")
    return df
    
    
def plot_estimated_rates(df, fixed_ps, fixed_pd, fixed_d, method_name, ax):
    """
    Plot estimated rates for two fixed parameters. Only two parameters can be fixed at a time.
    The third needs to be None. Given the method name, the proper estimated rates are selected.
    The parameter that is set as None is in the x-axis, and estimated rates of that parameter is
    in the y-axis.
    """
    
    # check if the method name is valid
    if method_name not in ["linear", "poly"]:
        raise ValueError("Method name should be either 'linear' or 'poly'")
    
    # check if the fixed parameters are valid. None means variable
    column_names = ["ps", "pd", "d"]
    est_column_names_poly = ["subst_rate_poly", "del_rate_poly", "ins_rate_poly"]
    est_column_names_lin = ["subst_rate_lin", "del_rate_lin", "ins_rate_lin"]
    x_labels = {
        'ps': 'True Substitution Rate',
        'pd': 'True Deletion Rate',
        'd' : 'True Insertion Rate'
    }
    y_labels = {
        'ps': 'Estimated Substitution Rate',
        'pd': 'Estimated Deletion Rate',
        'd' : 'Estimated Insertion Rate'
    }
    
    num_variables, variable_column_name = 0, None
    for i, param_value in enumerate([fixed_ps, fixed_pd, fixed_d]):
        if param_value is None:
            num_variables += 1
            variable_column_name = column_names[i]
            est_column_name = est_column_names_poly[i] if method_name == "poly" else est_column_names_lin[i]
            
    if num_variables != 1:
        raise ValueError("Only one parameter can be variable at a time. The others should be fixed.")
    
    # filter the dataframe based on the fixed parameters
    for column_name, value in zip(column_names, [fixed_ps, fixed_pd, fixed_d]):
        if value is not None:
            df = df[df[column_name] == value]
            
    # filter for ksize = 21
    df = df[df['ksize'] == 21]
            
    # get the true_rates and estimated_rates
    true_rates = df[variable_column_name]
    estimated_rates_fmm = df[est_column_name]
    estimated_rates_smm = df['subst_rate_smm']
    
    true_rates_to_fmm_estimates = {}
    for true_rate, estimate in zip(true_rates, estimated_rates_fmm):
        if true_rate not in true_rates_to_fmm_estimates:
            true_rates_to_fmm_estimates[true_rate] = []
        true_rates_to_fmm_estimates[true_rate].append(estimate)
        
    true_rates_to_smm_estimates = {}
    for true_rate, estimate in zip(true_rates, estimated_rates_smm):
        if true_rate not in true_rates_to_smm_estimates:
            true_rates_to_smm_estimates[true_rate] = []
        true_rates_to_smm_estimates[true_rate].append(estimate)
        
    # create boxplots
    #sns.boxplot(x=true_rates, y=estimated_rates_fmm, ax=ax, color='blue', linewidth=0.2, width=0.3)
    #sns.boxplot(x=true_rates, y=estimated_rates_smm, ax=ax, color='red', linewidth=0.2, width=0.3)
    
    # for every ps, compute the mean values
    mean_estimates_fmm = []
    mean_estimates_smm = []
    var_estimates_fmm = []
    var_estimates_smm = []
    for ps in true_rates_to_fmm_estimates.keys():
        mean_estimates_fmm.append(sum(true_rates_to_fmm_estimates[ps]) / len(true_rates_to_fmm_estimates[ps]))
        var_estimates_fmm.append(sum([(x - mean_estimates_fmm[-1]) ** 2 for x in true_rates_to_fmm_estimates[ps]]) / len(true_rates_to_fmm_estimates[ps]))
        mean_estimates_smm.append(sum(true_rates_to_smm_estimates[ps]) / len(true_rates_to_smm_estimates[ps]))
        var_estimates_smm.append(sum([(x - mean_estimates_smm[-1]) ** 2 for x in true_rates_to_smm_estimates[ps]]) / len(true_rates_to_smm_estimates[ps]))
        
    std_estimates_fmm = [var ** 0.5 for var in var_estimates_fmm]
    std_estimates_smm = [var ** 0.5 for var in var_estimates_smm]
        
    # plot the mean estimates
    ax.scatter(true_rates_to_fmm_estimates.keys(), mean_estimates_fmm, color='blue', label='Our estimates (k=21)', s=4)
    ax.scatter(true_rates_to_smm_estimates.keys(), mean_estimates_smm, color='red', label='SMM estimates (k=21)', s=4)
    
    #ax.plot(true_rates_to_fmm_estimates.keys(), mean_estimates_fmm, color='blue', linestyle='--', linewidth=0.5)
    #ax.plot(true_rates_to_smm_estimates.keys(), mean_estimates_smm, color='red', linestyle='--', linewidth=0.5)
    
    # manually show error bars
    for i, ps in enumerate(true_rates_to_fmm_estimates.keys()):
        ax.errorbar(ps, mean_estimates_fmm[i], yerr=std_estimates_fmm[i], color='blue', linewidth=0.2, capsize=5)
        ax.errorbar(ps, mean_estimates_smm[i], yerr=std_estimates_smm[i], color='red', linewidth=0.2, capsize=5)
    
    # plot a y = x line
    ax.plot([0.01, 0.05], [0.01, 0.05], color='black', linestyle='--', linewidth=0.8, label='True subst. rate')
    
    # set x ticks
    ax.set_xticks([0.01, 0.02, 0.03, 0.04, 0.05])
    
    # rotate x ticks 45
    ax.tick_params(axis='x', rotation=45)
    
    # legend
    ax.legend(loc='upper left', fontsize=6, title_fontsize=6)
    
    # set the x and y labels
    ax.set_xlabel(x_labels[variable_column_name])
    ax.set_ylabel(y_labels[variable_column_name])
    
    
    

def main():
    # set sizes of text
    plt.rcParams['axes.titlesize'] = 8
    plt.rcParams['axes.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    
    # make the legend smaller
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['legend.title_fontsize'] = 6
    
    # parse the arguments
    args = parse_arguments()
    
    # read the data
    df = parse_data(args.estimated_rates_file, args.method)
    
    # create a figure with three subplots
    fig, axs = plt.subplots(2, 2, figsize=(7, 6))
    
    # plot the estimated rates for different fixed parameters
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.01, fixed_d=0.01, method_name=args.method, ax=axs[0][0])
    axs[0][0].set_title("$p_d = 0.01, d = 0.01$")
    
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.01, fixed_d=0.05, method_name=args.method, ax=axs[0][1])
    axs[0][1].set_title("$p_d = 0.01, d = 0.05$")
    
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.05, fixed_d=0.01, method_name=args.method, ax=axs[1][0])
    axs[1][0].set_title("$p_d = 0.05, d = 0.01$")
    
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.05, fixed_d=0.05, method_name=args.method, ax=axs[1][1])
    axs[1][1].set_title("$p_d = 0.05, d = 0.05$")
    
    # tight layout to avoid overlap
    plt.tight_layout()
    
    # save the figure
    plt.savefig(args.output_filename)
    #plt.show()
    
    
if __name__ == "__main__":
    main()
    
    
    