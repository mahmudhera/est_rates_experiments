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
    
    ksizes = [21, 31, 41]
    colors = ['magenta', 'blue', 'green']
    
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
            
    # get the true_rates and estimated_rates
    true_rates = df[variable_column_name]
    estimated_rates = df[est_column_name]
    ksizes = df['ksize']
    
    # create a df for plotting
    plot_df = pd.DataFrame({
        'true'  : true_rates,
        'est'   : estimated_rates,
        'ksize' : ksizes
    })
    
    sns.boxplot(x='true', y='est', hue='ksize', data=plot_df, ax=ax, linewidth=0.5, showfliers=False)
    
    # set y ticks at 0.01, 0.02, 0.03, 0.04, 0.05
    ax.set_yticks([0.01, 0.02, 0.03, 0.04, 0.05])
    
    # set the x and y labels
    ax.set_xlabel(x_labels[variable_column_name])
    ax.set_ylabel(y_labels[variable_column_name])
    
    # show grid
    ax.grid(True)
    

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
    fig, axs = plt.subplots(2, 3, figsize=(8, 6))
    
    # plot the estimated rates for different fixed parameters
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.01, fixed_d=0.01, method_name=args.method, ax=axs[0][0])
    axs[0][0].set_title("$p_d = 0.01, d = 0.01$")
    
    plot_estimated_rates(df, fixed_ps=0.01, fixed_pd=None, fixed_d=0.01, method_name=args.method, ax=axs[0][1])
    axs[0][1].set_title("$p_s = 0.01, d = 0.01$")
    
    plot_estimated_rates(df, fixed_ps=0.01, fixed_pd=0.01, fixed_d=None, method_name=args.method, ax=axs[0][2])
    axs[0][2].set_title("$p_s = 0.01, p_d = 0.01$")
    
    plot_estimated_rates(df, fixed_ps=None, fixed_pd=0.05, fixed_d=0.05, method_name=args.method, ax=axs[1][0])
    axs[1][0].set_title("$p_d = 0.05, d = 0.05$")
    
    plot_estimated_rates(df, fixed_ps=0.05, fixed_pd=None, fixed_d=0.05, method_name=args.method, ax=axs[1][1])
    axs[1][1].set_title("$p_s = 0.05, d = 0.05$")
    
    plot_estimated_rates(df, fixed_ps=0.05, fixed_pd=0.05, fixed_d=None, method_name=args.method, ax=axs[1][2])
    axs[1][2].set_title("$p_s = 0.05, p_d = 0.05$")
    
    
    
    # tight layout to avoid overlap
    plt.tight_layout()
    
    # save the figure
    plt.savefig(args.output_filename)
    
    
if __name__ == "__main__":
    main()
    
    
    