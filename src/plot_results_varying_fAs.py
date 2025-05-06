import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# arguments: --estimated_rates_file, --output_filename

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot estimated rates from a file containing true rates and estimated rates.")
    
    # Add arguments
    parser.add_argument('--estimated_rates_file', type=str, required=True,
                        help='Path to the file containing estimated rates.')
    parser.add_argument('--output_filename', type=str, required=True, help='Output filename for the plot.')
    
    # float argument: true_rate
    parser.add_argument('--true_rate', type=float, default=0.01, help='True mutation rate to plot against.')
    
    # Parse the arguments
    args = parser.parse_args()
    return args


def parse_data(estimated_rates_filename):
    """
    Read the estimated rates from the file and return them as a dataframe
    """
    # the column names should be:
    # fA	simulation	ksize	subst_rate	del_rate	ins_rate
    
    # parse using pandas    
    df = pd.read_csv(estimated_rates_filename, sep="\t")
    
    # remove if fA < 0.2 or > 0.3
    df = df[(df['fA'] > 0.2) & (df['fA'] < 0.3)]
    
    return df


def plot(df, which_rate_to_plot, ax, true_rate):
    # plot true rate horizontal line
    #ax.axhline(y=true_rate, color='grey', linestyle='-.', zorder=1, linewidth=0.8)
    
    # plot boxplot
    sns.boxplot(x='fA', y=which_rate_to_plot, hue='ksize', data=df, ax=ax, linewidth=0.5, showfliers=False) 
    
    
    pretty_rate_names = {
        'subst_rate': '$p_s$',
        'del_rate': '$p_d$',
        'ins_rate': '$d$'
    }
    
    # x and y labels
    ax.set_xlabel('fraction of As in $S$')
    ax.set_ylabel('Estimated ' + pretty_rate_names[which_rate_to_plot])
    ax.legend(title='$k$-mer size')
    
    # set x ticks to 90 degrees, make sure two decimal points are shown
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    
    # add a y tick at true_rate
    if which_rate_to_plot == 'subst_rate':
        y_ticks = ax.get_yticks() + [true_rate]
        y_ticks = sorted(list(set(y_ticks)))
        ax.set_yticks(y_ticks)
    
    
    # set y limit to -1 to 1.0
    #if which_rate_to_plot == 'subst_rate':
        #ax.set_ylim(-0.1, 0.1)
    
    
    
    
def main():
    # Parse the arguments
    args = parse_arguments()
    
    # Read the data
    df = parse_data(args.estimated_rates_file)
    
    # set text to serif font
    #plt.rcParams['font.family'] = 'serif'
    
    plt.rcParams['axes.titlesize'] = 8
    plt.rcParams['axes.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    
    # make the legend smaller
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['legend.title_fontsize'] = 6
    
    # Create a figure with subplots
    fig, axs = plt.subplots(1, 3, figsize=(8, 2.5))
    
    # Plot estimated rates for each rate
    plot(df, 'subst_rate', axs[0], args.true_rate)
    plot(df, 'del_rate', axs[1], args.true_rate)
    plot(df, 'ins_rate', axs[2], args.true_rate)
    
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(args.output_filename)
    
    
if __name__ == "__main__":
    main()