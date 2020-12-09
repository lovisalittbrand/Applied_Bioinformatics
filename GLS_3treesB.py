'''Run visualization script like this:

    python3 GLS_5.py file_likelihoods.sitelh file_partition.nexus PLOT_STYLE

where PLOT_STYLE can be any or all of the following: 
-gs     - plots difference in log-likelihoods per gene, sorts them and shows them on one side of the x-axis (DEFAULT)
-ss     - plots difference in log-likelihoods per site
'''

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
import argparse
import numpy as np
import time

# Parser function that allows to have optional commands
def parseArguments():
    my_parser = argparse.ArgumentParser(description='Display log-likelihood differences obtained from 2 conflicting phylogenetic trees')
    # Mandatory arguments
    my_parser.add_argument("file_likelihood", help="One file with log-likelihoods from 2 trees", type=str)
    my_parser.add_argument("nexus_partition", help="One Nexus partition file with gene boundaries", type=str)
    my_parser.add_argument("-gs", "--style1_plot_gene_sort", help="Display sorted log-likelihood per gene", action='store_true')
    
    # Optional arguments
    my_parser.add_argument("-gt", "--style2_plot_gene_two_sides", help="Display separated log-likelihood per gene", action='store_true')
    my_parser.add_argument("-ss", "--style3_plot_site", help="Display log-likelihood per site", action='store_true')  
    
    args = my_parser.parse_args()
    return args

args = parseArguments()  
   
# Read the input agruments to variables
file_likelihood = args.file_likelihood
nexus_partition = args.nexus_partition
style1_sorted = args.style1_plot_gene_sort
style2_two_sides = args.style2_plot_gene_two_sides
style3_mode_site=args.style3_plot_site

######################## CODE ########################

# import log likelihood per site for 3 different topologies: T1, T2 and T3 
with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, index_col=0, lineterminator='\n', skiprows=1)

# calculate the intensity of the phylogenetic signal 
df_site_lk.loc['Diff_site']=np.absolute(df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values) + np.absolute(df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree3'].values) + np.absolute(df_site_lk.loc['Tree2'].values - df_site_lk.loc['Tree3'].values)
df_site_lk.loc['Support']= df_site_lk.loc[['Tree1','Tree2','Tree3']].idxmax(axis=0)
df_site_lk.head()

#export the dataframe with genomic position and difference in site likelihood to .csv file
out=df_site_lk.T
out = pd.DataFrame(out).reset_index()
out.columns = ['Genomic_position', 'Tree1', 'Tree2', 'Tree3','Diff_site', 'Support']
col_change = ['Tree1', 'Tree2', 'Tree3', 'Diff_site']
out[col_change] = out[col_change].apply(pd.to_numeric, downcast='float').fillna(0)
timestr2 = time.strftime("_%Y_%m_%d-%H_%M_%S")
out.to_csv(file_likelihood+timestr2+'_sitelk.csv', float_format='%.5f', header = True, index=False) #rounds nr to 2 decimals

# import nexus partition file and extract gene name, its beginning, end and length (new parser)
gene_data = dict()
# Check if partition file is available 
if (nexus_partition != None):
  # Open file and read each line in search for lines that starts with charstet 
  with open(nexus_partition, 'r') as file:
    for line in file:
      # Regex that fetch the name of teh gene and start and stop position and then puts the length and name in a dict
      m = re.search(r"charset (.+) = (\d+)-(\d+)", line)
      if m != None:
        gene_data[m[1]] = [int(m[2]), int(m[3]), int(m[3])-int(m[2])+1]
# If no partition file is specified, create a gene as long as the sequence 
else:
  gene_data = {'No genes specified': len(df_site_lk)}

df_genes = pd.DataFrame.from_dict(gene_data, orient='index').reset_index()
df_genes.columns = ['gene_name', 'gene_begin', 'gene_end', 'gene_length']
#########################

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS_Tree1 = []  #gene likelihoods
GLS_Tree2 = []
GLS_Tree3 = []

for k in range(0, len(df_genes)):
    GLS_Tree1.append(sum(df_site_lk.loc['Tree1'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
    GLS_Tree2.append(sum(df_site_lk.loc['Tree2'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
    GLS_Tree3.append(sum(df_site_lk.loc['Tree3'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
df_genes['Diff_gene_Tree1'] = GLS_Tree1
df_genes['Diff_gene_Tree2'] = GLS_Tree2
df_genes['Diff_gene_Tree3'] = GLS_Tree3

df_genes['Diff_gene'] = np.absolute(df_genes['Diff_gene_Tree1'].values - df_genes['Diff_gene_Tree2'].values) + np.absolute(df_genes['Diff_gene_Tree1'].values - df_genes['Diff_gene_Tree3'].values) + np.absolute(df_genes['Diff_gene_Tree2'].values - df_genes['Diff_gene_Tree3'].values)
df_genes['Support']= (df_genes[['Diff_gene_Tree1','Diff_gene_Tree2','Diff_gene_Tree3']].idxmax(axis=1)).str.split(pat='_', expand=True)[2]        #extract just tree name
#df_genes

timestr = time.strftime("_%Y_%m_%d-%H_%M_%S")
df_genes.to_csv(file_likelihood+timestr+'_genelk.csv', float_format='%.5f', header = True, index=False) #rounds nr to 2 decimals

# sorts the log likelihood per gene 
df_sort_GLS2 = df_genes.iloc[(-df_genes['Diff_gene']).argsort()]
#df_sort_GLS2.head()

######################## which graph gets printed ########################
# visualize log-likelihood per gene, sorted and one side of the x-axis = "-gs", "--style1_plot_gene_sort"
if (style1_sorted == False and style2_two_sides == False and style3_mode_site == False) or style1_sorted == True:
    clrs  = ['green' if i == "Tree1" else 'red' if i == 'Tree2' else 'blue' for i in df_sort_GLS2['Support']] # set color 'green' if Tree1 is supported, else set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_sort_GLS2['gene_name'], y=df_sort_GLS2['Diff_gene'],
                marker=dict(color=clrs),
                hovertemplate='Gene: %{x} <br> Signal: %{y:.3f} <extra></extra>'
            )])
    
    fig.update_layout(  
        title='Phylogenetic signal (difference in log-likelihood) per gene', title_x=0.5,
        xaxis_title="Gene",
yaxis_title=u"\u0394"+"GLS",
        annotations = [
        dict(x=1, y=1,
            text="Tree 1", 
            font=dict(size=20, color="green"),
            xref="paper", yref="paper", 
            showarrow=False),
        dict(x=1, y=0.95,
            text="Tree 2",
            font=dict(size=20, color="red"),
            xref="paper", yref="paper", 
            showarrow=False,
            xanchor='auto'
            ),
        dict(x=1, y=0.90,
            text="Tree 3",
            font=dict(size=20, color="blue"),
            xref="paper", yref="paper", 
            showarrow=False,
            xanchor='auto')],
        shapes=[
            dict(
                type= 'line',
                yref= 'y', y0= 0.5, y1= 0.5,
                xref= 'paper', x0= 0, x1= 1,
                line=dict(
                    color="black",
                    width=1,
                    dash="dashdot")
                )],
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    fig.update_yaxes(fixedrange=True)
    fig.show()
################################################################################

# visualize log likelihood per site, no binning = "-ss", "--style3_plot_site"    
if style3_mode_site == True:
    fig = go.Figure()

    # Create visulisation grid
    fig = make_subplots(
      rows=9, cols=1,
      shared_xaxes=True,
      vertical_spacing=0.03,
      specs=[[{"rowspan": 6}],
            [None],
            [None],
            [None],
            [None],
            [None],
            [{"rowspan": 1}],
            [None],
            [{}]]
    )
    
    clrs2  = ['green' if i == "Tree1" else 'red' if i == 'Tree2' else 'blue' for i in df_site_lk.loc['Support']] 
    
    fig.add_trace(
        go.Bar(
                x=df_site_lk.loc['Diff_site'].index, y=df_site_lk.loc['Diff_site'],
                marker=dict(color=clrs2),
                hovertemplate='Position: %{x} <br> Signal: %{y:.3f} <extra></extra>'),
            # Position in graph grid
            row=1, col=1)
    
    fig.update_layout(
        title='Phylogenetic signal (difference in log-likelihood) per site', title_x=0.5,
         annotations = [
            dict(x=1, y=1,
                text="Tree 1", 
                font=dict(size=20, color='green'),
                 xref="paper", yref="paper",
                showarrow=False),
            dict(x=1, y=0.95,
                text="Tree 2",
                font=dict(size=20, color='red'),
                 xref="paper", yref="paper",
                showarrow=False),
            dict(x=1, y=0.90,
                 text="Tree 3",
                 font=dict(size=20, color="blue"),
                 xref="paper", yref="paper", 
                 showarrow=False,
                 xanchor='auto')
            ],
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    
    fig.add_trace(go.Bar(
        x=[df_genes['gene_length'][i] for i in range(len(df_genes))],
        base=[df_genes['gene_begin'][i] for i in range(len(df_genes))],
        y=[df_genes['gene_begin'][0] for i in range(len(df_genes))],
        customdata=df_genes['gene_name'],
        orientation='h',
        hovertemplate='Gene: %{customdata} <extra></extra>', 
        text=df_genes['gene_name'],  
        textposition="inside",        
        marker=dict(color=[f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})' for _ in range(25)]),
            ),  
        row=7, col=1)
    fig.update_layout(uniformtext_minsize=10, uniformtext_mode = 'hide')
    fig.update_yaxes(fixedrange=True)
    fig.show()