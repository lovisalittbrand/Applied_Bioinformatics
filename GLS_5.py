'''Run visualization script like this:

    python3 GLS_5.py file_likelihoods.sitelh file_partition.nexus PLOT_STYLE

where PLOT_STYLE can be any or all of the following: 
-gs     - plots difference in log-likelihoods per gene, sorts them and shows them on one side of the x-axis (DEFAULT)
-gt     - plots difference in log-likelihoods per gene, sorts them and shows them on both sides of the x-axis 
-ss     - plots difference in log-likelihoods per site
'''

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
import argparse
import numpy as np

#file_likelihood = 'T1_and_T2_ML.sitelh'    # the likelihood file
#nexus_partition = 'tetra.part'       # the partition file

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

# import log likelihood per site for 2 different topologies: T1 and T2 
with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df_site_lk.loc['Diff_site']=df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values

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

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS = []  #gene likelihoods
for k in range(0, len(df_genes)):
    GLS.append(sum(df_site_lk.loc['Diff_site'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))

df_genes['Diff_gene'] = GLS

# sort the log likelihood per gene 
df_sort_GLS = df_genes.sort_values('Diff_gene')

# sorts the log likelihood per gene in absolute values
df_sort_GLS2 = df_genes.iloc[(-df_genes['Diff_gene'].abs()).argsort()]

######################## which graph gets printed ########################
# visualize log-likelihood per gene, sorted and one side of the x-axis = "-gs", "--style1_plot_gene_sort"
if (style1_sorted == False and style2_two_sides == False and style3_mode_site == False) or style1_sorted == True:
    clrs  = ['green' if i >= 0 else 'red' for i in df_sort_GLS2['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_sort_GLS2['gene_name'], y=df_sort_GLS2['Diff_gene'].abs(),
                marker=dict(color=clrs),
                hovertemplate='Gene: %{x} <br> Signal: %{y:.3f} <extra></extra>'
            )])
    
    fig.update_layout(  
        title='Phylogenetic signal (difference in log-likelihood) per gene', title_x=0.5,
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
                )],
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

# visualize log-likelihood per gene, sorted and on both sides of the x-axis = "-gt", "--style2_plot_gene_two_sides"
if style2_two_sides == True:
    clrs  = ['green' if i >= 0 else 'red' for i in df_sort_GLS['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_sort_GLS['gene_name'], y=df_sort_GLS['Diff_gene'],
                marker=dict(color=clrs),
                hovertemplate='Gene: %{x} <br> Signal: %{y:.3f} <extra></extra>'
            )])
    
    fig.update_layout(
        title='Phylogenetic signal (difference in log-likelihood) per gene', title_x=0.5,
        annotations = [
            dict(x=0, y=1,
                text="Tree 1", 
                font=dict(size=20, color='green'),
                xref="paper", yref="paper",
                showarrow=False),
            dict(x=0, y=0.95,
                text="Tree 2",
                font=dict(size=20, color='red'),
                xref="paper", yref="paper",
                showarrow=False)],
        shapes=[
            dict(
                type= 'line',
                yref= 'y', y0= 0.5, y1= 0.5,
                xref= 'paper', x0= 0, x1= 1,
                line=dict(
                    color="gray",
                    width=1,
                    dash="dashdot")),
            dict(
                type= 'line',
                yref= 'y', y0= -0.5, y1= -0.5,
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
#######################################################

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
    
    clrs2  = ['green' if i >= 0 else 'red' for i in df_site_lk.loc['Diff_site']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
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
            dict(x=0, y=1,
                text="Tree 1", 
                font=dict(size=20, color='green'),
                 xref="paper", yref="paper",
                showarrow=False),
            dict(x=0, y=0.95,
                text="Tree 2",
                font=dict(size=20, color='red'),
                 xref="paper", yref="paper",
                showarrow=False
                )],
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    
    fig.add_trace(go.Bar(
        x=df_genes['gene_length'],
        #y = df_genes.loc[:, row].values,
        y = [0 for i in range(len(df_genes))],
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