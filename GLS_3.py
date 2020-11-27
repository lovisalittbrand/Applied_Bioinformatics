'''Run file like:
python3 GLS_diff-1 file_likelihoods.sitelh file_partition.nexus PLOT_STYLE
where PLOT_STYLE can be any or all of the following: 
-ss     - plots difference in log-likelihoods per site
-gt     - plots difference in log-likelihoods per gene, sorts them and shows them on both sides of the x-axis 
-gs     - plots difference in log-likelihoods per gene, sorts them and shows them on one side of the x-axis (DEFAULT)
'''

import pandas as pd
import plotly.graph_objects as go
import re
import sys
import argparse
import os
#import plotly.io as pio
#pio.renderers.default='browser'	

#file_likelihood = 'T1_and_T2_ML.sitelh'    # the likelihood file
#file_genes_partition = 'tetra.part'       # the partition file


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
'''   
# Read the input agruments to variables
print('file_likelihood =', args.file_likelihood)
print('nexus_partition =', args.nexus_partition)
print('style_sorted =', args.style1_plot_gene_sort)
print('style_two_sides =', args.style2_plot_gene_two_sides)
print('mode_site =', args.style3_plot_site)
'''

file_likelihood = args.file_likelihood
nexus_partition = args.nexus_partition
style1_sorted = args.style1_plot_gene_sort
style2_two_sides = args.style2_plot_gene_two_sides
style3_mode_site=args.style3_plot_site

'''
if (style1_sorted == False and style2_two_sides == False and style3_mode_site == False) or style1_sorted == True: 
    print('graph1')
if style2_two_sides == True:
    print('graph2')
if style3_mode_site == True:
    print('graph3')
'''
######################## CODE ########################

# import log likelihood per site for 2 different topologies: T1 and T2 
with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df_site_lk.loc['Diff_site']=df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values
#print(df_site_lk.head(n=5))

# import nexus partition file and extract gene name, its beginning and end
with open(nexus_partition) as file_genes:
    df_genes = pd.DataFrame(columns=["gene_name", "gene_begin", "gene_end"])
    for line in file_genes:
        if line.startswith((' ')) and 'charset' in line:
            split_line = re.split(r'[\W_]+', line)       #splits on all non-word characters, '+' prevents empty strings 
            #print(split_line)
            df_genes = df_genes.append({"gene_name": split_line[2], "gene_begin": split_line[5], "gene_end": split_line[6]}, ignore_index=True)
        else:
            continue
#print(df_genes.head(n=5))

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS = []  #gene likelihoods
for k in range(0, len(df_genes)):
    GLS.append(sum(df_site_lk.loc['Diff_site'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
#print(GLS)
df_genes['Diff_gene'] = GLS

# sort the log likelihood per gene 
df_sort_GLS = df_genes.sort_values('Diff_gene')

# sorts the log likelihood per gene in absolute values
#df2 = df.iloc[-df.C.abs().argsort()]
df_sort_GLS2 = df_genes.iloc[(-df_genes['Diff_gene'].abs()).argsort()]
#print(df_sort_GLS2.head(n=10))


######################## which graph gets printed ########################
'''
if (style1_sorted and style2_two_sides) == True and (style3_mode_site == False):
    print('graphs1,2')
elif (style1_sorted and style2_two_sides and style3_mode_site) == True:
    print('graphs1,2,3')
elif (style2_two_sides and style3_mode_site) == True and (style1_sorted == False):
    print('graphs2,3')
elif (style1_sorted and style3_mode_site) == True and (style2_two_sides == False):
    print('graphs1,3')
elif (style2_two_sides == True) and (style3_mode_site and style1_sorted) == False:
    print('graph2')
elif (style3_mode_site == True) and (style1_sorted and style2_two_sides) == False:
    print('graph3')
elif (style1_sorted and style2_two_sides and style3_mode_site) == False  or (((style2_two_sides and style3_mode_site) == False) and (style1_sorted == True)):
    print('graph1')
'''

# visualize log-likelihood per gene, sorted and one side of the x-axis = "-gs", "--style1_plot_gene_sort"
if (style1_sorted == False and style2_two_sides == False and style3_mode_site == False) or style1_sorted == True:
    clrs  = ['green' if i >= 0 else 'red' for i in df_sort_GLS2['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_sort_GLS2['gene_name'], y=df_sort_GLS2['Diff_gene'].abs(),
                marker=dict(color=clrs)
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
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    fig.show()
################################################################################

# visualize log-likelihood per gene, sorted and on both sides of the x-axis = "-gt", "--style2_plot_gene_two_sides"
if style2_two_sides == True:
    clrs  = ['green' if i >= 0 else 'red' for i in df_sort_GLS['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_sort_GLS['gene_name'], y=df_sort_GLS['Diff_gene'],
                marker=dict(color=clrs)
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
                showarrow=False
                )],
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    fig.show()  
#######################################################

# visualize log likelihood per site, no binning = "-ss", "--style3_plot_site"    
if style3_mode_site == True:
    clrs2  = ['green' if i >= 0 else 'red' for i in df_site_lk.iloc[2]] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
    
    fig = go.Figure(data=[go.Bar(
                x=df_site_lk.iloc[2].index, y=df_site_lk.iloc[2],
                marker=dict(color=clrs2)
            )])
    
    fig.update_layout(
        title='Phylogenetic signal (difference in log-likelihood) per site', title_x=0.5,
         annotations = [
            dict(x=0, y=0.95,
                text="Tree 1", 
                font=dict(size=20),
                 xref="paper", yref="paper",
                showarrow=False),
            dict(x=0, y=0.95,
                text="Tree 2",
                font=dict(size=20),
                 xref="paper", yref="paper",
                showarrow=False
                )],
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray
    fig.show()
