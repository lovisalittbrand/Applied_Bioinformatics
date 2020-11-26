import pandas as pd
import plotly.graph_objects as go
import re
import sys
import argparse
import plotly.io as pio
pio.renderers.default='browser'	

#file_likelihood = sys.argv[1]       # turtle_site_lk.sitelh'    # the likelihood file
file_likelihood = 'T1_and_T2_ML.sitelh'    # the likelihood file
#file_genes = sys.argv[2]            #'turtle_genes_order.txt'       # the partition file
file_genes_partition = 'tetra.part'       # the partition file


########## CODE ########################
#### import file with log likelihood per site for 2 different topologies: T1 and T2, and visualize log likelihood per site ####
with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df_site_lk.loc['Diff_site']=df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values



#### import file with gene metadata (beginning, end) and visualize log likelihood per gene ####

# import the file that has the gene limits and extract the gene name, its beginning and the end
with open(file_genes_partition) as file_genes:
    df_genes = pd.DataFrame(columns=["gene_name", "gene_begin", "gene_end"])
    for line in file_genes:
        if line.startswith((' ')) and 'charset' in line:
            split_line = re.split(r'[\W_]+', line)       #splits on all non-word characters, '+' prevents empty strings 
            print(split_line)
            df_genes = df_genes.append({"gene_name": split_line[2], "gene_begin": split_line[5], "gene_end": split_line[6]}, ignore_index=True)
        else:
            continue
print(df_genes)

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS = []  #gene likelihoods
for k in range(0, len(df_genes)):
    GLS.append(sum(df_site_lk.loc['Diff_site'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
#print(GLS)

df_genes['Diff_gene'] = GLS
#df2 = df.iloc[-df.C.abs().argsort()]



#### sort the log likelihood per gene ####

df_sort_GLS = df_genes.sort_values('Diff_gene')

# visualize sorted gene log likelihood with gene name
clrs  = ['green' if i >= 0 else 'red' for i in df_sort_GLS['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'

fig = go.Figure(data=[go.Bar(
            x=df_sort_GLS['gene_name'], y=df_sort_GLS['Diff_gene'],
            marker=dict(color=clrs)
        )])

fig.update_layout(
    title='Phylogenetic signal (difference in log-likelihood) per gene', title_x=0.5,
    annotations = [
        dict(x=len(df_genes)/2, y=-4,
            text="Tree 1", 
            font=dict(size=20),
            showarrow=False),
        dict(x=len(df_genes)/2, y=4,
            text="Tree 2",
            font=dict(size=20),
            showarrow=False
            )],
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray

fig.show()


#######################################################

#### doesn't sort, prints all log likelihood per gene above the y axis and colorcoded####

# visualize gene log likelihood with gene name
clrs1  = ['green' if i >= 0 else 'red' for i in df_genes['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'

fig = go.Figure(data=[go.Bar(
            x=df_genes['gene_name'], y=df_genes['Diff_gene'].abs(),
            marker=dict(color=clrs1)
        )])
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)') # remove this line if you want your background to be gray
fig.show()

#######################################################

#### sorts the log likelihood per gene in absolute ####

#df2 = df.iloc[-df.C.abs().argsort()]
df_sort_GLS2 = df_genes.iloc[(-df_genes['Diff_gene'].abs()).argsort()]
print(df_sort_GLS2)
# visualize sorted gene log likelihood with gene name
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
            font=dict(size=20, color="red"),
            xref="paper", yref="paper", 
            showarrow=False),
        dict(x=1, y=0.95,
            text="Tree 2",
            font=dict(size=20, color="green"),
            xref="paper", yref="paper", 
            showarrow=False,
            xanchor='auto'
            )],
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=False) # remove 'plot_bgcolor' if you want your background to be gray

fig.show()
