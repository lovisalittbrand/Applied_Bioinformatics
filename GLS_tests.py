import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
import argparse
import numpy as np
import plotly.io as pio
pio.renderers.default='browser'


file_likelihood = 'T1_and_T2_ML.sitelh'    # the likelihood file
nexus_partition = 'tetra.part'       # the partition file


######################## CODE ########################

# import log likelihood per site for 2 different topologies: T1 and T2 
with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df_site_lk.loc['Diff_site']=np.absolute(df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values)
df_site_lk.loc['Support']= df_site_lk.loc[['Tree1','Tree2']].idxmax(axis=0)
df_site_lk.head()

df_site_lk.loc['Diff_site'][df_site_lk.loc['Support'] =='Tree2'] = df_site_lk.loc['Diff_site']*(-1)





#export the dataframe with genomic position and difference in site likelihood to .csv file
out=df_site_lk.T
out = pd.DataFrame(out).reset_index()
out.columns = ['Genomic_position', 'Tree1', 'Tree2', 'Diff_site', 'Support']
out.to_csv(r'genomePosition_siteLikelihood.csv', float_format='%.5f', header = True, index=False) #rounds nr to 2 decimals

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
print(df_genes)

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS_Tree1 = []  #gene likelihoods
GLS_Tree2 = []

for k in range(0, len(df_genes)):
    GLS_Tree1.append(sum(df_site_lk.loc['Tree1'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
    GLS_Tree2.append(sum(df_site_lk.loc['Tree2'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
df_genes['Diff_gene_Tree1'] = GLS_Tree1
df_genes['Diff_gene_Tree2'] = GLS_Tree2
df_genes['Diff_gene'] = np.absolute(df_genes['Diff_gene_Tree1'].values - df_genes['Diff_gene_Tree2'].values)
df_genes['Support']= (df_genes[['Diff_gene_Tree1','Diff_gene_Tree2']].idxmax(axis=1)).str.split(pat='_', expand=True)[2]        #extract just tree name
#df_genes

#########################


######################## which graph gets printed ########################
# visualize log likelihood per site, no binning

fig = go.Figure()

# Create visulisation grid
fig = make_subplots(
  rows=10, cols=1,
  shared_xaxes=True,
  vertical_spacing=0.03,
  specs=[[{"rowspan": 6}],
        [None],
        [None],
        [None],
        [None],
        [None],
        [None],
        [{"rowspan": 1}],
        [None],
        [{}]]
)

clrs2  = ['green' if i == 'Tree1' else 'red' for i in df_site_lk.loc['Support']] # set color 'green' tree1 is supported, else set color 'red'


fig.add_trace(
    go.Bar(
           x=[i  for j in range(len(df_genes)) for i in range(df_genes['gene_begin'][j], df_genes['gene_end'][j]+1)], 
            y=[df_site_lk.loc['Diff_site'][i]  for j in range(len(df_genes)) for i in range(df_genes['gene_begin'][j], df_genes['gene_end'][j]+1)],
            marker=dict(color=clrs2),
            hovertemplate='Position: %{x} <br> Signal: %{y:.3f} <extra></extra>'),
        # Position in graph grid
        row=1, col=1)

fig.update_layout(
    title='Phylogenetic signal (difference in log-likelihood) per site', title_x=0.5,
    xaxis_title="Position",
    yaxis_title=u"\u0394"+"SLS",
    #range=[df_genes['gene_begin'][0], 'gene_begin'][]],
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
    #x=[df_genes['gene_length'] for i in range(len(df_genes))],
    #y = df_genes.loc[:, row].values,
    #x=[i  for j in range(len(df_genes)) for i in range(df_genes['gene_begin'][j], df_genes['gene_end'][j]+1)], 
    #y=[df_site_lk.loc['Diff_site'][i]  for j in range(len(df_genes)) for i in range(df_genes['gene_begin'][j], df_genes['gene_end'][j]+1)],
    #y = [df_genes['gene_begin'] for i in range(len(df_genes))],
     x=df_genes['gene_length'],
     y = df_genes['gene_begin'],
    customdata=df_genes['gene_name'],
    orientation='h',
    hovertemplate='Gene: %{customdata} <extra></extra>', 
    text=df_genes['gene_name'],  
    textposition="inside",        
    marker=dict(color=[f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})' for _ in range(len(df_genes))]),
        ),  
    row=8, col=1)
fig.update_layout( xaxis2_title="Gene",
    uniformtext_minsize=10, uniformtext_mode = 'hide')
fig.update_yaxes(fixedrange=True)
fig.show()