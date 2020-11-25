import pandas as pd
import plotly.graph_objects as go
import re
import plotly.io as pio
pio.renderers.default='browser'	

file_likelihood ='turtle_site_lk.sitelh'
file_genes = 'turtle_genes_order.txt'


with open(file_likelihood) as log_likelihood:
    df_site_lk = pd.read_table(log_likelihood, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df_site_lk.loc['Diff_site']=df_site_lk.loc['Tree1'].values - df_site_lk.loc['Tree2'].values

# no binning - display each number per position
'''
# for all positions
clrs2  = ['green' if i >= 0 else 'red' for i in df_site_lk.iloc[2]] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'

fig = go.Figure(data=[go.Bar(
            x=df_site_lk.iloc[2].index, y=df_site_lk.iloc[2],
            marker=dict(color=clrs2)
        )])
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)') # remove this line if you want your background to be gray
fig.show()
'''

# import the file that has the gene limits and extract only the beginning and the end of the genes
with open(file_genes) as file_genes:
    df_genes = pd.DataFrame(columns=["gene_begin", "gene_end"])
    for line in file_genes:
        res = line.split(" = ", 1)[1]
        #gene_beg, gene_end = res.split("-", 1)      #splits on '-'
        gene_beg, gene_end, a = re.split('\W', res)  #splits on all non-word characters
        df_genes = df_genes.append({"gene_begin": gene_beg, "gene_end": gene_end}, ignore_index=True)

#print(sum(df_site_lk.loc['Diff'][1:499]))
#print(sum(df_site_lk.loc['Diff'][int(df_genes['gene_begin'][0]): int(df_genes['gene_end'][0])+1])) #sums items from 1st through stop-1

# calculate the log likelihoods differences per gene and add it as a column to the dataframe with gene data
GLS = []  #gene likelihoods
for k in range(0, len(df_genes)):
    GLS.append(sum(df_site_lk.loc['Diff_site'][int(df_genes['gene_begin'][k]): int(df_genes['gene_end'][k])+1]))
print(GLS)

df_genes['Diff_gene'] = GLS

# visualize gene log likelihood
clrs1  = ['green' if i >= 0 else 'red' for i in df_genes['Diff_gene']] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'

fig = go.Figure(data=[go.Bar(
            x=df_genes['Diff_gene'].index, y=df_genes['Diff_gene'],
            marker=dict(color=clrs1)
        )])
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)') # remove this line if you want your background to be gray
fig.show()