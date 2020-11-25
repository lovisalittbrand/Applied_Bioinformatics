import pandas as pd
import plotly.graph_objects as go
import sys

filename = sys.argv[1]

with open(filename) as f:
    df = pd.read_table(f, delim_whitespace=True, header=None, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df.loc['Diff']=df.loc['Tree1'].values - df.loc['Tree2'].values

# no binning - display each number per position
#####################################################

# for all 187026 positions - you can tell it worked if you change the background from gray to white:
clrs2  = ['green' if i >= 0 else 'red' for i in df.iloc[2]] # set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'

fig = go.Figure(data=[go.Bar(
            x=df.iloc[2].index, y=df.iloc[2],
            marker=dict(color=clrs2)
        )])
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)') # remove this line if you want your background to be gray
fig.show()
