import pandas as pd
import plotly.graph_objects as go

with open('turtle_site_lk.sitelh') as f:
    df = pd.read_table(f, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')

# calculate the T1-T2 difference between the log likelihoods for each position
df.loc['Diff']=df.loc['Tree1'].values - df.loc['Tree2'].values

# no binning - display each number per position
####################################################
df2=df.loc['Diff'][0:1000]        # I took out 1000 positions for displaying

# set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
clrs  = ['green' if i >= 0 else 'red' for i in df2]     

fig = go.Figure(data=[go.Bar(
            x=df2.index, y=df2,
            marker=dict(color=clrs)
        )])
fig.show()
#####################################################
# for all 187026 positions - you can tell it worked if you change the background from gray to white:
clrs2  = ['green' if i >= 0 else 'red' for i in df.iloc[2]]

fig = go.Figure(data=[go.Bar(
            x=df.iloc[2].index, y=df.iloc[2],
            marker=dict(color=clrs2)
        )])
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
fig.show()
