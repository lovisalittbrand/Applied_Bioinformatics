import pandas as pd
import plotly.graph_objects as go

with open('turtle_site_lk.sitelh') as f:
    df = pd.read_table(f, delim_whitespace=True, header=None, skiprows=1, index_col=0, lineterminator='\n')
'''
print(f'Shape is: {df.shape}')
print(f'The types are: {df.dtypes}') 
print(f'Information about the df: {df.info}')
'''
df.loc['Diff']=df.loc['Tree1'].values - df.loc['Tree2'].values
'''
print(f'Shape is: {df.shape}')
print(f'The types are: {df.dtypes}') 
print(f'Information about the df: {df.info}')
df.loc['Diff'].describe()
'''
# no binning - display each number per position
df2=df.loc['Diff'][0:50]        # I took out 50 positions for displaying; did not work for me with all 187026 positions

# set color 'green' if T1-T2 >= 0, else if T1-T2 < 0 set color 'red'
clrs  = ['green' if i >= 0 else 'red' for i in df2]     

fig = go.Figure(data=[go.Bar(
            x=df2.index, y=df2,
            marker=dict(color=clrs)
        )])
fig.show()
