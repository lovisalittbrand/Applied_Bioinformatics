import pandas as pd

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
print(df.index)

import plotly.express as px

fig = px.bar(x=df.loc['Diff'].index, y=df.loc['Diff'])
fig.show()
