# Tool to visualise the base frequency from a MSA
# input like this: python 3 align_nis.py MSA.fasta grouping_size DNA_or_AA

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

filename = sys.argv[1]
window = int(sys.argv[2])

# Check if data is DNA 
if sys.argv[3] == "DNA":
  # Count number of sequences
  no_seq = 0
  # Reads alignment file and counts bases in for each posisiton and creates a matrix containing the data
  with open(filename, 'r') as file:
    for line in file:
      if not (line.startswith('>')):
        if no_seq == 0:
          count_matrix = np.zeros((7, math.ceil(len(line)/window)))
        no_seq += 1
        for j in range(0, len(line)):
          i = math.floor(j/window)
          base = line[j]
          if (base.lower() == "a"):
            count_matrix[0, i] += 1
          elif (base.lower() == "t"):
            count_matrix[1, i] += 1
          elif (base.lower() == "c"):
            count_matrix[2, i] += 1
          elif (base.lower() == "g"):
            count_matrix[3, i] += 1
          elif (base.lower() == "-"):
            count_matrix[4, i] += 1
          elif (base.lower() == "n"):
            count_matrix[5, i] += 1
          else:
            count_matrix[6, i] += 1
  
  # Normalizes the data and creates dataframe of the data.
  freq_matrix = count_matrix/(no_seq*window)

  raw_data = {'A': freq_matrix[0,:], 'T': freq_matrix[1,:],'C': freq_matrix[2,:],'G': freq_matrix[3,:],'-': freq_matrix[4,:]
              ,'n': freq_matrix[5,:],'unknown': freq_matrix[6,:]}
  df = pd.DataFrame(raw_data)

  combined_data = {'A_T': df['A']+df['T'], 'C_G': df['C']+df['G'], 'u_n': df['n']+df['unknown'], 'gap': df['-']}
  df_combined = pd.DataFrame(combined_data)

  # Visualise the data using plotly 
  fig = go.Figure()

  for column in df_combined.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df_combined.index*window,
              y = df_combined[column],
              name = column
          )
      )

  fig.update_layout(barmode='stack')
  fig.update_layout(
      updatemenus=[go.layout.Updatemenu(
          active=0,
          buttons=list(
              [dict(label = 'All',
                    method = 'update',
                    args = [{'visible': [True, True, True, True]},
                            {'barmode':'stack', 'title': 'All',
                            'showlegend':True,}
                            ]),
              dict(label = 'A and T',
                    method = 'update',
                    args = [{'visible': [True, False, False, False]},
                            {'barmode':'group', 'title': 'A_T',
                            'showlegend':True, 'marker_color': 'rgb(26, 118, 255)'}]),
              dict(label = 'C and G',
                    method = 'update',
                    args = [{'visible': [False, True, False, False]},
                            {'barmode':'group', 'title': 'C_G',
                            'showlegend':True}]),
              dict(label = 'N and Unknown',
                    method = 'update',
                    args = [{'visible': [False, False, True, False]},
                            {'barmode':'group', 'title': 'u_n',
                            'showlegend':True}]),
              dict(label = 'Gap',
                    method = 'update',
                    args = [{'visible': [False, False, False, True]},
                            {'barmode':'group', 'title': '-',
                            'showlegend':True}])
              ])
          )
      ])
  fig.update_yaxes(range=[0, 1])

  fig.show()

# Check if AA
elif sys.argv[3] == "AA":
  no_seq = 0
  # Reads alignment file and counts aa in for each posisiton and creates a matrix containing the data
  aa_list = ["a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "-"]
  with open(filename, 'r') as file:
    for line in file:
      if not (line.startswith('>')):
        if no_seq == 0:
          count_matrix = np.zeros((21, math.ceil(len(line)/window)))
        no_seq += 1
        for j in range(0, len(line)):
          i = math.floor(j/window)
          obs_aa = line[j]
          for aa in aa_list:
            if obs_aa.lower() == aa:
              count_matrix[aa_list.index(aa), i] += 1
              break
  
  # Normalizes the data and creates a dataframe of it
  freq_matrix = count_matrix/(no_seq*window)
  aa_dict = dict()
  for i in range(0, len(aa_list)):
    aa_dict[aa_list[i].upper()] = freq_matrix[i,:]

  df = pd.DataFrame(aa_dict)

  # Visualise the data.
  fig = go.Figure()

  fig.update_layout(barmode='stack')

  for column in df.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df.index*window,
              y = df[column],
              name = column
          )
      )
  
  buttons = [dict(label = 'All',
                    method = 'update',
                    args = [{'visible': [True for i in range(len(aa_list))]},
                            {'barmode':'stack', 'title': 'All',
                            'showlegend':True,}
                            ])]

  for aa in aa_list:
    false_list =  [False for i in range(len(aa_list))] 
    false_list[aa_list.index(aa)] = True
    button = [dict(label = aa.upper(),
                    method = 'update',
                    args = [{'visible': false_list},
                            {'barmode':'group', 'title': aa.upper(),
                            'showlegend':True}])]
    buttons = buttons + button

  fig.update_layout(barmode='stack')
  fig.update_layout(
      updatemenus=[go.layout.Updatemenu(
          active=0,
          buttons=list(
              buttons
              )
          )
      ])
  fig.update_yaxes(range=[0, 1])

  fig.show()

else:
  print("Must specify data as DNA or AA")