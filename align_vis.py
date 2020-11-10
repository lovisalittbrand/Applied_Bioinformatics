# Tool to visualise the base frequency from a MSA
# input like this: python3 align_nis.py MSA.fasta window DNA_or_AA -pf partition_file
# MSA.fasta: The alignment file
# window: Size of window that calculates the frequency
# DNA_or_AA: Is the seaquence DNA or amin acid. DNA for DNA and AA for amino acid
# partition_file: The partition file in nexus format. Optional

import argparse
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re


def parseArguments():
    parser = argparse.ArgumentParser()
    # Mandatory arguments
    parser.add_argument("alignemnt_file", help="Alignment file", type=str)
    parser.add_argument("window", help="Window size", type=int)
    parser.add_argument("data", help="DNA or AA", type=str)
    # Optional arguments
    parser.add_argument("-pf", "--partfile", help="Partitionfile in nexus format", type=str, default=None)
    # Version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args
args = parseArguments()

filename = args.alignemnt_file
window = args.window
partition_file = args.partfile
data_type = args.data

# Check if data is DNA 
if data_type == "DNA":
  # Count number of sequences
  no_seq = 0
  seq_length = 0
  # Reads alignment file and counts bases in for each posisiton and creates a matrix containing the data
  with open(filename, 'r') as file:
    for line in file:
      if not (line.startswith('>')):
        if no_seq == 0:
          seq_length = len(line)
          count_matrix = np.zeros((7, math.ceil(seq_length/window)))
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

  gene_data = dict()
  if (partition_file != None):
    with open(partition_file, 'r') as file:
      for line in file:
        if (line.startswith('charset')):
          m = re.match(r"charset (.+) = (\d+)-(\d+);", line)
          gene_data[m[1]] = [int(m[3])-int(m[2])+1]
  else:
    gene_data = {'No genes specified': [seq_length]}
  
  df_gene_info = pd.DataFrame(gene_data)


  # Visualise the data using plotly 
  fig = go.Figure()

  fig = make_subplots(
    rows=6, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.03,
    specs=[[{"rowspan": 5}],
          [None],
          [None],
          [None],
          [None],
           [{}]]
  )

  base_colors = ['red', 'blue', 'green', 'gray']
  l = 0
  for column in df_combined.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df_combined.index*window +1,
              y = df_combined[column],
              name = column,
              marker_color = base_colors[l]
          ),
          row=1, col=1
      )
      l += 1

  for column in df_gene_info.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df_gene_info[column],
              y = [0 for i in range(len(gene_data))] ,
              name = column,
              text=column,
              textposition='auto',
              textfont_color="white",
              showlegend=False,
              orientation='h'
          ),
          row=6, col=1
      )

  gene_plot = [True for i in range(len(gene_data))] 
  fig.update_layout(barmode='stack')
  fig.update_layout(
      updatemenus=[go.layout.Updatemenu(
          active=0,
          buttons=list(
              [dict(label = 'All',
                    method = 'update',
                    args = [{'visible': [True, True, True, True] + gene_plot},
                            {'barmode':'stack', 'title': 'All',
                            'showlegend':True,}
                            ]),
              dict(label = 'A and T',
                    method = 'update',
                    args = [{'visible': [True, False, False, False] + gene_plot},
                            {'barmode':'stack', 'title': 'A_T',
                            'showlegend':True, 'marker_color': 'rgb(26, 118, 255)'}]),
              dict(label = 'C and G',
                    method = 'update',
                    args = [{'visible': [False, True, False, False] + gene_plot},
                            {'barmode':'stack', 'title': 'C_G',
                            'showlegend':True}]),
              dict(label = 'N and Unknown',
                    method = 'update',
                    args = [{'visible': [False, False, True, False] + gene_plot},
                            {'barmode':'stack', 'title': 'u_n',
                            'showlegend':True}]),
              dict(label = 'Gap',
                    method = 'update',
                    args = [{'visible': [False, False, False, True] + gene_plot},
                            {'barmode':'stack', 'title': '-',
                            'showlegend':True}])
              ])
          )
      ])
  fig.update_yaxes(range=[0, 1], row=1, col=1)
  fig.update_yaxes(range=[0, 0.4], row=2, col=1)
  fig.update_yaxes(fixedrange=True)
  fig.update_xaxes(range=[0, seq_length])
  fig.update_xaxes(showticklabels=True, row=2, col=1)

  fig.update_layout(
    xaxis2=dict(
        rangeslider=dict( 
            visible=True 
        )
    )
) 

  fig.show()

# Check if AA
elif data_type == "AA":
  no_seq = 0
  seq_length = 0
  # Reads alignment file and counts aa in for each posisiton and creates a matrix containing the data
  aa_list = ["a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w", "y", "v", "-"]
  with open(filename, 'r') as file:
    for line in file:
      if not (line.startswith('>')):
        if no_seq == 0:
          seq_length = len(line)
          count_matrix = np.zeros((21, math.ceil((seq_length/window))))
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

  gene_data = dict()
  if (partition_file != None):
    with open(partition_file, 'r') as file:
      for line in file:
        if (line.startswith('charset')):
          m = re.match(r"charset (.+) = (\d+)-(\d+);", line)
          gene_data[m[1]] = [int(m[3])-int(m[2])+1]
  else:
    gene_data = {'No genes specified': [seq_length]}
  
  df_gene_info = pd.DataFrame(gene_data)

  # Visualise the data.
  fig = go.Figure()

  fig = make_subplots(
    rows=6, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.03,
    specs=[[{"rowspan": 5}],
          [None],
          [None],
          [None],
          [None],
           [{}]]
  )

  fig.update_layout(barmode='stack')

  for column in df.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df.index*window,
              y = df[column],
              name = column
          ),
          row=1, col=1
      )

  for column in df_gene_info.columns.to_list():
      fig.add_trace(
          go.Bar(
              x = df_gene_info[column],
              y = [0 for i in range(len(gene_data))] ,
              name = column,
              text=column,
              textposition='auto',
              textfont_color="white",
              showlegend=False,
              orientation='h'
          ),
          row=6, col=1
      )

  gene_plot = [True for i in range(len(gene_data))] 
  buttons = [dict(label = 'All',
                    method = 'update',
                    args = [{'visible': [True for i in range(len(aa_list))] + gene_plot},
                            {'barmode':'stack', 'title': 'All',
                            'showlegend':True,}
                            ])]


  for aa in aa_list:
    false_list =  [False for i in range(len(aa_list))] 
    false_list[aa_list.index(aa)] = True
    button = [dict(label = aa.upper(),
                    method = 'update',
                    args = [{'visible': false_list + gene_plot},
                            {'barmode':'stack', 'title': aa.upper(),
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

  fig.update_yaxes(range=[0, 1], row=1, col=1)
  fig.update_yaxes(range=[0, 0.4], row=2, col=1)
  fig.update_yaxes(fixedrange=True)
  fig.update_xaxes(range=[0, seq_length])
  fig.update_xaxes(showticklabels=True, row=2, col=1)

  fig.update_layout(
    xaxis2=dict(
        rangeslider=dict( 
            visible=True 
        )
    )
) 

  fig.show()

else:
  print("Must specify data as DNA or AA")