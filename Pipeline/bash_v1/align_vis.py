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
import plotly.express as px
from plotly.subplots import make_subplots
import re
import ast


#################################### Functions #####################################

# Parser function that allows to have optional commands
def parseArguments():
    parser = argparse.ArgumentParser()
    # Mandatory arguments
    parser.add_argument("alignemnt_file", help="Alignment file", type=str)
    parser.add_argument("window", help="Window size", type=int)
    parser.add_argument("data", help="DNA or AA", type=str)
    # Optional arguments
    parser.add_argument("-pf", "--partfile", help="Partitionfile in nexus format", type=str, default=None)
    parser.add_argument("-g", "--gapexclude", help="Calculate frequency excluding gaps", type=str2bool, default=False)
    # Version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args

#function to enable boolean inputs in the parser, also handles other inputs than true or false
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Function to calculate entropy from a frequency martix
def calc_entropy(freq_mat):
    aln_len = len(freq_matrix[0,:])
    no_features = len(freq_matrix[:,0])
    entropy_per_site = np.zeros(aln_len)

    for i in range(0,aln_len):
        s = 0
        for j in range(0,no_features):
            s = s + freq_mat[j,i]*math.log(freq_mat[j,i]+1E-20)
        entropy = s * -1 / math.log(no_features)
        entropy_per_site[i] = entropy

    entropy_mod_vector = {'Entropy': entropy_per_site, 'Site':range(0, aln_len)}
    entropy_df = pd.DataFrame(entropy_mod_vector)
    return entropy_df

# Founction to count no of uccurences in a sequence.
def count_chars(filename, char_list):
  no_seq = 0
  seq_length = 0
  # Reads alignment file and counts bases in for each posisiton and creates a matrix containing the data
  # Read alignemnt file
  with open(filename, 'r') as file:
    for line in file:
      # Read only lines with sequence data
      if not (line.startswith('>')):
        # Declarate length and create a zero matrix with correct size in the first loop.
        if no_seq == 0:
          seq_length = len(line)-1
          count_matrix = np.zeros((len(char_list)+1, math.ceil(seq_length/window)))
        no_seq += 1
        # Loop through each position in the sequence
        for j in range(0, len(line)-1):
          # Calculate to what index the count should be grouped to based on the window value. 
          i = math.floor(j/window)
          obs_char = line[j]
          found = False
          # Loop over each possible state in the base_list, if a match is foud add 1 to the correct position in the matrix
          # If no match is found add to the unknown column in the matrix.
          for char in char_list:
            if obs_char.lower() == char.lower():
              count_matrix[char_list.index(char), i] += 1
              found = True
              break
          if found == False:
            count_matrix[len(char_list), i] += 1
  return [count_matrix, seq_length, no_seq]


# Reads a partitionfile ant transform the data to a dataframe for visualisation.
def gene_info_to_dict(partition_file, seq_length):
  gene_data = dict()
  # Check if partition file is awailable 
  if (partition_file != None):
    # Open file and read each line in search for lines that starts with charstet 
    with open(partition_file, 'r') as file:
      for line in file:
        if (line.startswith('charset')):
          # Regex that fetch the name of teh gene and start and stop position and then puts the length and name in a dict
          m = re.match(r"charset (.+) = (\d+)-(\d+);", line)
          gene_data[m[1]] = [int(m[3])-int(m[2])+1]
  # If no partition file is specified, create a gene as long as the sequence 
  else:
    gene_data = {'No genes specified': [seq_length]}
  
  return gene_data

def freq_and_entropy(count_matrix, char_list):
  # Normalizes the data by dividing each cell responding to one position by the sum of all cells on that position.
  freq_matrix = count_matrix/(count_matrix.sum(axis=0)[:,None]).T

  #  Calculate entropy of the freq matrix
  df_entropy = calc_entropy(freq_matrix)

  # Transforms the frequency matrix to a dictionary and then a dataframe, allso adds "u" to the base list so that the loop will parse through every column
  base_dict = dict()
  for i in range(0, len(base_list)):
    base_dict[base_list[i].lower()] = freq_matrix[i,:]
  df_freq = pd.DataFrame(base_dict)
  return [df_freq, df_entropy]

################################## Program #####################################

# Parse the input arguments
args = parseArguments()

# Read the iput agruments to varables
filename = args.alignemnt_file
window = args.window
partition_file = args.partfile
data_type = args.data
gapexclude = args.gapexclude

# Check if data is DNA 
if data_type == "DNA":
  base_list = ["a", "t", "c", "g", "-", "n"]
  count_matrix, seq_length, no_seq = count_chars(filename, base_list)

  base_list = base_list + ["u"]
  # Handles the gapexclude option.
  

  if gapexclude == True:
    # Removes counts for gaps.
    count_matrix_del = np.delete(count_matrix, base_list.index('-'), 0)
    # Removes "-" and adds "u" to the base list.
    base_list.remove("-") 
    # Normalizes the data by dividing each cell responding to one position by the sum of all cells on that position.
    freq_matrix = count_matrix_del/(count_matrix_del.sum(axis=0)[:,None]).T
    
    base_dict = dict()
    for i in range(0, len(base_list)):
      base_dict[base_list[i].lower()] = freq_matrix[i,:]
    df = pd.DataFrame(base_dict)

    combined_data = {'A and T': df['a']+df['t'], 'C and G': df['c']+df['g'], 'N': df['n'], "Unknown": df['u']}
    # Creates a list with names for the different combinations
    combined_bases_list = ['A and T', 'C and G', 'N', 'Unknown']
    df_combined = pd.DataFrame(combined_data)

  else:
    # Normalizes the data by dividing each cell responding to one position by the sum of all cells on that position.
    freq_matrix = count_matrix/(count_matrix.sum(axis=0)[:,None]).T

    # Transforms the frequency matrix to a dictionary and then a dataframe, allso adds "u" to the base list so that the loop will parse through every column
    base_dict = dict()
    for i in range(0, len(base_list)):
      base_dict[base_list[i].lower()] = freq_matrix[i,:]
    df = pd.DataFrame(base_dict)

    # Combines the base data in a meaningfull way.
    combined_data = {'A and T': df['a']+df['t'], 'C and G': df['c']+df['g'], 'N': df['n'], "Unknown": df['u'], 'Gap': df['-']}
    # Creates a list with names for the different combinations
    combined_bases_list = ['A and T', 'C and G', 'N', 'Unknown', 'Gap']
    df_combined = pd.DataFrame(combined_data)

  # Replace zeros to nan so that the wont show up in the plot 
  df_combined = df_combined.replace(0, np.nan)

  # Generate dataframe of gene data that willwe used in the visualisation
  gene_data = gene_info_to_dict(partition_file, seq_length)
  df_gene_info = pd.DataFrame(gene_data)

  # Visualise the data using plotly 
  fig = go.Figure()

  # Create the grid in which the plots will be put.
  fig = make_subplots(
    rows=9, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.03,
    specs=[[{"rowspan": 5}],
          [None],
          [None],
          [None],
          [None],
          [{"rowspan": 3}],
          [None],
          [None],
           [{}]]
  )

  # Define the colors of the bases
  file = open("dna_color.txt", "r")
  contents = file.read()
  color_dict = ast.literal_eval(contents)
  file.close()
  # Loop through the columns and create a barplot for each column and put them in the same position in the grid
  for column in df_combined.columns.to_list():
    if column in color_dict:
      fig.add_trace(
          # Define that it's a bar graph
          go.Bar(
              # position data
              x = df_combined.index*window +1,
              # Base freq data
              y = df_combined[column],
              # Define the name of the data based on the column name
              name = column,
              # Decide the color of bargraph 
              marker_color = color_dict[column],
              # Styling of hover info
              hovertemplate='Base: ' + column + ' Frequency: %{y:.3f} Position: %{x:.3f}',
          ),
          row=1, col=1
      )
    else:
      fig.add_trace(
          go.Bar(
              x = df.index*window+1,
              y = df[column],
              name = column,
              hovertemplate='Base: ' + column + ' Frequency: %{y:.3f} Position: %{x:.3f}',        
          ),
          row=1, col=1
      )

  #  Calculate entropy of the freq matrix
  df_entropy = calc_entropy(freq_matrix)

  # Plot the lineplot of the Entropy 
  fig = fig.add_trace(
          # Define that it's a line graph
          go.Scatter(
            y=df_entropy['Entropy'], 
            x=df_entropy.index*window+1, 
            name = 'Entropy',
            hovertemplate='Entropy: %{y:.3f} <br> Position: %{x:.3f}'),
          row=6, col=1)

  # Same principle as above 
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
              orientation='h',
              hoverinfo='skip',
          ),
          row=9, col=1
      )

######### Optional? ##########

  # Define list with true to use later to decide on their related graphs visibility 
  entropy_plot = [True]
  gene_plot = [True for i in range(len(gene_data))] 
  # Create the first option in the dropdown menu and the features that comes whith that option.
  buttons = [dict(label = 'All',
                    method = 'update',
                    args = [{'visible': [True for i in range(len(base_list))] + entropy_plot + gene_plot},
                            {'barmode':'stack', 'title': 'All',
                            'showlegend':True,}
                            ])]

  # Loop over the positions in the combined base list to decide on which bar plot that shuld be visible.
  for base in combined_bases_list:
    false_list =  [False for i in range(len(combined_bases_list))] 
    false_list[combined_bases_list.index(base)] = True
    # Create a button much like above
    button = [dict(label = base,
                    method = 'update',
                    args = [{'visible': false_list + entropy_plot + gene_plot},
                            {'barmode':'stack', 'title': base,
                            'showlegend':True}])]
    buttons = buttons + button

  
  # Show the dropdown menu.
  fig.update_layout(
      updatemenus=[go.layout.Updatemenu(
          active=0,
          buttons=list(
              buttons
              )
          )
      ])
######### ######### ######### ######### ######### ######### 

  

# Check if AA
elif data_type == "AA":
  aa_list = ["f", "i", "w", "l", "v", "m", "y", "c", "a", "g", "p", "h", "t", "s", "q", "n", "e", "d", "k", "r", "-"]
  count_matrix, seq_length, no_seq = count_chars(filename, aa_list)
  aa_list = aa_list +["u"]

  
   
  if gapexclude == True:
    # Removes counts for gaps.
    count_matrix_del = np.delete(count_matrix, aa_list.index('-'), 0)
    # Normalizes the data by dividing each cell responding to one position by the sum of all cells on that position.
    freq_matrix = count_matrix_del/(count_matrix_del.sum(axis=0)[:,None]).T
    # Removes "-" and adds "u" to the base list.
    aa_list.remove("-") 

    aa_dict = dict()
    for i in range(0, len(aa_list)):
      aa_dict[aa_list[i].upper()] = freq_matrix[i,:]
    df = pd.DataFrame(aa_dict)

    df_combined = pd.DataFrame(aa_dict)


  else:
  # Normalizes the data by dividing each cell in the matrix with the total number of bases for that site/grouping.
    freq_matrix = count_matrix/(count_matrix.sum(axis=0)[:,None]).T

    # Transforms the frequency matrix to a dictionary and then a dataframe
    
    aa_dict = dict()
    for i in range(0, len(aa_list)):
      aa_dict[aa_list[i].upper()] = freq_matrix[i,:]

    df = pd.DataFrame(aa_dict)

  # Replace zeros to nan so that the wont show up in the plot 
  df = df.replace(0, np.nan)

  # Generate dataframe of gene data that willwe used in the visualisation
  gene_data = gene_info_to_dict(partition_file, seq_length)
  df_gene_info = pd.DataFrame(gene_data)

  # Visualise the data.
  fig = go.Figure()

  fig = make_subplots(
    rows=9, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.03,
    specs=[[{"rowspan": 5}],
          [None],
          [None],
          [None],
          [None],
          [{"rowspan": 3}],
          [None],
          [None],
           [{}]]
  )

  file = open("aa_color.txt", "r")
  contents = file.read()
  color_dict = ast.literal_eval(contents)
  file.close()

  # Loop through the columns and create a barplot for each column and put them in the same position in the grid
  for column in df.columns.to_list():
    # If statement to handle aa that are inputed in the aa_list but arent in the aa_color.txt
    if column in color_dict:
      fig.add_trace(
          go.Bar(
              x = df.index*window+1,
              y = df[column],
              name = column,
              marker_color = color_dict[column],   
              hovertemplate='Amino acid:' + column + ' Frequency: %{y:.3f} Position: %{x:.3f}',     
          ),
          row=1, col=1
      )
    else:
      fig.add_trace(
          go.Bar(
              x = df.index*window+1,
              y = df[column],
              name = column,
              hovertemplate='Amino acid:' + column + ' Frequency: %{y:.3f} Position: %{x:.3f}',        
          ),
          row=1, col=1
      )


  df_entropy = calc_entropy(freq_matrix)

  fig = fig.add_trace(
          go.Scatter(
            y=df_entropy['Entropy'], 
            x=df_entropy.index*window+1, 
            name = 'Entropy',
            hovertemplate='Entropy: %{y:.3f} <br> Position: %{x:.3f}'),
          row=6, col=1)

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
              orientation='h',
              hoverinfo='skip',
          ),
          row=9, col=1
      )


######### Optional? ##########
  # gene_plot = [True for i in range(len(gene_data)+1)] 
  # buttons = [dict(label = 'All',
  #                   method = 'update',
  #                   args = [{'visible': [True for i in range(len(aa_list))] + gene_plot},
  #                           {'barmode':'stack', 'title': 'All',
  #                           'showlegend':True,}
  #                           ])]


  # for aa in aa_list:
  #   false_list =  [False for i in range(len(aa_list))] 
  #   false_list[aa_list.index(aa)] = True
  #   button = [dict(label = aa.upper(),
  #                   method = 'update',
  #                   args = [{'visible': false_list + gene_plot},
  #                           {'barmode':'stack', 'title': aa.upper(),
  #                           'showlegend':True}])]
  #   buttons = buttons + button

  # fig.update_layout(barmode='stack')
  # fig.update_layout(
  #     updatemenus=[go.layout.Updatemenu(
  #         active=0,
  #         buttons=list(
  #             buttons
  #             )
  #         )
  #     ])

  #############################################

# Settings for the axis in the visualisation
# General settings
fig.update_yaxes(fixedrange=True)
fig.update_xaxes(range=[0, seq_length])
# Stack the bar plots on each other in the initial visualisation
fig.update_layout(barmode='stack', hovermode='x')
# Frequency graph 
fig.update_yaxes(range=[0, 1], row=1, col=1)
# Entropy graph
fig.update_yaxes(range=[0, 1], row=6, col=1)
# Gene graph
fig.update_yaxes(range=[-0.4, 0.4], row=9, col=1)  
fig.update_yaxes(showticklabels=False, row=9, col=1)
fig.update_xaxes(showticklabels=True, row=9, col=1)

# Create the rangeslider based on the gene graph.
fig.update_layout(
  xaxis3=dict(
      rangeslider=dict( 
          visible=True 
      )
  )
) 

fig.show()