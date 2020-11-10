
import plotly.express as px
import plotly.graph_objects as go
import sys
import math
import numpy as np
import pandas as pd

filename = sys.argv[1]
window = int(sys.argv[2])

# Check if data is DNA
if sys.argv[3] == "DNA":
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
              print(aa)
              break

 # Normalizes the data and creates dataframe of the data.
freq_matrix = count_matrix/(no_seq*window)


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

    return entropy_per_site

no_sites = np.arange(len(freq_matrix[0,:]))
entropy_vector = calc_entropy(freq_matrix)
entropy_mod_vector = {'Entropy': entropy_vector, 'Site':no_sites}
entropy_df = pd.DataFrame(entropy_mod_vector)


fig = go.Figure()
fig = px.line(entropy_df, y="Entropy", x="Site", title="Entropy across alignment", range_y = [0,1])
fig.show()

#From nextstrain: return sum([v * math.log(v+1E-20) for v in vals]) * -1 / math.log(len(vals))
