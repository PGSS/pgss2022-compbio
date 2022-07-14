#Reading the file
import pandas as pd
Seq_Analy = pd.read_csv("./Sequence_Analyses.csv")
Seq_Analy.agg({'seq_length' : ['min', 'max', 'mean', 'std']})