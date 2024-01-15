### Python script to run CN data (features) through the CNV_chromothripsis model
# Paper: https://doi.org/10.1093/bioinformatics/btad422
# GitHub: https://github.com/luvyfdawnYu/CNV_chromothripsis

# Import the required classes and libraries
import itertools
import torch
import pandas as pd
import numpy as np
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from torch_geometric.loader import DataLoader as GraphDataLoader
from GECNVNet import GECNVNet
from utils import prepare_data, to_categorical
from create_graph import TrainDataset, ValDataset, TestDataset

# Read in the data and prepare it for model ingestion; re-scale between 0 and 1 and standardize
df = pd.read_csv("./cn_calls/30kb_aCN_features.csv", index_col = 0).fillna(0)
df_org = df.iloc[:, :-1]
df_pred = df['chromothripsis'] # This is a dummy ground truth value; the actual chromothripsis state of each sample is unkown

# Re-scale between 0 and 1
mm = MinMaxScaler((0,1))
data_org = np.array(df_org.values, dtype = np.float32)
data_org = mm.fit_transform(data_org)

# Standardize
ss = StandardScaler(with_mean = True, with_std = True)
data_org = ss.fit_transform(data_org)

# Restore back to dataframe and array of labels
df_org = pd.DataFrame(data = data_org, columns = df_org.columns)
data_all = prepare_data(df_org)

label_all = np.array(df_pred.values, dtype = np.float32)

# Deserialize and initialize the models
paths = ["./ml_models/GECNVNet_Split1/model_weights/", "./ml_models/GECNVNet_Split2/model_weights/", "./ml_models/GECNVNet_Split3/model_weights/", 
"./ml_models/GECNVNet_Split4/model_weights/", "./ml_models/GECNVNet_Split5/model_weights/"]
model_paths = []

for path in paths:
	for i in range(10):
		model_paths.append(path + "model" + str(i) + ".pkl")

test_data = TestDataset(root = "./graph_data/cn_sigs_all_models", feat = data_all, label = to_categorical(label_all, 2))
data_iter = GraphDataLoader(test_data, batch_size = 64)
prob_list = []

for path in model_paths:
	model = GECNVNet([16, 32], 32, [1024,24*32,16*32,8*32], [24,16,8,4], 32*4, 32*4, 0.5).to("cpu")
	model.load_state_dict(torch.load(path))
	model.eval()

	# Generate predictions
	out_list = []

	with torch.no_grad():
		for data in data_iter:
			out = model(data)
			out_prob=F.softmax(out, dim = -1)
			out_list.append(out_prob)

	out_prob = torch.cat(out_list, 0)
	prob_list.append(out_prob.detach().cpu().numpy())

# Get the average probability across the models
ens_pred = np.zeros((label_all.shape[0], 2))
for i in prob_list:
	ens_pred += i

ens_pred = ens_pred / len(model_paths)

df_prob = pd.DataFrame()
df_prob["sample"] = df.index
df_prob["pred_prob"] = ens_pred[:, 1] # Second column contains the probabiltiy of chromothripsis
df_prob.to_csv(path_or_buf = "./output/30kb_aCN_CNV_chromothripsis_predictions_all_models.csv", sep = ",", header = True)
