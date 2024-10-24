import tskit
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import hopkins
from sklearn.manifold import MDS
from sklearn.decomposition import PCA

input_path = "runs/r9"
mts = tskit.load(input_path)
with open(input_path + "_ds", 'rb') as file:
    distance_list = pickle.load(file)
    distance_matrix = pickle.load(file)

nsample = mts.num_samples

### MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
mds_coords = mds.fit_transform(distance_matrix)

plt.figure(figsize = (9,9))
sns.scatterplot(x=mds_coords[:, 0], y=mds_coords[:, 1])
plt.xlabel("MDS 1")
plt.ylabel("MDS 2")
plt.show()

### PCA
genotype_matrix = np.array(list(mts.genotype_matrix()))
pca = PCA()
pca_result = pca.fit_transform(genotype_matrix)

plt.figure(figsize=(9, 9))
plt.scatter(pca_result[:, 0], pca_result[:, 1])
plt.title("PCA of Genotype Matrix")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.show()

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, bins = np.linspace(0, 1, 100));
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
plt.show()

### HOPKINS STATISTIC
print(hopkins.hopkins_statistic(distance_matrix))

### COEFFICIENT OF VARIATION
print(np.std(distance_list) / np.mean(distance_list))

### KMEANS
# n_clusters = 4
# kmeans = KMeans(n_clusters=n_clusters, random_state=42)
# kmeans.fit(pairwise_dist)
#
# print(kmeans.inertia_)
