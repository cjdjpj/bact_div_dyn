import json
import tskit
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import hopkins
from sklearn.manifold import MDS
from sklearn.decomposition import PCA

###
save_fig = False
run_index = "r002"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

mts = tskit.load(input_path)
nsample = mts.num_samples
with open(input_path + "_ds", "rb") as file:
    distance_list = pickle.load(file)
    distance_matrix = pickle.load(file)

print("Average pi:" + str(sum(distance_list)/len(distance_list)))
### MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
mds_coords = mds.fit_transform(distance_matrix)

plt.figure(figsize = (9,9))
sns.scatterplot(x=mds_coords[:, 0], y=mds_coords[:, 1])
plt.xlabel("MDS 1")
plt.ylabel("MDS 2")
if save_fig == True:
    plt.savefig(run_index + "b.png", dpi=300)
else:
    plt.show()


### PCA
# genotype_matrix = np.array(list(mts.genotype_matrix()))
# pca = PCA()
# pca_result = pca.fit_transform(genotype_matrix)
#
# plt.figure(figsize=(9, 9))
# plt.scatter(pca_result[:, 0], pca_result[:, 1])
# plt.title("PCA of Genotype Matrix")
# plt.xlabel("Principal Component 1")
# plt.ylabel("Principal Component 2")
# plt.show()

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, bins = np.linspace(0, 0.15, 100));
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
if save_fig == True:
    plt.savefig(run_index + "a.png", dpi=300)
else:
    plt.show()

### FRACTION OF IDENTICAL BLOCKS
# pairs = [(i, j) for i in range(nsample) for j in range(i)]
# prop_identical_blk = []
# for (i, j) in pairs:
#     identical_sites = 0
#     total_sites = mts.num_sites
#
#     for variant in mts.variants():
#         allele_i = variant.genotypes[i]
#         allele_j = variant.genotypes[j]
#
#         if allele_i == allele_j:
#             identical_sites += 1
#
#     fraction_identical = identical_sites / total_sites
#     prop_identical_blk.append(fraction_identical)
#
# print(prop_identical_blk)

### HOPKINS STATISTIC
# print("hopkins: " + str(hopkins.hopkins_statistic(distance_matrix)))

### COEFFICIENT OF VARIATION
# print("cv: " + str(np.std(distance_list) / np.mean(distance_list)))

### KMEANS
# n_clusters = 4
# kmeans = KMeans(n_clusters=n_clusters, random_state=42)
# kmeans.fit(pairwise_dist)
#
# print(kmeans.inertia_)
