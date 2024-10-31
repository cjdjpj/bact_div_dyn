import numpy as np
import msprime
import seaborn as sns
import matplotlib.pyplot as plt
import hopkins
from sklearn.manifold import MDS
from sklearn.decomposition import PCA

l = 1000  # number of genes
t = 1  # tract length
nsample = 500  # the number of genomes sampled
mu = 3e-6 # mutation rate
r_m = 0.003

r = r_m * mu

a = 1.00001
ts = msprime.sim_ancestry(nsample,
                          model=msprime.BetaCoalescent(alpha=a),
                          ploidy=1,
                          sequence_length=l,
                          gene_conversion_rate=r,
                          gene_conversion_tract_length=t,
                          discrete_genome=False)

mts = msprime.sim_mutations(ts, rate=mu)

### compute distance list and distance matrix
pairs = [(i, j) for i in range(nsample) for j in range(i)]
distance_list = mts.diversity(pairs, mode='site')
distance_matrix = np.zeros((nsample,nsample))

for pair_index, distance in enumerate(distance_list):
    (i,j) = pairs[pair_index]
    distance_matrix[i, j] = distance
    distance_matrix[j, i] = distance

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

### MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
mds_coords = mds.fit_transform(distance_matrix)

plt.figure(figsize = (9,9))
sns.scatterplot(x=mds_coords[:, 0], y=mds_coords[:, 1])
plt.xlabel("MDS 1")
plt.ylabel("MDS 2")

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, bins = np.linspace(0, 0.1, 100));
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
plt.show()

### HOPKINS STATISTIC
print("hopkins: " + str(hopkins.hopkins_statistic(distance_matrix)))

### COEFFICIENT OF VARIATION
print("cv: " + str(np.std(distance_list) / np.mean(distance_list)))

### KMEANS
# n_clusters = 4
# kmeans = KMeans(n_clusters=n_clusters, random_state=42)
# kmeans.fit(pairwise_dist)
#
# print(kmeans.inertia_)
