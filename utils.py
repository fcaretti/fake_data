from sklearn.manifold import MDS
import numpy as np

fake_points_MDS(similarity_matrix):
    # perform MDS
    mds = MDS(n_components=N, dissimilarity='precomputed')
    X_mds = mds.fit_transform(matrix)
    return X_mds