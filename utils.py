from sklearn.manifold import MDS
import numpy as np

def fake_points_MDS(similarity_matrix,n_components):
    # perform MDS
    mds = MDS(n_components=n_components, dissimilarity='precomputed')
    X_mds = mds.fit_transform(similarity_matrix)
    return X_mds