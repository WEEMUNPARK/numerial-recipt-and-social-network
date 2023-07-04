from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import numpy.matlib as mat

def PCA_component(data_Mat, n_comps =3):
    [n, p] = np.shape(data_Mat)  # shape of data

    pca = PCA(n_components = n_comps)  # n_components can be integer or float in (0,1)
    pca.fit(data_Mat)  # fit the model
    print('\n PCA by Scikit-learn:')

    explained_ratios = pca.explained_variance_ratio_.cumsum()
    print(explained_ratios)
    print(pca.explained_variance_)

    print('After PCA transformation, data becomes:')
    PCA_data_Mat = pca.fit_transform(data_Mat)
    # df_pca_data = pd.DataFrame(PCA_data_Mat, index=range(n), columns=['PCA_1', 'PCA_2', 'PCA_3'])
    if n_comps < p:
        U = np.zeros((n, p - n_comps), dtype=float)
        Y = np.append(PCA_data_Mat, U, axis=1)
        df_pca_data = pd.DataFrame(Y, index=range(n))
    else:  # full rank
        df_pca_data = pd.DataFrame(PCA_data_Mat, index=range(n))

    # restore back to raw data
    data_mean = pca.mean_
    PCA_coeff_A = np.transpose(pca.components_)
    new_data_Mat = np.dot(PCA_data_Mat, np.transpose(PCA_coeff_A)) + mat.repmat(data_mean, n, 1)  # be consistent with the text

    return df_pca_data, new_data_Mat, PCA_coeff_A, data_mean
