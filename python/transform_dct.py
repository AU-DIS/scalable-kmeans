##DCT and iDCT transformation codes

import torch
import numpy as np
from itertools import product
import math




def dct(x, norm=None):

    x_shape = x.shape
    N = x_shape[-1]
    x = x.contiguous().view(-1, N)

    v = torch.cat([x[:, ::2], x[:, 1::2].flip([1])], dim=1)

    # Vc = torch.rfft(v, 1, onesided=False)           # comment this line
    Vc = torch.view_as_real(torch.fft.fft(v, dim=1))  # add this line

    k = - torch.arange(N, dtype=x.dtype, device=x.device)[None, :] * np.pi / (2 * N)
    W_r = torch.cos(k)
    W_i = torch.sin(k)

    V = Vc[:, :, 0] * W_r - Vc[:, :, 1] * W_i

    if norm == 'ortho':
        V[:, 0] /= np.sqrt(N) * 2
        V[:, 1:] /= np.sqrt(N / 2) * 2

    V = 2 * V.view(*x_shape)

    return V


def idct(X, norm=None):

    x_shape = X.shape
    N = x_shape[-1]

    X_v = X.contiguous().view(-1, x_shape[-1]) / 2

    if norm == 'ortho':
        X_v[:, 0] *= np.sqrt(N) * 2
        X_v[:, 1:] *= np.sqrt(N / 2) * 2

    k = torch.arange(x_shape[-1], dtype=X.dtype, device=X.device)[None, :] * np.pi / (2 * N)
    W_r = torch.cos(k)
    W_i = torch.sin(k)

    V_t_r = X_v
    V_t_i = torch.cat([X_v[:, :1] * 0, -X_v.flip([1])[:, :-1]], dim=1)

    V_r = V_t_r * W_r - V_t_i * W_i
    V_i = V_t_r * W_i + V_t_i * W_r

    V = torch.cat([V_r.unsqueeze(2), V_i.unsqueeze(2)], dim=2)

    
    # v = torch.irfft(V, 1, onesided=False)                             # comment this line
    v= torch.fft.irfft(torch.view_as_complex(V), n=V.shape[1], dim=1)   # add this line

    x = v.new_zeros(v.shape)
    x[:, ::2] += v[:, :N - (N // 2)]
    x[:, 1::2] += v.flip([1])[:, :N // 2]

    return x.view(*x_shape)



def dct_2d(x, norm=None):
    """
    2-dimentional Discrete Cosine Transform, Type II (a.k.a. the DCT)
    For the meaning of the parameter `norm`, see:
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html
    :param x: the input signal
    :param norm: the normalization, None or 'ortho'
    :return: the DCT-II of the signal over the last 2 dimensions
    """
    X1 = dct(x, norm=norm)
    X2 = dct(X1.transpose(-1, -2), norm=norm)
    return X2.transpose(-1, -2)


def idct_2d(X, norm=None):
    """
    The inverse to 2D DCT-II, which is a scaled Discrete Cosine Transform, Type III
    Our definition of idct is that idct_2d(dct_2d(x)) == x
    For the meaning of the parameter `norm`, see:
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html
    :param X: the input signal
    :param norm: the normalization, None or 'ortho'
    :return: the DCT-II of the signal over the last 2 dimensions
    """
    x1 = idct(X, norm=norm)
    x2 = idct(x1.transpose(-1, -2), norm=norm)
    return x2.transpose(-1, -2)


def dct_2d_flattened(x, norm=None):
    """
    2-dimentional Discrete Cosine Transform, Type II (a.k.a. the DCT)
    For the meaning of the parameter `norm`, see:
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html
    :param x: the input signal
    :param norm: the normalization, None or 'ortho'
    :return: the DCT-II of the signal over the last 2 dimensions
    """
    X1 = dct(x, norm=norm)
    X2 = dct(X1.transpose(-1, -2), norm=norm)
    return (X2.transpose(-1, -2)).flatten()


def dct_truncated(x, trunc):

    x_shape = x.shape
    N = x_shape[-1]
    x = x.contiguous().view(-1, N)

    v = torch.cat([x[:, ::2], x[:, 1::2].flip([1])], dim=1)

    # Vc = torch.rfft(v, 1, onesided=False)           # comment this line
    Vc = torch.view_as_real(torch.fft.fft(v, dim=1))  # add this line

    k = - torch.arange(N, dtype=x.dtype, device=x.device)[None, :] * np.pi / (2 * N)
    W_r = torch.cos(k)
    W_i = torch.sin(k)

    V = Vc[:, :, 0] * W_r - Vc[:, :, 1] * W_i

    V = 2 * V.view(*x_shape)

    return V[:trunc]


def dct_2d_truncated_flattened(x, trunc, norm=None):
    """
    2-dimentional Discrete Cosine Transform, Type II (a.k.a. the DCT)
    For the meaning of the parameter `norm`, see:
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html
    :param x: the input signal
    :param norm: the normalization, None or 'ortho'
    :return: the DCT-II of the signal over the last 2 dimensions
    """
    X1 = dct(x, norm=norm)
    X2 = dct(X1.transpose(-1, -2), norm=norm)
    return (X2.transpose(-1, -2))[:trunc, :trunc].flatten()


def make2d(array, cols=None, dtype=None):
    '''
    Make a 2D array from an array of arrays.  The `cols' and `dtype'
    arguments can be omitted if the array is not empty.
    '''
    if (cols is None or dtype is None) and not len(array):
        raise RuntimeError("cols and dtype must be specified for empty "
                           "array")

    if cols is None:
        cols = len(array[0])

    if dtype is None:
        dtype = array[0].dtype

    return np.fromiter(array, [('_', dtype, (cols,))],
                        count=len(array))['_'] 



##(Way1) calculate DCT features from raw dataset
def cal_dct_features(dct_len, dataset):
    coordinate_grid = list(product(range(dataset.shape[0]),range(dataset.shape[1])))
    #print(len(coordinate_grid))
    dct_image = list(
        map(dct_2d_truncated_flattened, iter(torch.tensor((dataset[i,j,:,:])) for i, j in coordinate_grid), [dct_len for _ in range(len(coordinate_grid))]))
    dct_image = make2d(dct_image, dtype=float)
    print(dct_image.shape)
    return dct_image



##(Way2) calculate DCT features from larger DCT features
##This requires source_dct_len is multiple of dest_dct_len, e.g. source_dct_len is 256 when dest_dct_len is 128
def trunc_dct_features(source_dct_image, dest_dct_len):
    dct_image = [[] for _ in range(len(source_dct_image))]

    for i in range(len(source_dct_image)):
        for j in range(dest_dct_len):
            for k in range(dest_dct_len):
                dct_image[i] += [source_dct_image[i][j*len(source_dct_image) + k]]
    dct_image = np.array(dct_image)
    print(dct_image.shape)
    return dct_image


## Make the squared sum pre-proc arrays
## copied from higher cell, remember to check consistency
def get_square_squared_sums(data):
    d_sqrt = int(math.sqrt(len(data[0])))
    log_d_sqrt = int(math.log(d_sqrt, 2))
    
    res = [[0 for _ in range(log_d_sqrt + 1)] for _ in range(len(data))]
    
    for i in range(len(data)):
        res[i][0] = (data[i][0] ** 2)
        for level in range(1, log_d_sqrt+1):
            
            if level > 0: res[i][level] += res[i][level-1]
            two_p_level_m1 = int(2**(level-1))
            two_p_level = int(2**level)
            for l in range(two_p_level_m1):
                res[i][level] += sum([folan**2 for folan in data[i][(d_sqrt*(two_p_level_m1 + l)): (d_sqrt*(two_p_level_m1 + l)) + two_p_level]])
                res[i][level] += sum([folan**2 for folan in data[i][l*d_sqrt + two_p_level_m1: l*d_sqrt + two_p_level]])
    res = np.array(res)
    print(res.shape)
    return res



def main():
    #Running DCT for c++ call
    filename = sys.argv[1]
    f = h5py.File(filename, 'r')
    dataset = f['data']

    dct_len = 1024
    image_dct = cal_dct_features(dct_len, dataset)
    #image_dct = trunc_dct_features(_image_dct,512)

    
    arr = np.array(image_dct)
    df = pd.DataFrame(data=arr)
    #min_ = df.min().min()
    #max_ = df.max().max()
    #mean_ = df.mean().mean()
    #df_norm = (df - min_) / (max_ - min_)

    df.to_csv(filename[:-3].replace("1024",str(dct_len))+"_h5_dct.txt", header=False, index=False, sep=" ", line_terminator=" \n")


if __name__ == "__main__":
    import sys
    import h5py
    import pandas as pd
    main()
