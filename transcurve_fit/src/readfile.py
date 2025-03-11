from scipy.io import loadmat


#!---------------  读取MATLAB的.mat格式文件  ---------------!#
def read_matlab(mat_name, X_key, Y_key ):
    """
    X_key, Y_key 分别为需要的数据列的名称, 返回是一个字典, 键值与输入的键值同名, 对应数据格式为nparray
    """
    matdata = loadmat(mat_name)
    
    if (X_key not in matdata.keys()) or (Y_key not in matdata.keys()) :
        raise KeyError(f"Key '{X_key}' or '{Y_key}' not found in .mat file. Available keys: {list(matdata.keys())}")
    else :
        dataset = {
            X_key : matdata[X_key][0],
            Y_key : matdata[Y_key][0]
        }

    return dataset




if __name__ == "__main__":

    mat_name = "9V_201926_ScanB_I7-1Vxx5-3&8-10Vxy10-3&8-5ground2-11_hengya100-10k_RI1k__T2.0011K_Vbg9V_Vbias0V_ac1V.mat"
    #data1 = read_matlab(mat_name, 'B', "Rxx2")
    #print(data1['B'], data1['Rxx2'])
    print("-----------------------")
    read_matlab(mat_name, 'B', "rxx2")

    
