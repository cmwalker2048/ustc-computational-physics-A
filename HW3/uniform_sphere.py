import numpy as np # 数据处理
import matplotlib.pyplot as plt # 作图
from mpl_toolkits.mplot3d import Axes3D # 3D作图工具
import csv # 读取数据的工具


# 数据抽样与作图主函数

    # 参数为:随机数文件路径,写入文件路径,产生随机点数的两倍(2n)
def uniform_sphere(filepath,outpath,n):
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
   
    theta = np.arccos(orign[0:n:2]) # theta的样本点
    varphi = 2*np.pi*orign[1:n+1:2] # varphi的样本点
    
    # 写入文件的主函数
    with open(outpath,'w') as out: #参数为:写入文件路径,写入模式
        for j in range(0,int((n-1)/2)):
            # 写入数据
            out.write("%.4f,%.4f\n"%(theta[j],varphi[j]))
            
    # 作图时需要的x,y,z坐标
    x = np.sin(theta)*np.cos(varphi)
    y = np.sin(theta)*np.sin(varphi)
    z = np.cos(theta)  
    
    # 以下为作图的函数,不需要阅读!
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title("2500 points on the sphere")
    #ax1.get_proj = lambda: np.dot(Axes3D.get_proj(ax1), np.diag([1, 1, 1, 1]))
    ax1.scatter(x,y,z,marker='.')
    #ax1.view_init(elev=90,azim=-90)
    ax1.view_init(elev=0,azim=-90)
    plt.show()
    
# 运行以上的程序
uniform_sphere('./16807.csv','./out.csv',5000)