import numpy as np # 数据处理
import matplotlib.pyplot as plt # 作图
from mpl_toolkits.mplot3d import Axes3D # 3D作图工具
import csv # 读取数据的工具

# 直接抽样主函数

    # 参数为:随机数文件路径,写入文件路径,生成的随机点数的2倍
def xy(filepath,outpath,n):
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    
    # 生成的theta,varphi随机数
    theta = np.arccos(2*orign[0:n:2]-1)
    varphi = 2*np.pi*orign[1:n+1:2]
    
    # 将其转换为x,y的随机数
    x = np.sin(theta)*np.cos(varphi)
    y = np.sin(theta)*np.sin(varphi)
    
    # 写入文件循环
    with open(outpath,'w') as out: # 参数为:写入文件路径,写入模式
        for j in range(0,int((n-1)/2)):
            out.write("%.4f,%.4f\n"%(x[j],y[j])) #写入文件
    
    # 以下为画图函数,不需要阅读!
    plt.xlabel('$X$')
    plt.ylabel('$Y$')
    plt.title('Direct sampling of 2500 points projection on $XY$ plane ')
    plt.scatter(x,y,marker='.')
    plt.show()


# Marsaglia抽样法主函数

    # 参数为:随机数文件路径,写入文件路径,生成的随机点数数目
def xy_marsaglia(filepath,outpath,n):
    # 将CSV文件中的随机数导入数组orign中
    orign = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    # 变量的转换
    transform = orign*2-1
    
    # 用于计点数的寄存器
    counter = 0
    # 用于遍历随机数的数
    j = 0
    
    x = [] # 符合要求的x坐标
    y = [] # 符合要求的y坐标

    # Marsaglia主循环函数
    while(counter <= n):
        u = transform[j]
        v = transform[j+1]
        z = u**2+v**2
        if(z <= 1):
            x.append(2*u*np.sqrt(1-z))
            y.append(2*v*np.sqrt(1-z))
            counter+=1
        j+=2
    
    #写入文件循环
    with open(outpath,'w') as out:
        for j in range(0,counter-1):
            out.write("%.4f,%.4f\n"%(x[j],y[j]))
    
    #以下为作图函数,不需要阅读!
    x_p = np.array(x)
    y_p = np.array(y)
    plt.xlabel('$X$')
    plt.ylabel('$Y$')
    plt.title('Marsaglia sampling of 2500 points projection on $XY$ plane ')
    plt.scatter(x_p,y_p,marker='.')
    plt.show()
    
# 利用直接抽样法
xy('./16807.csv','./out1.csv',5000)

# 利用Marsaglia抽样法
xy_marsaglia('./16807.csv','./marsaglia.csv',2500)