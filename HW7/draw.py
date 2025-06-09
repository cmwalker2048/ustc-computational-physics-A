import numpy as np #数据处理的库
import matplotlib.pyplot as plt #作图的库
import csv # CSV文件的读写库

def plot_dots(filepath,x_lable,y_lable,_title):
    # 将CSV文件中的随机数导入数组orign中
    orign1 = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    orign2 = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    # 作为x轴坐标的随机数点
    x = orign1
    # 作为y轴坐标的随机数点
    y = orign2

    plt.scatter(x,y,marker='.',c='red')
    plt.plot([2900,2994],[5700-0.001,5700+0.001],'--',c='red')
    plt.plot([2994,3004],[38000-0.001,38000+0.001],'--',c='red')
    plt.plot([3004,3013],[5700-0.001,5700+0.001],'--',c='red')
    plt.plot([2994-0.001,2994+0.001],[5700,38000],'--',c='red')
    plt.plot([3004-0.001,3004+0.001],[5700,38000],'--',c='red')
    plt.plot(x,y,'-')
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.title(_title)
    plt.show()
    #plt.savefig(_title+".png")
    
plot_dots('./data.csv','Energy($eV$)','Intensity(counts)','Data')