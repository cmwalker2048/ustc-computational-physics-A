import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm

# 画出 Markov chain 的所有点，并用不同颜色给它们上色
def markov_2d_color(filepath):
    x = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    
    n = len(x)
    
    point_numbers=list(range(len(x[0:n:1])))
    fig = plt.scatter(x[0:n:1],y[0:n:1],c=point_numbers,cmap=plt.cm.summer,edgecolors='none',s=1)
    plt.scatter(x[0],y[0],c='purple',s=50)
    plt.scatter(x[-1],y[-1],c='red',s=50)
    cbar = plt.colorbar(fig)
    cbar.set_label('$N \\ \\ $(Markov Chain Step)')
    plt.xlabel("$X$")
    plt.ylabel("$Y$")
    plt.title("Markov Chain on $X$-$Y$ plain ($\\beta = 5,\Delta r = 1$)")
    plt.show()
    pass

# 哈密顿量的表达式
def hamilton(x,y):
    return -2*(x**2+y**2) + 0.5*(x**4+y**4) + 0.5*(x-y)**4

# 3d空间中的 Markov Chain 的绘制
def markov_3d(filepath):
    
    x_p = np.linspace(-10,10,1000)
    y_p = np.linspace(-10,10,1000)
    
    x_p,y_p = np.meshgrid(x_p,y_p)
    
    H = hamilton(x_p,y_p)
    #'''
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(x_p,y_p,H,cmap=cm.Blues,label='$H(x,y)$',alpha = 0.6)
    c = fig.colorbar(surf, shrink=0.5, aspect=5)
    c.set_label("$H$")
    plt.xlabel("$X$")
    plt.ylabel("$Y$")
    ax.set_zlabel("$H$")
    plt.title("$H(x,y) \\ \\  x,y \in \left[-10,10\\right]^2$")
    
    x = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(0),encoding='utf-8')
    y = np.loadtxt(filepath,dtype=np.float32,delimiter=",",usecols=(1),encoding='utf-8')
    
    n = len(x)
    h = hamilton(x,y)
    
    point_numbers=list(range(len(x[0:n:1])))
    fig = ax.scatter(x[0:n:1],y[0:n:1],h[0:n:1],c=point_numbers,cmap=plt.cm.plasma,edgecolors='none',s=1)
    #plt.scatter(x[0],y[0],c='purple',s=50)
    #plt.scatter(x[-1],y[-1],c='red',s=50)
    cbar = plt.colorbar(fig)
    cbar.set_label('$N \\ \\ $(Markov Chain Step)')
    plt.xlabel("$X$")
    plt.ylabel("$Y$")
    plt.title("Markov Chain on $X$-$Y$ plain ($\\beta = 5,\Delta r = 1$)")
    plt.show()
    pass
    

if __name__ == "__main__":
    markov_2d_color("./m_rosenbluth_1.csv")
    markov_3d("./m_rosenbluth_1.csv")
