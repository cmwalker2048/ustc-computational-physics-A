# ustc-computational-physics-A
This is a repository to record my HW &amp; Project for Computational Physics A

## HW1
用 Schrage 方法编写随机数子程序，用指定间隔（非连续 $l >1$）两个随机数作为点的坐标值绘出若干点的平面分布图。再用 $<x^k>$ 测试均匀性（取不同量级的 $N$ 值，讨论偏差与 $N$ 的关系）、$C(l)$测试其 2 维独立性（总点数 $N > 10^7$）。

## HW2 
用 16807 产生器测试随机数序列中满足关系 $X_{n−1} > X_{n+1} > X_n$ 的比重。讨论 Fibonacci 延迟产生器中出现这种关系的比重。

## HW3
在球坐标系 $(\rho, \theta, \phi)$ 下，产生上半球面上均匀分布的随机坐标点，给出其直接抽样方法。

## HW4 
设 pdf 函数满足关系式

$$ p'(x) = a \delta(x) + b \exp (-cx) , x \in \left[ -1,1 \right] , a \neq 0 $$

讨论该函数性质并给出抽样方法。

## HW5
对于球面上均匀分布的随机坐标点，给出它们在 $(x,y)$ 平面上投影的几率分布函数。并由此验证Marsaglia抽样方法 $x=2u \sqrt{1-r^2},y=2v\sqrt{1-r^2},z=1-2r^2$ 确为球面上均匀分布的随机抽样。

## HW6
对两个函数线型（Gauss分布和类Lorentz 型分布），设其一为  $p(x)$ ，另一为
 $F(x)$ ，其中常数 $a \neq b \neq 1$ , 用舍选法对 $p(x)$ 抽样。将计算得到的归一化频数分布直方图
与理论曲线 $p(x)$ 进行比较，讨论差异，讨论抽样效率。

$$
    Gaussian: \sim  exp{-ax^2}; \quad Lorentzian \ like: \sim \frac{1}{1+bx^4}
$$

## HW7
对一个实验谱数值曲线 $p(x)$ ，自设 $F(x)$ ，分别用直接抽样和舍选法
对 $p(x)$ 抽样。比较原曲线和抽样得到的曲线以验证。讨论抽样效率。

## HW8
用Monte Carlo方法计算如下定积分，并讨论有效数字位数。

$$
    I_1 = \int_{0}^{5} \,dx \sqrt{x^2+2\sqrt{x}}
$$

$$
    I_2 = \int_{0}^{\frac{7}{10}} \int_{0}^{\frac{4}{7}} \int_{0}^{\frac{9}{10}} \int_{0}^{2} \int_{0}^{\frac{13}{11}}  \,dx  \,dy  \,dz  \,du  \,dv
    (5+x^2-y^2+3xy-z^2+u^3-v^3) 
$$

## HW9
考虑泊松分布、指数分布，并再自设若干个随机分布（它们有相同或不同
的 $\mu$ 和 $\sigma^2$ ），通过Monte Carlo模拟，验证中心极限定理成立（N =2、5、10）。

## HW10
Monte Carlo方法研究二维平面上荷电粒子在正弦外电场（ $\sim \sin{\omega t}$ ）中的
随机行走。推导速度自相关函数的表达式，它随时间的变化是怎样的行为？能否模拟
得到该自相关函数的曲线？是的话与理论曲线进行比较，否的话讨论理由。

## HW11
模拟 2 维 DLA 以及介电击穿（DBM）图案并讨论

## HW12
推导正方格子点阵上键逾渗的重整化群变换表达式 $p' = R(p)$ ，求临界点 $p_c$ 与临
界指数
 $\nu$ ，与正确值（表1.6.1.3-1）相比较。

| 维数 | 点阵  | 座逾渗 $p_c$ | 键逾渗 $p_c$ | 配位数 |
|----|-----|----------|----------|-----|
| 2  | 正方形 | 0.592746 | 0.50000  | 4   |

## HW13
用Metropolis-Hasting抽样方法计算积分：

$$
    I = \int_{0}^{\infty} (x-\alpha \beta)^2 f(x) \,dx = \alpha \beta^2
$$

$$
    f(x) = \frac{1}{\beta \Gamma(\alpha)} (\frac{x}{\beta})^{\alpha-1} \exp{-\frac{x}{\beta}}
$$

设积分的权重函数为： $p(x)=f(x)$ 和 $p(x)=(x-\alpha\beta)^2 f(x)$

给定参数 $\alpha,\beta$ ，并用不同的 $\gamma$ 值，分别计算积分，讨论计算精度和效率

其中，设 $T_{ij} = T(x \rightarrow x') = T(x') = 0.5 \exp{-\frac{x'}{\gamma}}$

## HW14
**苏格拉底：** 诘问法是发现真理和明确概念的有效方法，请同学们以Ising经典自旋模型为例，
论述相空间、Liouville定理、正则系综、Markov链等概念。

**学生A：** 相空间是以 $N$ 个粒子的位置坐标 $q$ 和动量 $p$ 展开的 $6N$ 维空间。Ising模型中的
Hamiltonian仅与自旋变量有关，与坐标和动量无关， $\frac{\partial H}{\partial q} = \frac{\partial H}{\partial p} = 0$ ，因此： $\left[ \rho , H \right] = 0$  ,
即Liouville定理成立,  $\frac{d \rho}{d t} = \left[ \rho,H \right] = 0$ ，几率密度分布因此也为$H$的函数，因此它就是正则
系综中的Boltzmann分布： $\rho \varpropto \exp{-\beta H}$

**学生B：** 非也。将自旋作为广义坐标，则同样得到自旋也是广义动量。相空间是以物理问
题中的自由度为坐标展开的高维空间，对 $N$ 个自旋体系展开的则是 $N$ 维空间，空间的每
一维坐标只有两个取值： $+1$ 和 $-1$ 。如对 $2$ 个自旋的相空间，代表点只能取 $\left(+1,+1\right)$ 、 $\left(+1,-1\right)$ 、
 $\left(-1,+1\right)$ 、 $\left(-1,-1\right)$ 这 $4$ 个点。类似地，多自旋情况下代表点也只能位于多维相空间立方盒子
的顶点上。不同于坐标 $q$ 和动量 $p$ 组成的相空间中代表点是流动的情况，现在这些代表点
是与时间无关的，即密度不随时间改变的，因此 $\frac{d \rho}{d t}=0$ 。

**学生A：** 我不能同意你的观点。如果相空间是这样的话，由于代表点只能取在顶点上，连
几率密度分布本身都是离散的，而不是在该相空间中连续分布的。另外，
 $\frac{d \rho}{d t} = \sum_{i} (\frac{d \rho}{d \sigma_i})(\frac{d \sigma_i}{d t})$ ，在无穷小的时间变化 $dt$ 内，自旋的变化 $\Delta \sigma$ 则是有限的，不
能得到Liouville定理。更何况系综理论推导时基于的也是 $\left(q , p\right)$ 变量。

**学生C：**（请以学生C的身份参与辩论）

## HW15
设体系的能量为  $H(x,y)=-2(x^2+y^2)+\frac{1}{2}(x^4+y^4)+\frac{1}{2}(x-y)^4$ ，取 $\beta = 0.2,1,5$ ，
采用Metropolis抽样法计算 $\left\langle x^2 \right\rangle$ , $\left\langle y^2 \right\rangle$ , $\left\langle x^2+y^2 \right\rangle$ 。抽样时在2维平面上依次标出Markov链点
分布，从而形象地理解Markov链。

## HW16
以 $x_{n+1} = \lambda \sin{ \pi x_n}$ 为迭代方程进行迭代：

1. 画出系统状态随参数 $\lambda$ 的变化图，要求在图中体现出定值状态、
倍周期分叉和混沌状态；
2. 列出各个倍周期分叉处的 $\lambda$ 值，求相应的 $Feigenbaum$ 常数。

## HW17
进行单中心 DLA 模型的模拟 (可以用圆形边界，也可以用正方形边界)，并用两种方法计算模
拟得到的 DLA 图形的分形维数，求分形维数时需要作出双对数图。
**注：由于在第 11 次作业中已经对 DLA 模型做过一些讨论，在这里会参考一部分第 11 次作业
中的内容。**
