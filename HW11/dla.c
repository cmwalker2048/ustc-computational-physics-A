#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define PI 3.1415926 //圆周率
#define L 10000     //正方形边界的长度
#define K 3        //种子边界的长度
#define HALF_L L/2  //正方形边界长度的一半
#define HALF_K K/2  //种子边界长度的一半
//#define ALPHA       //逃逸半径系数，已弃用，用主函数中的参数

#define SAND_R K/2   //sandbox计数法 起始长度
#define STEP K/2     //sandbox计数法 间隔步长
#define SAND_N HALF_L/STEP  //sandbox计数法 采样点总步长数
#define THRUSTER 0.01   //sandbox阈值判断
#define TOL 50000   //总模拟点数

#define FILE_DATA "./dla_3_10_50000.csv"
#define FILE_SAND "./dla_sand_3_10_50000.csv"

//随机数产生器
float rn(){
    // 16807产生器的参数
    int A = 16807; //参数a
    int B = 0; //参数b
    long M = 2147483647; //参数m
    static long i = 23942907; //起始值i_0的大小(种子值)
    // Schrage法的参数
    long Q = 127773; //参数q（除数）
    int R = 2836; //参数r（余数）
    float x;
    while(1){
        x = (i*1.0) / M; //得到随机数
        i = A*(i % Q) - R*(i / Q);
        if(i < 0){
            i = i + M; //Schrage法的一个步骤
        }
        break;
    }
    return x;//返回随机数
}

//用于检查格点有无被占有的函数，若被占有，则返回1.
int check(int **dot,int x,int y){
    if(dot[x][y]==1){ //dot[x][y]代表第x行j列的占有状态
        return 1;
    }
    else{
        return 0;
    }
}

//自定义乘幂函数
double pow_2(long cor,double index){
    double temp = (double) cor;
    return pow(temp,index);//调用了系统乘幂函数
}

//sandbox计数法
    //传入占有状态矩阵，和两个储存数组
void sandbox(int **dot,int *sandbox_r,int *sandbox_n){
    int i,j;        //计数器
    int r = SAND_R; //sandbox计数法 起始长度
    int step = STEP;    //sandbox计数法 间隔步长
    int len = 0;    //N截至长度，记录N不再变大时的步长
    long counter=0; //计数器
    double delta=0; //误差记录器
    int k=0;        //计数器
    while(r<HALF_L){
        counter=0;
        //记录于r矩形内的点数
        for(i=HALF_L-r;i<HALF_L+r;i++){
            for(j=HALF_L-r;j<HALF_L+r;j++){
                if(dot[i][j]==1){
                    counter++;
                }
            }
        }
        sandbox_r[k]=r;         //储存数据
        sandbox_n[k]=counter;   //储存数据
        r+=step;                //增加r
        //判断结束条件
        if(k!=0){
            delta = sandbox_n[k]-sandbox_n[k-1];
            //若差值变化小的时候，则退出循环
            if(delta/sandbox_n[k]<THRUSTER){
                len = k;
                break;
            }
        }
        k++;
    }

    //写入文件
    FILE *fp =NULL;
    fp = fopen(FILE_SAND,"w");

    for(i=0;i<len;i++){
        fprintf(fp,"%d,%d\n",sandbox_r[i],sandbox_n[i]);
    }
    fclose(fp);
}

//主函数
void main(){
    long R;     //随机数粒子产生的位置半径
    float r;    //暂时储存随机数
    int alpha = 5;  //每次增加的小半径
    int ALPHA = 10;  //逃逸半径系数
    int thruster = alpha*10;
    long row = L;   //行数
    long col = L;   //列数
    long i;         //计数器
    long j;         //计数器
    int counter;    //计数器
    long x;         //寄存器,用于随机游走时储存x坐标
    long y;         //寄存器,用于随机游走时储存y坐标
    long cor[2];    //寄存器,用于随机游走时储存前一时刻的x,y坐标
    int **dot;      //dot[x][y]代表第x行y列的占有状态
    int *sandbox_r; //sandbox_r[x]用于记录R的值
    int *sandbox_n; //sandbox_n[x]用于记录N的值
    int flag=0;     //判断器
    sandbox_r = (int *)malloc(SAND_N*sizeof(int));  //给数组赋值长度
    sandbox_n = (int *)malloc(SAND_N*sizeof(int));  //给数组赋值长度
    dot = (int **)malloc(row*(sizeof(int*)));   //给数组赋值长度
    for(i=0;i<row;i++){
        dot[i] = (int *)malloc(col*(sizeof(int)));  //给数组赋值长度
    }

    //initialize
    printf("initializing...\n");
    for(i=HALF_L-HALF_K;i<=HALF_L+HALF_K;i++){
        for(j=HALF_L-HALF_K;j<=HALF_L+HALF_K;j++){
            dot[i][j] = 1;  //将种子值标记为1
        }
    }

    printf("begin\n");
    R = (long) sqrt(2)*HALF_K+alpha; //给出随机数初始半径
    while(counter<TOL){ //总点数约束
        flag =0; //初始化判读标志
        x = (long) (HALF_L+R*cos(2*PI*rn()));   //随机生成随机数x坐标
        y = (long) (HALF_L+R*sin(2*PI*rn()));   //随机生成随机数y坐标
        cor[0] = x;     //记录前一时刻x坐标
        cor[1] = y;     //记录前一时刻y坐标

        //当粒子没有离开逃逸半径时
        while(pow_2(cor[0]-HALF_L,2)+pow_2(cor[1]-HALF_L,2)<pow_2(ALPHA*R,2)){
            //若粒子的下一步走进种子，则将前一步记作占有状态，计数，并退出循环
            if(check(dot,x,y)){
                dot[cor[0]][cor[1]] = 1;
                counter++;
                flag =1;
                printf("order=%dth\n",counter);//指出占有状态数
                break;
            }
            else{           //若没有走进种子中，则继续随机游走
                cor[0]=x;
                cor[1]=y;
                r=rn();
                if(r<0.25){
                    x+=1;
                }
                else if(0.25<r && r<0.5){
                    x-=1;
                }
                else if(0.5<r && r<0.75){
                    y+=1;
                }
                else{
                    y-=1;
                }
            }
        }
        //若粒子撞进了种子，并且随机数半径以不够大时，更新半径
        if(flag==1 && pow_2(cor[0]-HALF_L,2)+pow_2(cor[1]-HALF_L,2)>pow_2(R,2)){
            R = R + alpha;  //更新随机数初始半径
        }
    }

    //将原始数据写入文件
    FILE *fp_data =NULL;
    fp_data = fopen(FILE_DATA,"w");

    for(i=0;i<L;i++){
        for(j=0;j<L;j++){
            if(dot[i][j]==1){
                fprintf(fp_data,"%d,%d\n",i,j);
            }
        }
    }
    fclose(fp_data);

    //进行Sandbox计数
    printf("calculating...\n");
    sandbox(dot,sandbox_r,sandbox_n);
}