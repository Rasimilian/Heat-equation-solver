import numpy as np
from numpy import sin, pi, cos,exp
import matplotlib.pyplot as plt
from cmath import phase
from mpl_toolkits.mplot3d import Axes3D
import sys
z=open('temperature.txt','w')        ## Температура - функция
z1=open('temperature1.txt','w')      ## х - координата
z2=open('temperature2.txt','w')      ## у - координата
z3=open('temperature3.txt','w')      ## Температура - функция
z4=open('temperature4.txt','w')      ## х - координата
z5=open('temperature5.txt','w')      ## у - координата
z6=open('temperature6.txt','w')      ## t - время


def func(x):                         ## гран условие на границе диска
    return 1+sin(3*x)
l=200
d_fi=2*pi/l                      ## l-кол-во слагаемых в разложении в ряд Фурье (в 2 раза меньше на самом деле)
signal=[func(x*d_fi) for x in range(0,l)]        
fourier_plus=np.fft.rfft(signal)*2/l           ## Нормировка неверная (2/l)
real_fourier=[i.real for i in fourier_plus]
imag_fourier=[i.imag for i in fourier_plus]

# x={0,...,R}
# t={0,...,inf}
R=1
# R1=1
u_x_0=[0 for i in range(0,int(l/2)+1)]
u_x_0_1=[0 for i in range(0,int(l/2)+1)]
# u_x_0=[i for i in real_fourier]
# u_x_0=[u_x_0[i]/R**(i) for i in range(0,int(l/2)+1)]          
# u_x_0_1=[i for i in imag_fourier]
# u_x_0_1=[u_x_0_1[i]/R**(i) for i in range(0,int(l/2)+1)]

u_x_n=[i for i in real_fourier]         ## Коэффициенты фурье для гран условия
u_x_n=[u_x_n[i]/R**(i) for i in range(0,int(l/2)+1)]
u_x_n1=[i for i in imag_fourier]        ## Коэффициенты фурье для гран условия
u_x_n1=[u_x_n1[i]/R**(i) for i in range(0,int(l/2)+1)]
# u_x_n=[0]

## Функция собирает комплексное число из реал и мним частей
def angle(x):
    return complex(cos(x),sin(x))

n=200
## Матрицы без учета времени (стационарное решение)
solution=np.zeros((int(l/2)+1,n+1))
solution1=np.zeros((int(l/2)+1,n+1))

## Матрицы с учетом времени (Распределение температуры в каждый момент времени на двумерной сетке)
solution_time=np.zeros((int(l/2)+1,2001,n+1))
print(solution_time.shape)
solution_time1=np.zeros((int(l/2)+1,2001,n+1))

## Произвольное гран условие при t=0
def u_t_0(x):
    return 1+x

h=R/n    ## Шаг по сетке радиуса
t=0.001  ## Шаг по сетке времени
setka=[i*h for i in range(0,n+1)] 
m1=0

result=[]
for m in range(0,int(l/2)+1):
    matrix_d=[u_t_0(h*i)/t for i in range(1,n)]
    matrix_d[0]=u_t_0(h*1)/t+u_x_0[m]*(1/h**2)
    matrix_d[n-2]=u_t_0(h*(n-1))/t+u_x_n[m1]*(1/(h*(n-1))/h+1/h**2)
    
    matrix_d1=[u_t_0(h*i)/t for i in range(1,n)]
    matrix_d1[0]=u_t_0(h*1)/t+u_x_0[m]*(1/h**2)
    matrix_d1[n-2]=u_t_0(h*(n-1))/t+u_x_n1[m1]*(1/(h*(n-1))/h+1/h**2)

    time=[]
#     roots=[]
    roots1=np.zeros((n+1))
    real_roots_last=0
    real_roots=1
    chetchik=0
    time_dependence=[]

#     for w in range(0,601):
    while abs(abs(real_roots_last)-abs(real_roots))>=sys.float_info.epsilon:    
## Метод установления. Норма вектора решения в предпоследней итерации. Сумма модулей элементов (каждый узел сетки) комплексного числа roots[i]+j*roots1[i]
        real_roots_last=0
        for i in roots1:
            real_roots_last=real_roots_last+i

        # Tridiagonal matrix input
        matrix=np.zeros((n-1,n-1))
        for i in range(0,n-2):
            matrix[i][i+1]=-1/h**2-1/(h*(i+1))/h                       # c-coefficients
            matrix[i+1][i]=-1/h**2                                         # a-coefficients
        for i in range(0,n-1):
            matrix[i][i]=1/t+2/(h**2)+1/(h*(i+1))/h+(m**2)/((h*(i+1))**2)                        # b-coefficients
        for i in range(0,n-2):
            c=matrix[i+1][i]/matrix[i][i]
            matrix[i+1][i]=0
            matrix[i+1][i+1]=matrix[i+1][i+1]-c*matrix[i][i+1]
            matrix_d[i+1]=matrix_d[i+1]-c*matrix_d[i]
        x=matrix_d[n-2]/matrix[n-2][n-2]
        roots=[x]
        for i in range(n-3,-1,-1):
            x=(matrix_d[i]-matrix[i][i+1]*x)/matrix[i][i]
            roots.append(x)
        roots.reverse()

        matrix_d=[roots[i]/t for i in range(1,n-2)]
        matrix_d.append(roots[n-2]/t+u_x_n[m1]*(1/(h**2)+1/(h*(n-1))/(h)))
        matrix_d.reverse()
        matrix_d.append(roots[0]/t+u_x_0[m]/h**2)
        matrix_d.reverse()
        roots.append(u_x_n[m1])
        roots.reverse()
        roots.append(u_x_0[m])
        roots.reverse()
        
        

         # Tridiagonal matrix input
        matrix=np.zeros((n-1,n-1))
        for i in range(0,n-2):
            matrix[i][i+1]=-1/h**2-1/(h*(i+1))/h                       # c-coefficients
            matrix[i+1][i]=-1/h**2                                         # a-coefficients
        for i in range(0,n-1):
            matrix[i][i]=1/t+2/(h**2)+1/(h*(i+1))/h+(m**2)/((h*(i+1))**2)                        # b-coefficients
        for i in range(0,n-2):
            c=matrix[i+1][i]/matrix[i][i]
            matrix[i+1][i]=0
            matrix[i+1][i+1]=matrix[i+1][i+1]-c*matrix[i][i+1]
            matrix_d1[i+1]=matrix_d1[i+1]-c*matrix_d1[i]
        x=matrix_d1[n-2]/matrix[n-2][n-2]
        roots1=[x]
        for i in range(n-3,-1,-1):
            x=(matrix_d1[i]-matrix[i][i+1]*x)/matrix[i][i]
            roots1.append(x)
        roots1.reverse()
    
        matrix_d1=[roots1[i]/t for i in range(1,n-2)]
        matrix_d1.append(roots1[n-2]/t+u_x_n1[m1]*(1/(h**2)+1/(h*(n-1))/(h)))
        matrix_d1.reverse()
        matrix_d1.append(roots1[0]/t+u_x_0_1[m]/h**2)
        matrix_d1.reverse()
        roots1.append(u_x_n1[m1])
        roots1.reverse()
        roots1.append(u_x_0_1[m])
        roots1.reverse()
## Метод установления. Норма вектора решения в последней итерации. Сумма модулей элементов (каждый узел сетки) комплексного числа roots[i]+j*roots1[i]         
        real_roots=0
        for i in roots1:
            real_roots=real_roots+i

        chetchik=chetchik+1
#         print(chetchik)
#         time.append(chetchik)

        if chetchik<201:
            for i in range(0,n+1):
                solution_time[m][chetchik][i]=roots[i]
                solution_time1[m][chetchik][i]=roots1[i]
#                 solution_time[m][w][i]=roots[i]
#                 solution_time1[m][w][i]=roots1[i]
    # assembling of Um(r,t) without m_angular part    
    for i in range(0,n+1):
        solution[m1][i]=roots[i]
        solution1[m1][i]=roots1[i]
    m1=m1+1
    print(m1,chetchik)

    


## angular part adding
## Для r=R=1 зависимость T(fi)
# setka_fi=[]
# delta_fi=2*pi/100
# function=[]
# for k in range(0,101):
#     setka_fi.append(delta_fi*k)
#     m1=0 #m1=-l
#     x=0
#     for m in range(0,int(l/2)+1):
#         x=x+(complex(solution[m][n],solution1[m][n])*angle(delta_fi*k*m1)).real
# 
#         m1=m1+1
#     function.append(x)
# 
# plt.plot(setka_fi,function)
# plt.show()

## Для угла 2pi/3 зависимость T(r) 
# function=[]
# for k in range(0,101):
#     m1=-l
#     x=0
#     for m in range(0,2*l+1):
#         x=x+(complex(solution[m][k],solution1[m][k])*angle(2*pi/3*m1)).real
# 
#         m1=m1+1
#     function.append(x)
# plt.plot(setka,function)
# plt.show()

# setka_fi=[]
# ugol=20
# delta_fi=2*pi/ugol
# function=[]
# x_axis=[]
# y_axis=[]
# for j in range(0, ugol+1):
#     setka_fi.append(delta_fi*j)
#     for k in range(0,n+1):
#         m1=0
#         x=0
#         for m in range(0,int(l/2)+1):
#             x=x+(complex(solution[m][k],solution1[m][k])*angle(delta_fi*j*m1)).real
#             m1=m1+1
#         function.append(x)
#         x_axis.append((h*k)*cos(delta_fi*j))
#         y_axis.append((h*k)*sin(delta_fi*j))
        
# print(len(x_axis))
# print(len(y_axis))
# print(len(function))

# l=[str(i) for i in function]
# for i in l:
#     z.write(i+'\n')
# l=[str(i) for i in x_axis]
# for i in l:
#     z1.write(i+'\n')
#     
# l=[str(i) for i in y_axis]
# for i in l:
#     z2.write(i+'\n')


# fig = plt.figure()
# ax3D = fig.add_subplot(111, projection='3d')
# ax3D.scatter(x_axis,y_axis,function, s=10, marker='o')
# plt.show()

l=200
ugol=50
delta_fi=2*pi/ugol
function=np.zeros(((ugol+1)*(n+1),200))
time=[]

t_chetchik=0
# ## Цикл времени
for t in range(0,200):
    print(t)
    stroka=0
    x_axis=[]
    y_axis=[]
    ## Цикл угла
    for j in range(0,ugol+1):
        ## Цикл радиуса
        for k in range(0,n+1):
            m1=0
            x=0
            ## Цикл гармоник
            for m in range(0,int(l/2)+1):
                x=x+(complex(solution_time[m][t][k],solution_time1[m][t][k])*angle(delta_fi*m*j)).real
                m1=m1+1
            function[stroka][t_chetchik]=x
            x_axis.append((h*k)*cos(delta_fi*j))
            y_axis.append((h*k)*sin(delta_fi*j))
            stroka=stroka+1
#     absolute=np.zeros((len(function[:][t_chetchik])+1))
#     for i in range(0,len(function[:][t_chetchik])+1):
#         absolute[i]=abs(function[i][t_chetchik])
#     normirovka=max(absolute)
#     for i in range(0,len(function[:][t_chetchik])+1):
#         function[i][t_chetchik]=function[i][t_chetchik]/normirovka
#     print(max(function[:][t_chetchik]))
#     print(min(function[:][t_chetchik]))
    time.append(t)
    t_chetchik=t_chetchik+1


# l=[str(i) for i in function]
# for i in l:
#     z3.write(i+'\n')
np.savetxt(z3,function)

l=[str(i) for i in x_axis]
for i in l:
    z4.write(i+'\n')
    
l=[str(i) for i in y_axis]
for i in l:
    z5.write(i+'\n')
    
l=[str(i) for i in time]        
for i in l:
    z6.write(i+'\n')
# 
z.close()
z1.close()
z2.close()
z3.close()
z4.close()
z5.close()
z6.close()


