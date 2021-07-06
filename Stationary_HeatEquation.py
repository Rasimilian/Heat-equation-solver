import numpy as np
from numpy import sin, pi, cos,exp
import matplotlib.pyplot as plt
from cmath import phase
from mpl_toolkits.mplot3d import Axes3D
z=open('temperature.txt','w')
z1=open('temperature1.txt','w')
z2=open('temperature2.txt','w')



def func(x):                         # гран условие
    return 2+sin(3*x)+cos(x)
l=200
d_fi=2*pi/l                      # l-кол-во слагаемых в разложении в ряд Фурье
signal=[func(x*d_fi) for x in range(0,l)]
fourier_plus=np.fft.rfft(signal)*2/l
# fourier_minus=fourier_plus.conjugate()
# fourier_minus=np.delete(fourier_minus,0)
# fourier_minus=np.flip(fourier_minus)
# fourier=np.concatenate((fourier_minus,fourier_plus))
real_fourier=[i.real for i in fourier_plus]
imag_fourier=[i.imag for i in fourier_plus]
# phase_fourier=[phase(i) for i in fourier]

# x={0,...,R}
# t={0,...,inf}
# Задаем толщину диска через внутренний и внешний радиусы
R=1  
R1=2
u_x_0=0 # Гран условие на внутреннем радиусе
# u_x_0=[i for i in real_fourier]
# u_x_0_1=[i for i in imag_fourier]
# Гран условие на внешнем радиусе
u_x_n=[i for i in real_fourier]         # считать тут коэффициенты
u_x_n=[u_x_n[i]/R1**(i) for i in range(0,int(l/2)+1)] 
u_x_n1=[i for i in imag_fourier]         # считать тут коэффициенты
u_x_n1=[u_x_n1[i]/R1**(i) for i in range(0,int(l/2)+1)] 
# u_x_n=[0]

def angle(x):
    return complex(cos(x),sin(x))
n=100
solution=np.zeros((int(l/2)+1,n+1))
solution1=np.zeros((int(l/2)+1,n+1))



def u_t_0(x):
    return sin(x)
n=100
h=(R1-R)/n
t=0.0000001
setka=[i*h for i in range(0,n+1)]

m1=0


for m in range(0,int(l/2)+1):
    matrix_d=np.zeros(n-1)
    matrix_d[0]=0
    matrix_d[n-2]=-u_x_n[m1]*(1/h+(2*m+1)/(R+h*(n-1)))

    roots=[]
    # Tridiagonal matrix input
    matrix=np.zeros((n-1,n-1))
    for i in range(0,n-2):
        matrix[i][i+1]=1/h+(2*m+1)/(R+h*(i+1))                       # c-coefficients
        matrix[i+1][i]=1/h                                        # a-coefficients
    for i in range(0,n-1):
        matrix[i][i]=-2/h-(2*m+1)/(R+h*(i+1))                        # b-coefficients
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
    roots.append(u_x_n[m1])
    roots.reverse()
    roots.append(u_x_0)
    roots.reverse()

    
    
    matrix_d=np.zeros(n-1)
    matrix_d[0]=0
    matrix_d[n-2]=-u_x_n1[m1]*(1/h+(2*m+1)/(R+h*(n-1)))
    
    roots1=[]
    # Tridiagonal matrix input
    matrix=np.zeros((n-1,n-1))
    for i in range(0,n-2):
        matrix[i][i+1]=1/h+(2*m+1)/(R+h*(i+1))                       # c-coefficients
        matrix[i+1][i]=1/h                                        # a-coefficients
    for i in range(0,n-1):
        matrix[i][i]=-2/h-(2*m+1)/(R+h*(i+1))                        # b-coefficients
    for i in range(0,n-2):
        c=matrix[i+1][i]/matrix[i][i]
        matrix[i+1][i]=0
        matrix[i+1][i+1]=matrix[i+1][i+1]-c*matrix[i][i+1]
        matrix_d[i+1]=matrix_d[i+1]-c*matrix_d[i]
    x=matrix_d[n-2]/matrix[n-2][n-2]
    roots1=[x]
    for i in range(n-3,-1,-1):
        x=(matrix_d[i]-matrix[i][i+1]*x)/matrix[i][i]
        roots1.append(x)
    roots1.reverse()
    roots1.append(u_x_n1[m1])
    roots1.reverse()
    roots1.append(u_x_0)
    roots1.reverse()    
    
    # assembling of Um(r,t) without m_angular part
    for i in range(0,n+1):
        solution[m1][i]=(R+h*i)**m*roots[i]
        solution1[m1][i]=(R+h*i)**m*roots1[i]
    m1=m1+1

## angular part adding
## Для r=R=1 зависимость T(fi)
setka_fi=[]
delta_fi=2*pi/100
function=[]
for k in range(0,101):
    setka_fi.append(delta_fi*k)
    m1=0
    x=0
    for m in range(0,int(l/2)+1):
        x=x+(complex(solution[m][n],solution1[m][n])*angle(delta_fi*k*m1)).real

        m1=m1+1
    function.append(x)

plt.plot(setka_fi,function)
plt.show()

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

setka_fi=[]
delta_fi=2*pi/100
function=[]
x_axis=[]
y_axis=[]
for j in range(0, 101):
    setka_fi.append(delta_fi*j)
    for k in range(0,101):
        m1=0 #m1=-l
        x=0
        for m in range(0,int(l/2)+1):
            x=x+(complex(solution[m][k],solution1[m][k])*angle(delta_fi*j*m1)).real
            m1=m1+1
        function.append(x)
        x_axis.append((R+h*k)*cos(delta_fi*j))
        y_axis.append((R+h*k)*sin(delta_fi*j))
        
print(len(x_axis))
print(len(y_axis))
print(len(function))

l=[str(i) for i in function]
for i in l:
    z.write(i+'\n')
l=[str(i) for i in x_axis]
for i in l:
    z1.write(i+'\n')
    
l=[str(i) for i in y_axis]
for i in l:
    z2.write(i+'\n')


fig = plt.figure()
ax3D = fig.add_subplot(111, projection='3d')
ax3D.scatter(x_axis,y_axis,function, s=10, marker='o')
plt.show()

z.close()
z1.close()
z2.close()



