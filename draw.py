import matplotlib.pyplot as plt
import numpy as np
from math import *

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_zlim(-10,10)
plt.xlim(-10,10)
plt.ylim(-10,10)
ax.set_box_aspect([1,1,1])
plt.axis('off')

def reset():
    global ax
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_zlim(-10,10)
    plt.xlim(-10,10)
    plt.ylim(-10,10)
    ax.set_box_aspect([1,1,1])
    plt.axis('off')

def change_lim(x):
    ax.set_zlim(-x,x)
    plt.xlim(-x,x)
    plt.ylim(-x,x)

def get_x(p1,p2):
    return [p1[0],p2[0]]
def get_y(p1,p2):
    return [p1[1],p2[1]]
def get_z(p1,p2):
    return [p1[2],p2[2]]
def draw_line(A,B):
    ax.plot(get_x(A,B),get_y(A,B),get_z(A,B),color='#000000')
def draw(O,ver_list):
    for temp2 in ver_list:
        draw_line(O,temp2)
def end():
    global ax
    plt.show()
    plt.close()
    reset()

def R(u,cos,sin):
    ux,uy,uz = list(u)
    return np.array([[cos+ux**2*(1-cos),ux*uy*(1-cos) - uz*sin,ux*uz*(1-cos)+uy*sin],
                     [ux*uy*(1-cos)+uz*sin, cos+uy**2 * (1-cos), uy*uz*(1-cos) - ux*sin],
                     [uz*ux*(1-cos)-uy*sin, uz*uy*(1-cos)+ux*sin, cos+uz**2*(1-cos)]])
def rotation(u,cos,sin,A):
    return R(u,cos,sin)@A

def angle(ceta):
    return ceta * pi/180
def vec(A,O):
    return np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])

#line, tri, tetah, tri_pyra, bond
def find_pos_tetah(A,O):
    #평행이동, 회전행렬 생성
    AO = np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])
    lenght = np.linalg.norm(AO)
    A_O = np.array([lenght,0,0])
    if abs(AO[0])==abs(A_O[0]) and AO[1]==0 and AO[2]==0:
        u=np.array([0,0,1])
    else:
        u = np.cross(AO,A_O)
        u /= np.linalg.norm(u)

    cos1 = AO@np.array([[A_O[i]] for i in range(3)]) / lenght**2
    sin1 = sqrt(1-cos1[0]**2)

    #해당 도형 작도
    B= R(np.array([0,1,0]),1/3,sqrt(8)/3)@A_O + A_O
    b = np.array([B[0],0,0])
    bB = vec(b,B)
    Ob = vec(np.array([0,0,0]),b)
    C = R(np.array([1,0,0]),-1/2,sqrt(3)/2)@bB + Ob
    D = R(np.array([1,0,0]),-1/2,-sqrt(3)/2)@bB + Ob

    #평행 및 회전 이동
    B = (R(u,cos1[0],-sin1)@B) +A
    C = (R(u,cos1[0],-sin1)@C) +A
    D = (R(u,cos1[0],-sin1)@D) +A
    return [B,C,D]

def find_pos_line(A,O):
    #평행이동, 회전행렬 생성
    AO = np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])
    lenght = np.linalg.norm(AO)
    A_O = np.array([lenght,0,0])
    if abs(AO[0])==abs(A_O[0]) and AO[1]==0 and AO[2]==0:
        u=np.array([0,0,1])
    else:
        u = np.cross(AO,A_O)
        u /= np.linalg.norm(u)
    cos1 = AO@np.array([[A_O[i]] for i in range(3)]) / lenght**2
    sin1 = sqrt(1-cos1[0]**2)

    #해당 도형 작도
    B= A_O + A_O
    #평행 및 회전 이동
    B = (R(u,cos1[0],-sin1)@B) +A
    return [B]

def find_pos_tri(A,O):
    #평행이동, 회전행렬 생성
    AO = np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])
    lenght = np.linalg.norm(AO)
    A_O = np.array([lenght,0,0])
    if abs(AO[0])==abs(A_O[0]) and AO[1]==0 and AO[2]==0:
        u=np.array([0,0,1])
    else:
        u = np.cross(AO,A_O)
        u /= np.linalg.norm(u)
    cos1 = AO@np.array([[A_O[i]] for i in range(3)]) / lenght**2
    sin1 = sqrt(1-cos1[0]**2)

    #해당 도형 작도
    B= R(np.array([0,1,0]),1/2,sqrt(3)/2)@A_O + A_O
    C = R(np.array([0,1,0]),1/2,-sqrt(3)/2)@A_O + A_O

    #평행 및 회전 이동
    B = (R(u,cos1[0],-sin1)@B) + A
    C = (R(u,cos1[0],-sin1)@C) + A

    return [B,C]

def find_pos_bond(A,O):
    #평행이동, 회전행렬 생성
    AO = np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])
    lenght = np.linalg.norm(AO)
    A_O = np.array([lenght,0,0])
    if abs(AO[0])==abs(A_O[0]) and AO[1]==0 and AO[2]==0:
        u=np.array([0,0,1])
    else:
        u = np.cross(AO,A_O)
        u /= np.linalg.norm(u)
    cos1 = AO@np.array([[A_O[i]] for i in range(3)]) / lenght**2
    sin1 = sqrt(1-cos1[0]**2)

    #해당 도형 작도
    B= R(np.array([0,1,0]),1/4,sqrt(15)/4)@A_O + A_O

    #평행 및 회전 이동
    B = (R(u,cos1[0],-sin1)@B) +A
    return [B]

def find_pos_tri_pyra(A,O):
    #평행이동, 회전행렬 생성
    AO = np.array([O[0]-A[0],O[1]-A[1],O[2]-A[2]])
    lenght = np.linalg.norm(AO)
    A_O = np.array([lenght,0,0])
    if abs(AO[0])==abs(A_O[0]) and AO[1]==0 and AO[2]==0:
        u=np.array([0,0,1])
    else:
        u = np.cross(AO,A_O)
        u /= np.linalg.norm(u)
    cos1 = AO@np.array([[A_O[i]] for i in range(3)]) / lenght**2
    sin1 = sqrt(1-cos1[0]**2)

    #해당 도형 작도
    B= R(np.array([0,1,0]),cos(angle(68.83)),sin(angle(68.83)))@A_O + A_O
    b = np.array([B[0],0,0])
    bB = vec(b,B)
    Ob = vec(np.array([0,0,0]),b)
    C = R(np.array([1,0,0]),-1/2,sqrt(3)/2)@bB + Ob
    D = R(np.array([1,0,0]),-1/2,-sqrt(3)/2)@bB + Ob

    #평행 및 회전 이동
    B = (R(u,cos1[0],-sin1)@B) +A
    C = (R(u,cos1[0],-sin1)@C) +A
    D = (R(u,cos1[0],-sin1)@D) +A
    return [B,C,D]

def change_distance(A,O,lenght):
    AO = np.array([A])-np.array([O])
    ori_lenght = np.linalg.norm(AO)
    return (((lenght/ori_lenght)*(AO))+np.array([O]))[0]

def sphere(r,O,color):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = r * np.outer(np.cos(u), np.sin(v)) +O[0]
    y = r * np.outer(np.sin(u), np.sin(v)) +O[1]
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) +O[2]
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color=color, linewidth=0, alpha=1)

#길이 늘리기 = 2*(B-O)+O

#필요한 것 : 원소 종류 -> 길이
#계산할 것 : 길이를 기반으로 좌표 계산
#각 원자가 들고있을 정보 : 좌표, 길이