from operator import is_not
from functools import partial
import numpy as np
import copy
import itertools


DNA_list=[]

with open('Overlap_Graphs_data.txt','r') as data_Rosa:
    for line in data_Rosa:
        DNA_list.append(line.strip())

def join_string(L):
    a=0
    M=[]
    for id, elem in enumerate(L):
        if elem.isalnum() is not True:
            a=0
            M.append(elem)
        else:
            a+=1
            if a>1:
                M[len(M)-1]=M[len(M)-1]+elem
            else: M.append(elem)
    return M

def max_overlap(s1,s2):
    p=0
    ml=min(len(s1),len(s2))
    for k in range(ml):
        if s1[-ml+k::]==s2[:ml-k:]:
            p=ml-k
            break
    return p

DNA_list_join=join_string(DNA_list)

agj_matrix=[[0]*(len(DNA_list_join)//2) for i in range(len(DNA_list_join)//2)]

def gener_agj_matr(lst):
    for num1, elem1 in enumerate(DNA_list_join):
        if elem1.isalnum():
            for num2, elem2 in enumerate(DNA_list_join):
                if elem2.isalnum() and elem1!=elem2:
                    agj_matrix[num1//2][num2//2]=max_overlap(elem1,elem2)
    return agj_matrix

M=[0]*(len(DNA_list_join)//2)
agj_matrix_np=np.array(agj_matrix)

def gluing(s1,s2):
    p=max_overlap(s1,s2)
    return s1+s2[p::]

def recurse(id,M):
    k=0
    x=id
    while M[id]!=x:
        k+=1
        id=M[id]
    return [id,k]

def recurse1(id,M):
    return M[id]

def recurse2(id,k,M):
    J=[id]*(k+1)
    p=1
    while k>0:
        J[p]=recurse1(id,M)
        p+=1
        k-=1
        id=M[id]
    return J

def process_matr(matr,M):

    for i in range(len(M)):
        if np.max(agj_matrix_np[i])==0 and np.max(agj_matrix_np[:,i])==0:
            agj_matrix_np[i] = [-1] * len(M)
            agj_matrix_np[:, i] = [-1] * len(M)
            M[i]=i
       # print(agj_matrix_np,M,len(M))
        else:
            j, k = np.where(agj_matrix_np == np.max(np.max(agj_matrix_np, axis=0)))
        #print(j, k)
            M[j[0]] = k[0]
            agj_matrix_np[j[0]] = [-1] * len(M)
            agj_matrix_np[:, k[0]] = [-1] * len(M)