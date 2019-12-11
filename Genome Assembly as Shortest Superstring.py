import difflib
from operator import is_not
from functools import partial
import numpy as np
import copy
import itertools

DNA_list=[]

with open('Overlap_Graphs_data1.txt','r') as data_Rosa:
    for line in data_Rosa:
        DNA_list.append(line.strip())

def list_filter_None(L):
    return list(filter(partial(is_not, None), L))

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

# заполнение матрицы смежности

agj_matrix=[[0]*(len(DNA_list_join)//2) for i in range(len(DNA_list_join)//2)]

for num1, elem1 in enumerate(DNA_list_join):
    if elem1.isalnum():
        for num2, elem2 in enumerate(DNA_list_join):
            if elem2.isalnum() and elem1!=elem2:
                agj_matrix[num1//2][num2//2]=max_overlap(elem1,elem2)
           
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

M=[0]*(len(DNA_list_join)//2)
agj_matrix_np=np.array(agj_matrix)

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
        #print(agj_matrix_np)
        #print(M)
        
 # разобьем последовательность склеек на циклы
K=[]
for id in M:
    if id==recurse(M[id],M)[0] and recurse(id,M)[1]>0:
        K.append(recurse2(id,recurse(id,M)[1],M))

    elif id==M[id]:
        K.append(DNA_list_join[id*2+1])

#склеиваем по выделенным циклам
G=copy.deepcopy(K)
for num, elem in enumerate(K):
    if type(elem)==list:
        H = copy.deepcopy(DNA_list_join)
        for i in range(len(elem)-1):
            s=gluing(H[elem[i]*2+1],H[M[elem[i]]*2+1])
            H[elem[i+1]*2+1]=s
        G[num][i]=s

K=list(map(lambda x: sorted(x) if type(x)==list else x, K))
K_unic=sorted(K, key=lambda x: str(x), reverse=True)
K_unic=list(k for k,_ in itertools.groupby(K_unic))
U=[]

for id, elem in enumerate(K_unic):
    if type(elem) == list:
        U.append(list(G[t][len(elem)-2] for t in elem))
    else: U.append(elem)

#FINAL!!!!!!!!!!

STRING_FOR_GLUING=list(map(lambda x: min(x, key =lambda x: len((x))) if type(x)==list else x ,U))

#ОСТАЛОСЬ СКЛЕИТЬ

# заполнение матрицы смежности

agj_matrix=[[0]*len(STRING_FOR_GLUING) for i in range(len(STRING_FOR_GLUING))]

for num1, elem1 in enumerate(STRING_FOR_GLUING):
    for num2, elem2 in enumerate(STRING_FOR_GLUING):
        if elem2.isalnum() and elem1!=elem2:
            agj_matrix[num1][num2]=max_overlap(elem1,elem2)

### склейка строк
## полное назначение

M=[0]*(len(STRING_FOR_GLUING))
agj_matrix_np=np.array(agj_matrix)
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
   
# разобьем последовательность склеек на циклы
K=[]
for id in M:
    if id==recurse(M[id],M)[0] and recurse(id,M)[1]>0:
        K.append(recurse2(id,recurse(id,M)[1],M))

    elif id==M[id]:
        K.append(STRING_FOR_GLUING[id])

G=copy.deepcopy(K)
for num, elem in enumerate(K):
    if type(elem)==list:
        H = copy.deepcopy(STRING_FOR_GLUING)
        for i in range(len(elem)-1):
            s=gluing(H[elem[i]],H[M[elem[i]]])
            H[elem[i+1]]=s
        G[num][i]=s

K=list(map(lambda x: sorted(x) if type(x)==list else x, K))
K_unic=sorted(K, key=lambda x: str(x), reverse=True)
K_unic=list(k for k,_ in itertools.groupby(K_unic))
U=[]
for id, elem in enumerate(K_unic):
    if type(elem) == list:
        U.append(list(G[t][len(elem)-2] for t in elem))
    else: U.append(elem)

STRING_FOR_GLUING=list(map(lambda x: min(x, key =lambda x: len((x))) if type(x)==list else x ,U))

print(STRING_FOR_GLUING)
string=''
for i in STRING_FOR_GLUING:
    string+=i
print(string)
