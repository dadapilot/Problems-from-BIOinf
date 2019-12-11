from functools import partial
from operator import is_not
import numpy as np

DNA_list=[]

with open('Overlap_Graphs_data1.txt','r') as data_Rosa:
    for line in data_Rosa:
        DNA_list.append(line.strip())

def join_string(L):
    a=0
    M=[None]
    for id, elem in enumerate(L):
        if elem.isalnum() is not True:
            a=0
            M.append(elem)
        else:
            a+=1
            if a>1:
                M[len(M)-1]=M[len(M)-1]+elem
            else: M.append(elem)
    return list(filter(partial(is_not, None), M))
def cons(nplst):
    j=np.where(nplst == np.max(nplst))
    return int(np.array(j[0][0], dtype=np.float))

DNA_list_join=join_string(DNA_list)

for i in DNA_list_join:
    if i[0]=='>':
        DNA_list_join.remove(i)

for id, i in enumerate(DNA_list_join):
    DNA_list_join[id]=list(i)

DNA_list_join=np.array(DNA_list_join)
print(DNA_list_join)

profile_matr=[[0]*(len(DNA_list_join[1,:])+1) for _ in range(4)]
l=['A','C','G','T']

for i in range(4):
    for j in range(len(DNA_list_join[0,:])):
        profile_matr[i][j]=DNA_list_join[:,j].tolist().count(l[i])

profile_matr=np.array(profile_matr)

Consensus=''
for i in range(len(DNA_list_join[1,:])+1):
    k=cons(profile_matr[:,i])
    Consensus+=l[k]

L=['A:','C:','G:','T:']

for i in range(4):
    print(L[i],*profile_matr[i,:])
