
# coding: utf-8

# # Tugas Bioinformatika
# ## Needleman Wunsch Algorithm
# 
# ### Nadya Avirianta S
# ### 16/394096/PA/17187

# Link github: https://github.com/nadyavirianta/bioinformatika

# Needleman Wunsch algorithm digunakan untuk protein alignment. 

# In[1]:

import pandas as pd 
import numpy as np


# Mengimport data BlOSUM62 dari CSV file

# In[2]:

blosum = pd.read_csv('blosum62.csv')
blosum


# Fungsi needleman dengan parameter strings protein yang akan dialign 

# In[6]:

def needleman(a,b):
    #buat matrix dengan dimensi strings a+1 dan b+1
    row = len(b)+1
    col = len(a)+1
    #Matrix F untuk menyimpan value, T untuk menyimpan arah, B untuk menyimpan Blosum Score.
    F = np.zeros(shape=(row,col), dtype=np.int)
    T = np.full(shape=(row,col), fill_value=" ", dtype=np.str)
    B = np.zeros(shape=(row,col), dtype=np.int)
    #d adalah score gap
    d = 6
    #pengisian matrix 
    F[0][0]=0
    for i in range(1,row):
        F[i][0] = F[i-1][0]-d
    for j in range(1,col):
        F[0][j] = F[0][j-1]-d    
    for i in range (1,row):
        for j in range (1,col):
            match = F[i-1][j-1]+blosum[a[j-1]][b[i-1]]
            gapa = F[i][j-1]-d
            gapb = F[i-1][j]-d
            #memilih value maksimum dari vertical, horizontal, dan diagonal dan beri panah yang sesuai pada matrix T
            F[i][j] = max(match,gapa,gapb) 
            if max(match,gapa,gapb) ==match:
                T[i][j] = '↖'
            elif max(match,gapa,gapb) == gapa:
                T[i][j] ='←'
            elif max (match,gapa,gapb)==gapb:
                T[i][j] ='↑'
    #inisialisasi matrix T
    T[0][0] = "-"
    for i in range(1, row):
        T[i][0] = "↑"
    for j in range(1, col):
        T[0][j] = "←"
    #pengisian matrix B dengan blosumscore
    B[0][0]=0
    for i in range(1,row):
        B[i][0] = B[i-1][0]-d
    for j in range(1,col):
        B[0][j] = B[0][j-1]-d 
    for j in range(1,col):
        for i in range(1,row):
            B[i][j] = blosum[a[j-1]][b[i-1]]
    #traceback alignment mengikuti arah panah
    aalign = ''
    balign = ''
    i = len(b)
    j = len(a)
    score = 0
    while i>0 or j>0:
        if T[i][j] == '↖':
            aalign = a[j-1]+aalign
            balign = b[i-1]+balign
            score=score+B[i][j]
            i = i-1
            j = j-1         
        elif T[i][j] == '←':
            aalign = a[j-1]+aalign
            balign = '_'+balign
            score=score-d
            j = j-1
        elif T[i][j] == '↑':
            aalign = '_'+aalign
            balign = b[i-1]+balign
            score=score-d
            i=i-1
        
        elif j>0:
                while j>0:
                    aalign = a[j-1]+aalign
                    balign = '_'+balign
                    score=score-d
                    j=j-1
        elif i>0:
                while i>0:
                    aalign = '_'+aalign
                    balign = b[i-1]+balign
                    score=score-d
                    i=i-1
        else:
            break
        
    print(aalign)
    print(balign)
    print("SCORE:",score)
    print("\n")
    print(F)
    print("\n")
    print(T)
    print("\n")
    print(B)


# In[7]:

needleman('MNALQM','NALMSQA')


# ## Contoh penggunaan fungsi

# Mengalign protein QALVAYA dan NALWVAYMA

# In[9]:

needleman('QALVAYA','NALWVAYMA')


# Contoh lain membaca protein dari text, yaitu protein HEAGAWGHEE dan PAWHEAE

# In[14]:

with open('protein.txt', 'r') as myfile:
    data=myfile.readlines()


# In[15]:

data = [x.strip() for x in data]


# In[16]:

x =data[0]
y=data[1]


# In[17]:

needleman(x,y)

