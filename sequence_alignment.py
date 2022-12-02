# The Needleman-Wunsch Algorithm/ Sequence Alignment
# ==============================
# CS4436 Final Project
# --- Sanjit Sharma --- Peter Wu --- Stefan Pisic ---


from ctypes import alignment
from re import I
from xml.etree.ElementTree import tostring

import numpy as np

def nw(x, y, blosum62_dict, gap = 1, gap_initiate = -10, gap_extend = -2):
    nx = len(x)
    ny = len(y)
    # Optimal score at every possible character pair
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx * gap, nx + 1) # first col
    F[0,:] = np.linspace(0, -ny * gap, ny + 1) # first row

    # Pointers to run through an optimal aligment.
    # same  size as F
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3 # first col
    P[0,:] = 4 # first row
    

    # Record of whether a gap initiation or extension should happen
    # 1: need gap extension, 0: gap extension
    G  = np.zeros((nx + 1, ny + 1))
    G[0:0] = 1 # first gap


    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            
            # match/mistch cost, according to blosum62
            t[0] = F[i,j] + blosum62(x[i], y[j], blosum62_dict) 

            # gap cost, according to gap_extend/initiate
            if G[i,j+1] == 1:
                 t[1] = F[i,j+1] + gap_extend
            else:
                 t[1] = F[i,j+1] + gap_initiate 

            if G[i+1,j] == 1:
                 t[2] = F[i+1,j] + gap_extend
            else:
                 t[2] = F[i+1,j] + gap_initiate 
   
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2 # gap did not happen, G[i+1, j+1] remains 0
            if t[1] == tmax:
                G[i+1, j+1] = 1 # gap happens
                P[i+1,j+1] += 3
            if t[2] == tmax:
                G[i+1, j+1] = 1 # gap happens
                P[i+1,j+1] += 4
                    
    # retreive the total match score
    alignment_score = F[nx, ny]        
    # Run through optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # String reversal
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]

    # # verify alignment_score
    # rx_gap = 0
    # ry_gap = 0
    # score = 0
    # for i in range(len(rx)):
    #     if rx[i] == '-':
    #         if rx_gap == 0:
    #             score += -10
    #             rx_gap = 1
    #         else:
    #             score += -2
    #     elif ry[i] == '-':
    #         if ry_gap == 0:
    #             score += -10
    #             ry_gap = 1
    #         else:
    #             score += -2
    #     else:
    #         score += blosum62(rx[i], ry[i])
    #         rx_gap = 0
    #         ry_gap = 0

    # print(score)

    return [alignment_score,  '\n'.join([rx, ry])]

# store biosum62 matrix into a hash table (dict)
def score_matrix(proteins=''):

    blosum62_dict = dict() # dict to store matrix   

    file = open('BLOSUM62.txt', 'r')
    line = file.readline()

    # store each line of matrix
    for i in range(20):
        line = file.readline()
        splitted = line.split()

        col = 0 # col index

        for score in splitted[1:]:
  
            row = proteins[i] # row char
            blosum62_dict[row + proteins[col]] = int(score)
            col += 1

    return blosum62_dict

# return score in blosum62
def blosum62(c1, c2, blosum62_dict):

    # pairs are non-ordered, so check twice
    pair1 = c1 + c2 
    pair2 = c2 + c1

    if pair1 in blosum62_dict:
        return blosum62_dict[pair1]
    else:
        return blosum62_dict[pair2]