# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 20:55:32 2022

"""

import numpy as np
import sys
import argparse



def compare(a, b, m, d):
    # a: integer, the first value to be compared
    # b: integer, the second value to be compared
    # m: positive integer, the weight being returned when a == b
    # d: negative integer, the weight being returned when a != b
    return m * (a == b) + d * (a != b) 



def enumerate_alignments(A, B, F, W, Z, WW, ZZ):
    # A: string or list of strings, the first sequence to be aligned
    # B: string or list of strings, the second sequence to be aligned
    # F: 2D-ndarray, the alignment matrix
    # W: string, the alignment of A being constructed
    # Z: string, the alignment of B being constructed
    # WW: list, the list of currently constructed alignments of A
    # ZZ: list, the list of currently constructed alignments of B
    
    
    # i = |A|
    # j = |B|
    i = len(A)
    j = len(B)

    # If |A| == 0 and |B| == 0 the construction of the current alignment is finished;
    # Add them to their respective lists and return
    if (i == 0) and (j == 0):
        WW.append(W)
        ZZ.append(Z)
        return 
    
    
    if (i > 0) and (j > 0):
        
        # Check wether there is a match or a difference between the last elements
        # of the sequences.
        md = compare(A[i - 1], B[j - 1], m, d)
    
        # Check if the current value F[i, j]  from its upper-left diagonal 
        # element
        if F[i, j] == (F[i - 1, j - 1] + md):
            W_0 = A[i - 1] + W
            Z_0 = B[j - 1] + Z
            enumerate_alignments(A[0 : i - 1], B[0 : j - 1], F, W_0, Z_0, WW, ZZ)
    
    # If only A has non zero length check the element that is located above
    # the current position in the matrix F
    if (i > 0) and (F[i , j] == (F[i - 1, j] + g)):
        
        # Append the current element of A at the start of the current alignment of A
        W_1 = A[i - 1] + W
        
        # Append a '-' at the start of the current alignment of B since it has 
        # no more elements.
        if l == True:
            Z_1 = '-\n' + Z
        else:
            Z_1 = '-' + Z
            
        enumerate_alignments(A[0 : i - 1], B, F, W_1, Z_1, WW, ZZ)
    
    
    if (j > 0) and (F[i , j] == (F[i, j - 1] + g)):
        if l == True:
            W_2 = '-\n' + W
        else:
            W_2 = '-' + W
        Z_2 = B[j - 1] + Z
        enumerate_alignments(A, B[0 : j - 1], F, W_2, Z_2, WW, ZZ)



def compute_alignment_score(A, B, compare, g):
    # A: string or list of strings, the first sequence to be aligned
    # B: string or list of strings, the second sequence to be aligned
    # compare: function, a function calculating the match or difference score
    # between sequence items.
    # g: negative integer, the negative gap score
    
    # i = |A|
    # j = |B|
    i = len(A)
    j = len(B)
    
    L = np.zeros(j + 1, np.int64)
    for n in range(0, j + 1):
        L[n] = n * g
    
    for n in range(1, i + 1):
        K = L.copy()
        L[0] = n * g
        for k in range(1, j + 1):
            md = compare(A[n - 1], B[k - 1], m, d)
            L[k] = max(L[k - 1] + g, K[k] + g, K[k - 1] + md)
    return L



def compute_alignment_matrix(A, B, compare, g):
    # A: string or list of strings, the first sequence to be aligned
    # B: string or list of strings, the second sequence to be aligned
    # compare: function, a function calculating the match or difference score
    # between sequence items.
    # g: negative integer, the negative gap score
    
    # i = |A|
    # j = |B|
    i = len(A)
    j = len(B)
    
    F = np.zeros((i + 1, j + 1), np.int64)
    
    L = np.zeros(j + 1, np.int64)
    for n in range(0, j + 1):
        L[n] = n * g
    
    F[0] = L
    
    for n in range(1, i + 1):
        K = L.copy()
        L[0] = n * g
        for k in range(1, j + 1):
            md = compare(A[n - 1], B[k - 1], m, d)
            L[k] = max(L[k - 1] + g, K[k] + g, K[k - 1] + md)
            
        F[n] = L
    return F



def NeedlemanWunsch(A, B):
    # A: string or list of strings, the first sequence to be aligned
    # B: string or list of strings, the second sequence to be aligned
    
    # WW: list, the list of currently constructed alignments of A
    # ZZ: list, the list of currently constructed alignments of B
    WW = []
    ZZ = []
    
    # F: 2D-ndarray, the alignment matrix
    F = compute_alignment_matrix(A, B, compare, g)
    
    # Constructs the alignments of A and B
    enumerate_alignments(A, B, F, '', '', WW, ZZ)
    return (WW, ZZ)



def update_alignments(WW, ZZ, w_element, z_element):
    # Ensure that the pair of alignments is not already included in the list
    for n in range(len(WW)):
        if (WW[n] == w_element) and (ZZ[n] == z_element):
            return
      
    WW.append(w_element)
    ZZ.append(z_element)
    


def Hirschberg(A, B):
    # A: string or list of strings, the first sequence to be aligned
    # B: string or list of strings, the second sequence to be aligned    
    
    # i = |A|
    # j = |B|
    i = len(A)
    j = len(B)
    
    # A is empty
    if i == 0:
        if l == False:
            # For each character of B append a dash to WW
            WW = ['-' * j]
            ZZ = [B]
        else:
            # For every line of B append a line with a single dash to WW
            WW = ['-\n' for n in range(j)]
            ZZ = B
    
    # B is empty
    elif j == 0:
        if l == False:
            # For each character of A append a dash to ZZ
            WW = [A]
            ZZ = ['-' * i]
        else:
            # For every line of A append a line with a single dash to ZZ
            WW = A
            ZZ = ['-\n' for n in range(i)]
        
    # Sequence A or sequence B has a single element
    elif (i == 1) or (j == 1):
        # The size of the alignment matrix is small and NeedlemanWunch algorithm
        # can be used without requiring too many resources 
        (WW, ZZ) = NeedlemanWunsch(A, B)
    
    # Both sequence A and B have multiple elements
    else:
        i = len(A) // 2
        S_l = compute_alignment_score(A[0 : i], B, compare, g)
        S_r = compute_alignment_score(A[i : len(A)][::-1], B[::-1], compare, g)
        S = S_l + S_r[::-1]
        J = np.where(S == max(S))
        WW = []
        ZZ = []
        for j in J[0]:
            if t == True:
                print('%d, %d' %(i, j))
            (WW_l, ZZ_l) = Hirschberg(A[0 : i], B[0 : j])
            (WW_r, ZZ_r) = Hirschberg(A[i : len(A)], B[j : len(B)])
            for l_ in range(len(WW_l)):
                for r_ in range(len(WW_r)):
                    update_alignments(WW, ZZ, WW_l[l_] + WW_r[r_], ZZ_l[l_] + ZZ_r[r_])
    
    return (WW, ZZ)



# Auxiliary function for the visualisation of the results
def print_result(WW, ZZ, l):
    if l == False:
        for n in range(len(WW)):
            print(WW[n])
            print(ZZ[n])
            print()
    
    else:
        for n in range(len(WW)):
            for line_w, line_z in zip(WW[n].splitlines(), ZZ[n].splitlines()):
                
                # Format the output according to the exercise's requirements
                if (line_w == line_z):
                    symbol_w = '='
                    symbol_z = '='
                else:
                    symbol_w = '<'
                    symbol_z = '>'
                    
                print('%s %s' %(symbol_w, line_w))
                print('%s %s' %(symbol_z, line_z))
            
            
            
# Get the number of arguments given in the command line by the user
argc = len(sys.argv)

# Give a simple reminder to the user about the program's usage
if (argc < 6) or (argc > 9):
    print('Incorrect usage')
    print('Usage: %s [-t] [-f] [-l] g m d a b' %sys.argv[0])
    print('\t -t : (optional argument), show the (i, j) pairs in the Hirschberg algorithm')
    print('\t -f : (optional argument), specifies that the arguments a, b are file names')
    print('\t -l : (optional argument), specifies that the sequences to', end = '')
    print('be aligned are the lines of the files a and b')
    print('\t g : negative integer, gap')
    print('\t m : positive integer, match')
    print('\t d : negative integer, differ')
    print('\t a : string, the first sequence to be aligned')
    print('\t b : string, the second sequence to be aligned')
    exit(1)

# Initialize the optional parameters
t = False
f = False
l = False

# Initialize the required parameters as they are given by the user
g = int(sys.argv[argc - 5])
m = int(sys.argv[argc - 4])
d = int(sys.argv[argc - 3])
A = sys.argv[argc - 2]
B = sys.argv[argc - 1]


# Check for any optional parameters
if (argc - 6) == 3:
    t = True
    f = True
    l = True
   
if (argc - 6) == 2:
    if sys.argv[1] == '-t':
        t = True
    else:
        l = True
    f = True
    
if (argc - 6) == 1:
    if sys.argv[1] == '-t':
        t = True
    else:
        f = True



def main():
    WW = []
    ZZ = []
    
    if f == False:
        (WW, ZZ) = Hirschberg(A, B)
        print_result(WW, ZZ, l)
        
    else:
        f_A = open(A, 'r')
        f_B = open(B, 'r')
        
        if l == False:
            a = f_A.read()
            b = f_B.read()
            (WW, ZZ) = Hirschberg(a, b)
            print_result(WW, ZZ, l)
            
        else:
            a = [A_line for A_line in f_A]
            if a[len(a) - 1][len(a[len(a) - 1]) - 1] != '\n':
                a[len(a) - 1] = a[len(a) - 1] + '\n'

            b = [B_line for B_line in f_B]
            if b[len(b) - 1][len(b[len(b) - 1]) - 1] != '\n':
                b[len(b) - 1] = b[len(b) - 1] + '\n'
                
            (WW, ZZ) = Hirschberg(a, b)
            print_result(WW, ZZ, l)
            
        f_A.close()
        f_B.close()

        
if __name__ == "__main__":
    main()  