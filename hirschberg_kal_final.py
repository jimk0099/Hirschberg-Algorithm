#===== Imports =====
import numpy as np
import argparse


#===== FUNCTIONS =====
def EnumerateAlignments(A, B, F, W, Z, WW, ZZ):
    # A, B: input strings
    # F: Matrix
    # W, Z: current states of strings
    # WW, ZZ: strings to be returned
    
    i = len(A)
    j = len(B)
    
    # end of recursion -> append W and Z
    if i == 0 and j == 0:
        WW.append(W)
        ZZ.append(Z)
        
    # diagonal step
    if i > 0 and j > 0:
        md = Compare(A[i-1], B[j-1])
        if F[i][j] == (F[i-1][j-1] + md):
            EnumerateAlignments(A[:i-1], B[:j-1], F, A[i-1] + W, B[j-1] + Z, WW, ZZ)  
    
    # upward step
    if i > 0 and F[i][j] == (F[i-1][j] + args.gap):
        if args.l == False:
            EnumerateAlignments(A[:i-1], B, F, A[i-1] + W, '-' + Z, WW, ZZ)
        else:
            EnumerateAlignments(A[:i-1], B, F, A[i-1] + W, '-\n' + Z, WW, ZZ)
        
    # leftward step
    if j > 0 and F[i][j] == (F[i][j-1] + args.gap):
        if args.l == False:
            EnumerateAlignments(A, B[:j-1], F, '-' + W, B[j-1] + Z, WW, ZZ)
        else:
            EnumerateAlignments(A, B[:j-1], F, '-\n' + W, B[j-1] + Z, WW, ZZ)
        
    return (WW, ZZ)
    

def Compare(x,y):
    # just a simple compare function xd
    if(x == y):
        return args.match
    else:
        return args.differ
    

def ComputeAlignmentScore(A, B, Compare, g):
    # Get the last row of FMatrix without making the whole table
    # A, B: input strings
    # Compare is the function defined above
    # g is an integer given by the user
    
    L = np.zeros(len(B)+1)
    for j in range(len(L)):
        L[j] = j * g
    K = np.zeros(len(B)+1)
    for i in range(1, len(A)+1):
        L,K = K,L                                        #swap current last line with the previous                                    
        L[0] = i * g
        for j in range(1,len(B)+1):
            md = Compare(A[i-1], B[j-1])                    
            L[j] = max(L[j-1]+g, K[j]+g, K[j-1]+md)
    return L


def ComputeMatrix(A, B, Compare, g):
    # Compute whole table
    # A, B: input strings
    # Compare is the function defined above
    # g is an integer given by the user
    
    F = np.ndarray((len(A)+1, len(B)+1))
    L = np.zeros(len(B)+1)
    for j in range(len(L)):
        L[j] = j * g
    F[0] = L
    K = np.zeros(len(B)+1)
    for i in range(1, len(A)+1):
        L,K = K,L
        L[0] = i * g
        for j in range(1,len(B)+1):
            md = Compare(A[i-1], B[j-1])
            L[j] = max(L[j-1]+g, K[j]+g, K[j-1]+md)
        F[i] = L
    return F


def NeedlemanWunsch(A, B):
    # A, B: input strings
    # With the use of Hirschberg method this function 
    # is called only for len(A) == 1 or len(B) == 1
    F_matrix = ComputeMatrix(A, B, Compare, args.gap)
    return EnumerateAlignments(A, B, F_matrix, '', '', [], [])
 

def Hirschberg(A, B, WW, ZZ):
    # A, B: input strings
    # WW, ZZ: strings to be returned
    
    W = ''          
    Z = ''
    if len(A) == 0:                                     # base of recursion
        for k in range(len(B)):                         # fill with - untill the end of the smaller string
            W = W + '-'
            Z = Z + B[k]
        return([W], [Z])
    elif len(B) == 0:                                   # base of recursion
        for k in range(len(A)):                         # fill with - untill the end of the smaller string
            Z = Z + '-'
            W = W + A[k]
        return([W], [Z])
    elif len(A) == 1 or len(B) == 1:                    # we can call Needleman-Wunsch as we have small input
        #print('hello mum', A, B)
        temp = NeedlemanWunsch(A, B)
        for k in range(len(temp[0])):
            WW.append(temp[0][k])
            ZZ.append(temp[1][k])
    else:
        i = len(A) // 2                                                             # round down
        ScoreL = ComputeAlignmentScore(A[:i], B, Compare, args.gap)
        ScoreR = ComputeAlignmentScore(A[i:][::-1], B[::-1], Compare, args.gap)[::-1]
        Score = ScoreL + ScoreR
        J = np.where(Score == max(Score))               
        J = J[0]                                                                    # get list
        for j in J:                                                                 # we need a for loop as we might have more than one equivalent solutions  
            if args.t == True:
                print('%d, %d' %(i, j))                                             # print the t's as well  
            (WWl, ZZl) = Hirschberg(A[:i],B[:j], [], [])                            # call Hirschberg - divide and conquer
            (WWr, ZZr) = Hirschberg(A[i:len(A)],B[j:len(B)], [], [])
            for m in range(len(WWr)):
                for n in range(len(WWl)):   
                    WWtemp = WWl[n] + WWr[m]
                    ZZtemp = ZZl[n] + ZZr[m]
                    indices1 = [index for index, x in enumerate(WW) if x == WWtemp]         # indices are used to prevent printing out the same solutions multiple times  
                    indices2 = [index for index, x in enumerate(ZZ) if x == ZZtemp]         # that is why a function UpdateAlignments is defined 
                    (WW, ZZ) = UpdateAlignments(indices1, indices2, WW, ZZ, WWtemp, ZZtemp)
    return (WW, ZZ)       # tuple of str lists


def UpdateAlignments(list1, list2, WW, ZZ, WWtemp, ZZtemp):
    # Recursion function that updates the output properly 
    # list1, list2: the indices1 and indices2 lists that are made above
    
    if(list1 == [] or list2 == []):                                                 # base of recursion                             
        WW.append(WWtemp)
        ZZ.append(ZZtemp)
    elif(len(list1) < len(list2)):                                                  # if list1 is smaller that list2
        if list1[0] not in list2:                                                   # search for current item in list 1 and 'repeat' for the other elements with recursion 
            UpdateAlignments(list1[1:], list2, WW, ZZ, WWtemp, ZZtemp)
    else:
        if list2[0] not in list1:                                                   # if list2 bigger or equal with list1
            UpdateAlignments(list1, list2[1:], WW, ZZ, WWtemp, ZZtemp)
    return (WW, ZZ)
    

# Function to convert  
def listToString(s): 
    
    # initialize an empty string
    str1 = '' 
    
    # return string  
    return (str1.join(s))


def PrintRows(WW, ZZ):
    # prints the results when l is true
    for n in range(len(WW)):
        for line1, line2 in zip(WW[n].splitlines(), ZZ[n].splitlines()):
            if (line1 == line2):                                             # exactly same sentence
                W_symbol = '='
                Z_symbol = '='
            else:
                W_symbol = '<'
                Z_symbol = '>'
                
            print('%s %s' %(W_symbol, line1))
            print('%s %s' %(Z_symbol, line2))
            
            
#===== ARGUMENTS HANDLING =====
# manual for user
parser = argparse.ArgumentParser(description='Compute the optimal alignment of two sequences according to given weights')

parser.add_argument('-t', action='store_true', help="show the (i, j) pairs in the Hirschberg algorithm's execution (default: do not show)")
parser.add_argument('-f', action='store_true', help='a and b are names of files containing the sequences (default: a and b are the sequences)')
parser.add_argument('-l', action='store_true', help='interpret the lines of the files specified by a and b as the input sequences (default: the sequences may span multiple lines)')

parser.add_argument('gap', type=int, help='gap')
parser.add_argument('match', type=int, help='match')
parser.add_argument('differ', type=int, help='differ')

parser.add_argument('a', type=str, help='the first sequence to be aligned')
parser.add_argument('b', type=str, help='the second sequence to be aligned')

args = parser.parse_args()


#===== MAIN =====
def main():
    
    # we have two sequencies as input
    if args.l == False:
        
        # A and B are input strings
        if args.f == False:
            A_content = args.a
            B_content = args.b
        
        # A and B are names of the files that contains the input strings
        else:
            fd1 = open(args.a, "r")                                      # open files for reading
            fd2 = open(args.b, "r")
            A_content = fd1.read()                                  # store content to variables
            B_content = fd2.read()
            fd1.close()                                             # close files
            fd2.close()
      
        results = Hirschberg(A_content, B_content, [], [])          # results is a tuple with lists that contain the possible arrangements
        for i in range(len(results[0])):                            # for all the possible arrangements   
            print(results[0][i])                                    # print the fixed sequences
            print(results[1][i])
            if i != len(results[0]):                                # print an endline at the end for aesthetic reasons
                print('')
    
    else:
        fd1 = open(args.a, "r")
        fd2 = open(args.b, "r")
        
        A_content = [A_line for A_line in fd1]                      # read file line by line
        B_content = [B_line for B_line in fd2]
        
        lenA = len(A_content)                                       # number of rows            
        lenB = len(B_content)
        
        if A_content[lenA - 1][len(A_content[lenA - 1]) - 1] != '\n':       # if last char is not end of line then add it
            A_content[lenA - 1] = A_content[lenA - 1] + '\n'

        if B_content[lenB - 1][len(B_content[lenB - 1]) - 1] != '\n':
            B_content[lenB - 1] = B_content[lenB - 1] + '\n'

        (WW, ZZ) = NeedlemanWunsch(A_content, B_content)                    # Ideally here we would call Hirschberg function but it did not work very well
        PrintRows(WW, ZZ)

       
# execute main on call  
if __name__ == "__main__":
    main()
