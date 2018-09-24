from matrix_library import matrix as matrix

def full_matrix(rows, columns, fill = 0, rand = False, randrange = (0,10)):
    """
    Gives a matrix of a set size
    N x M and can be filled with a value
    or a random integer value for each entry 
    """
    mtrx = [[fill for j in range(columns)] for n in range(rows)]
    if rand:
        from random import randint
        mtrx = [[randint(randrange[0],randrange[1]) for j in range(columns)] for n in range(rows)]
    return matrix(mtrx)

def matrix_range(beg, end, step = 1):
    """
    Makes a matrix with a range
    of beg to end stepping by step 
    with 1 x M size
    """
    i = beg
    mtrx = []
    while i < end:
        mtrx.append(i)
        i += step
    return matrix([mtrx])

def spaced_matrix(beg,end,space):
    """
    Matrix with values evenly spaced from
    beg to end including end
    with values spaced space apart
    """
    step = (end - beg)/(space-1)
    i = beg
    mtrx = []
    while i <= end:
        mtrx.append(i)
        i += step
    return matrix([mtrx])