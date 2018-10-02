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

def matrix_solve(expressions, equals, variables):
    """
    Solves a system of equations put into Ax + By + ... Cz = D form
    All variables must be in the same order and all expression and variables must be strings.
    The solutions should be a float or int, and the variables must be single-character
    Returns the solutions in order of variables parameter
    """
    expressions = list(expressions)
    equals = list(equals)
    variables = list(variables)
    new_exprs = []
    for i, expr in enumerate(expressions):
        expr=expr.replace(" ","")
        expressions[i] = expr
        terms = []
        pre_index = 0
        for var in variables:
            index = expr.index(var)
            terms.append(expr[pre_index:index+1])
            pre_index = index+1
        new_exprs.append(terms)

    coefficients = []
    for expression in new_exprs:
        expr_coeff = []
        for term in expression:
            for var in variables:
                if var in term:
                    num = term.replace(var,"")
                if  num == '' or  num == '+':
                    num = 1
                elif num == '-':
                    num = -1
            num = float(num)
            expr_coeff.append(num)
        coefficients.append(expr_coeff)
    
    matrixForm = matrix(coefficients)
    equalsMatrix = []
    for a in equals:
        equalsMatrix.append([a])
    equalsMatrix = matrix(equalsMatrix)
    solutions = matrixForm.inv().dot(equalsMatrix)
    final = []
    for sol in solutions:
        final.append(sol[0])
    
    return tuple(final)