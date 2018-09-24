from matrix_functions import ShapeError as ShapeError
from matrix_functions import permutations as permutations
from matrix_functions import parity as parity

class matrix(object):
    """
    An rectangular array of numbers.
    Class for mathematical matrix
    Basically a 2D array
    Format example:
    matrix([[a_11, a_12]
            [a_21, a_22]
            [a_31, a_32]])
    """
    def __init__(self,elements,is_identity = False):
        # Declared as: matrix([[row 1],[row 2], [row 3]])
        self.elements = elements
        self.rows = len(self.elements)
        self.columns = len(self.elements[0])

        # Each row must have the same length
        for num, lst in enumerate(self):
            # Iterates through each row to check that they have the same number of columns
            
            col_len = len(lst)
            
            if col_len < self.columns:
                # Raises ShapeError if column has too few elements
                raise ShapeError('Column ' + str(num) + ' does not have enough elements. (Expected '+str(self.columns)+ ')')
           
            elif col_len > self.columns:
                # Raises ShapeError if column has too many elements
                raise ShapeError('Column ' + str(num) + ' has too many elements. (Expected '+str(self.columns)+ ')')

        self.dim = (self.rows, self.columns)

        self.is_square = self.dim[0] == self.dim[1]
        if self.is_square and not is_identity:
            # Creates identity matrix for all square matrices

            self.identity = [[0 for o in range(self.columns)] for p in range(self.rows)]
            for l in range(self.rows):
                self.identity[l][l] = 1

            self.identity = matrix(self.identity,is_identity=True)
            # Sets is_identity = True to prevent recursion,
            # As ID matrices do not need ID matrices

    def __mul__(self, c):
        """
        Scalar Multiplication
        Returns self * c
        """
        # Must be self * c, cannot be c * self
        if type(c) == int or type(c) == float:
            return matrix([[(int(c*1e16) * int(j*1e16))/1e32 for j in a] for a in self])
        else:
            raise TypeError("unsupported operand type(s) for *: 'matrix' and "+"'"+str(type(c).__name__)+"'")

    def __imul__(self,c):
        """
        Scalar Multiplication in place
        Returns self * c
        """
        return self * c
    
    def __add__(self, b):
        """
        Matrix Addition
        Returns self + b
        """

        if type(b).__name__ != 'matrix':
            raise TypeError("unsupported operand type(s) for +: 'matrix' and "+"'"+str(type(c).__name__)+"'")
        
        elif self.dim != b.dim:
            # Must operate on matrices of the same size
            raise ShapeError('Matrix dimensions not equal. expected '+'{}'.format(self.dim))
        
        else:
            return matrix([[(int(self[i][j]*1e16) + int(b[i][j]*1e16))/1e16 for j in range(self.columns)] for i in range(self.rows)])
    
    def __iadd__(self,b):
        """
        Matrix Addition in place
        Returns self + b
        """
        return self + b

    def __sub__(self, b):
        """
        Matrix Subtraction
        Returns self - b
        """
        if type(b).__name__ != 'matrix':
            raise TypeError("unsupported operand type(s) for -: 'matrix' and "+"'"+str(type(c).__name__)+"'")
        
        elif self.dim != b.dim:
            # Must operate on matrices of the same size
            raise ShapeError('Matrix dimensions not equivalent. expected '+'{}'.format(self.dim))
        
        else:
            return matrix([[(int(self[i][j]*1e16) - int(b[i][j]*1e16))/1e16 for j in range(self.columns)] for i in range(self.rows)])

    def __isub__(self,b):
        """
        Matrix Subtraction in place
        Returns self - b
        """
        return self - b

    def Hadamard(self, b):
        """
        Returns the Hadamard/Schur/entrywise product of self and b
        Takes two matrices of the same dimensions, and produces 
        another matrix where each element i,j is the product 
        of elements i,j of the original two matrices.
        Wikipedia article: https://en.wikipedia.org/wiki/Hadamard_product_(matrices)
        """
        if type(b).__name__ != 'matrix':
            raise TypeError("unsupported operand type(s) for +: 'matrix' and "+"'"+str(type(c).__name__)+"'")

        if self.dim != b.dim:
            # Must operate on matrices of the same size
            raise ShapeError('Matrix dimensions not equal. expected '+'{}'.format(self.dim))

        else:
            return matrix([[(int(self[i][j]*1e16) * int(b[i][j]*1e16))/1e32 for j in range(self.columns)] for i in range(self.rows)])

    def dot(self, c):
        """
        Returns the dot product of two matrices.
        Also known as Matrix Multiplication
        Wikipedia article: https://en.wikipedia.org/wiki/Dot_product#Algebraic_definition
        Uses rounding to get more accurate results
        """
        if type(c).__name__ != 'matrix':
            TypeError("unsupported operand type(s) for matrix.dot: 'matrix' and "+"'"+str(type(c).__name__)+"'")
        
        if self.dim[1] != c.dim[0]:
            # Must operate on matrices of the same size
            raise ShapeError('Matrix dimensions not compatible. (' + str(self.dim[1]) + ' != ' + str(self.dim[0]) + ')')
        
        new_matrix = [[0 for j in range(c.columns)] for i in range(self.rows)]
        for a in range(self.rows):
            for b in range(c.columns):
                element_sum = 0
                for d in range(c.rows):
                    
                    element = (int(self[a][d]*1e16)*int(c[d][b]*1e16))
                    # Counteracts rounding error by turning values into integers
                    element_sum += element

                element_sum = float(element_sum)
                # And then back to floats
                element_sum /= 1e32

                new_matrix[a][b] = element_sum
        
        return matrix(new_matrix)

    def det(self):
        """
        Finds the determinant of a N x N matrix using the Leibniz formula.
        Wikipedia article: https://en.wikipedia.org/wiki/Leibniz_formula_for_determinants
        """
        if not self.is_square:
            raise ShapeError('Matrix is not square')

        if self.dim[0] == 2:
            # Easy special case of 2 x 2 matrix
            det = (int(self[0][0]*1e16)*int(self[1][1]*1e16)) - (int(self[0][1]*1e16)*int(self[1][0]*1e16))
            # Prevents round off error
            det = det/100000000000000000000000000000000
            if int(det) == det:
                det = int(det)
            return det
            
        elif self.dim[0] > 2:
            # Finds all possible permutations in the permutation group S_n, where n = dimensions of the matrix
            base_set = list(range(self.rows))
            permutation_lst = permutations(base_set)

            determinant = 0
            for p in permutation_lst:
                sgn = parity(list(p)) # Find sign of permutation from parity
                
                # Computes determinant using formula
                run_mul = sgn # Prevents roundoff error
                for i in range(self.rows):
                    run_mul *= int(self[i][p[i]]*1e16)
                determinant += run_mul
            
            determinant /= 1e16**self.rows
            if int(determinant) == determinant:
                determinant = int(determinant)

            return determinant
    
    def inv(self):
        """
        Uses Cramer's Rule to compute the inverse of a matrix.
        Uses set rule for 2 x 2 matrices
        Uses soft code for N x N matrices
        Cramer's Rule:  1/det(A) * adj(A) = A^-1
        where A is a is_square matrix that has a non-zero determinant.
        adj(A) is the transpose of the matrix of cofactors of A, known as the adjugate matrix
        Wikipedia article: https://en.wikipedia.org/wiki/Invertible_matrix#Analytic_solution
        """
        if not self.is_square:
            raise ShapeError('Matrix is not square')

        elif self.det() == 0:
            raise ZeroDivisionError('Matrix determinant = 0')

        elif self.rows == 1:
            # Special case
            inverse = matrix([[1/self[0][0]]])

        elif self.rows == 2:
            # Special case
            inverse = matrix([[self[1][1],-1*self[0][1]],[-1*self[1][0],self[0][0]]])*(1/self.det())
        
        elif self.rows > 2:
            # Find matrix of minors
            # Wikipedia article: https://en.wikipedia.org/wiki/Minor_(linear_algebra)
            matrix_of_minors = [[0 for j in range(self.columns)] for n in range(self.rows)]
            for r in range(self.rows):
                minus_row = self.delrow(r) # Delete row
                for c in range(self.columns):
                    minus_row_and_col = minus_row.delcol(c) #Delete column
                    matrix_of_minors[r][c] = minus_row_and_col.det()

            # Calculate comatrix
            comatrix = matrix(matrix_of_minors) # Create Matrix with the matrix of minors
            for a, row in enumerate(comatrix):
                for b, col in enumerate(row):
                    if a%2 == 0:
                        if (b+1)%2 == 0:
                            comatrix[a][b] = -col
                    if a%2 != 0:
                        if (b)%2 == 0:
                            comatrix[a][b] = -col

            # Find transpose
            adjugate = comatrix.transpose()

            # Adjugate * 1/det
            # Wikipedia article: https://en.wikipedia.org/wiki/Adjugate_matrix
            det = self.det()
            inverse = matrix([[(int(j*1e16)/int(det*1e16)) for j in a] for a in adjugate])
       
        return inverse

    def transpose(self):
        """
        Find transpose of matrix
        Denoted as A^t
        Matrix is flipped over primary diagonal axis
        matrix([[2,3,4],
                [1,2,3]]).transpose() = 
        matrix([[2,1]
                [3,2]
                [4,3]])
        Wikipedia article: https://en.wikipedia.org/wiki/Transpose
        """
        new_matrix = [[0 for j in range(self.rows)] for n in range(self.columns)]       
        for i in range(self.rows):
            for j in range(self.columns):
                new_matrix[j][i] = self[i][j]
        return matrix(new_matrix)

    def delrow(self,n):
        """
        Deletes a row from a matrix
        """
        return matrix(self.elements[:n] + self.elements[n+1:]) 

    def delcol(self,n):
        """
        Deletes a column from a matrix
        """
        return matrix([r[:n] + r[n+1:] for r in self])
    
    def getcol(self,i,j = None):
        """
        Get a column from a matrix
        Necesary because matrix is stored
        as rows in lists
        """
        if j == None:
            j = i+1
        return matrix([r[i:j] for r in self])

    def __str__(self):
        """
        Returns:
        [row 1]
        [row 2]
        [row 3]
        """
        for i, row in enumerate(self): # Turns unnecesary floats into integers
            for j, element in enumerate(row):
                if int(element) == element:
                    self[i][j] = int(element)
        
        string = []
        for row in self:
            string.append('{}'.format(row).replace(',',' '))
            string.append('\n')
        string.pop(-1) # Removes final new line
        return ''.join(string)
    
    def __repr__(self): 
        # Python object method
        # Computes the “official” string representation of a matrix
        return "matrix(" + '{}'.format(self.elements)+")"

    def __iter__(self):
        # Python object method
        # Returns a new iterator object that can iterate over all the objects in the matrix
        self.n = 0
        return self

    def __next__(self):
        if self.n < self.rows:
            result = self.elements[self.n]
            self.n += 1
            return result
        else:
            raise StopIteration

    def __getitem__(self, key):
        # Python object method
        # Called to implement evaluation of self[key]
        if type(key)==slice:
            return matrix(self.elements[key])
        else:
            return self.elements[key]
    
    def __len__(self):
        # Python Object method
        # Called to implement the built-in function len()
        return self.rows
    
    def __setitem__(self, key, val):
        # Python object method
        # Called to implement assignment of self[key]
        if type(val) != int and type(val) != float and type(val) != list:
            raise TypeError('Item must be of type integer of float')
        if type(val) == list and len(val) != self.columns:
            raise ShapeError('Row must have same columns as original matrix')
        else:
            if type(key) == slice:
                if key.step != None:
                    raise IndexError('Index step is not supported by matrix')
                new_row = 0
                dist = key.stop - key.start - 1
                new_row = self.rows - dist
            else:
                new_row = self.rows
            self.elements[key] = [val]
            self.rows = new_row
            return self

    def __contains__(self,val):
        # Python Object method
        # Called to implement membership test operators. Should return true if item is in self, false otherwise
        result = False
        for a,b in enumerate(self):
            if val in b:
                result = True
                break
        return result

    def __neg__(self):
        # -self
        for i in range(self.rows):
            for j in range(self.columns):
                self[i][j] *= -1
        return self

    def __abs__(self):
        # abs(self)
        for i in range(self.rows):
            for j in range(self.columns):
                self[i][j] = abs(self[i][j])
        return self
    
    def __eq__(self,v):
        # self == 
        return self.elements == v.elements
    
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