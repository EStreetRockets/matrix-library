class ShapeError(Exception):
    """
    Exception raised for matrix errors concerning the matrix shape
    """
    def __init__(self,message):
        self.message = message

def permutations(lst):
    """
    Finds all permutations of a set
    """
    # If lst is empty then there are no permutations
    if len(lst) == 0:
        return []

    # If there is only one element in lst then, only one permuatation is possible
    if len(lst) == 1:
        return [lst]

    # Find the permutations for lst if there are more than 1 characters

    l = [] # empty list that will store current permutation

    # Iterate the input(lst) and calculate the permutation
    for i in range(len(lst)):
        m = lst[i]

        # Extract lst[i] or m from the list.  remLst isremaining list
        remLst = lst[:i] + lst[i+1:]

        # Generating all permutations where m is first element
        for p in permutations(remLst):
            l.append([m] + p)
    return l

def parity(lst):
    """
    Finds parity of a list of numbers
    """
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity