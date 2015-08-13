
def notNull(X):
    """
    REAL FUNCTION notNull (X). This function can be used to avoid "divide
    by zero" errors in division. The function result is defined as notNull =
    X if X is not 0 and notNull = 1 when X is 0.
    """
    return X if (X <> 0.) else X



def INSW(X1, X2, X3):
    """
    input switch relay
    """
    return X2 if X1 < 0 else X3

       
       
def REAAND(X1, X2):
    """
    Returns 1.0 if both input values are positive, otherwise Y=0.0.
    """
    return 1. if (X1 > 0 and X2 > 0) else 0.   


    
