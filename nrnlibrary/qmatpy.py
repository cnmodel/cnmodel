"""
Pure Python implementations of some Q-Matrix algorithms.

C.Schmidt-Hieber, University College London

2010-08-23

Equation numbers refer to:
Colquhoun D, Hawkes AG (1995a) 
The principles of the stochastic interpretation of ion channel mechanisms. 
In: Single-channel recording. 2nd ed. (Eds: Sakmann B, Neher E) 
Plenum Press, New York, pp. 397-482.

Colquhoun D, Hawkes AG (1995b) 
A Q-Matrix Cookbook. 
In: Single-channel recording. 2nd ed. (Eds: Sakmann B, Neher E) 
Plenum Press, New York, pp. 589-633.
"""

import numpy as np

def init_matrix(Q):
    """Updates the diagonal elements of Q in-place,
    i.e. Q will be changed by this function."""
    # Make sure that Q has square shape:
    if (not(Q.shape[0]==Q.shape[1])):
        print "Q doesn't have square shape in init_matrix(); aborting now."
        return
    for d in range(0,Q.shape[0]):
        Q[d,d]=0
        Q[d,d]=-np.sum(Q[d])

def p_inf(Q, debug=False):
    """Calculates p_inf for a Q matrix. Eq. 16, 17"""
    
    if (not(Q.shape[0]==Q.shape[1])):
        print "Q doesn't have square shape in init_matrix(); aborting now."
        return

    # add a unit vector:
    u = np.ones((Q.shape[0],1))
    S = np.concatenate((Q, u), axis=1)
    
    # Note that NumPy uses matrix multiplication for np.matrix,
    # but element-wise multiplication for np.array, so we have
    # to use np.dot here for matrix multiplication.
    p_mat = np.dot(np.transpose(u), np.linalg.inv((np.dot(S,np.transpose(S)))))
    return p_mat[0]

def mat_solve(Q, debug=False):
    """Returns the solutions to the equations given by
    Colquhoun & Hawkes, chapter 20, as a tuple:
    lambda, A. Eq. 88, 89."""

    lambda_r, X = np.linalg.eig(-Q)
    Y = np.linalg.inv(X)
    A = list()
    for i in range(0, Y.shape[0]):
        X_mat = np.empty((X.shape[0],1))
        Y_mat = np.empty((1,Y.shape[1]))
        X_mat[:,0] = X[:,i]
        Y_mat[0] = Y[i]
        A.append(np.dot(X_mat,Y_mat))

    # sort lambda and A. np.sort won't work here because
    # we need to synchronize lambda and A, so we use a
    # simple bubble sort.
    n_b = len(lambda_r)
    exchanged = True
    while exchanged and n_b >= 0:
        exchanged = False 
        for i in range(0, n_b-1):
            if (lambda_r[i] > lambda_r[i+1]):
                temp_l = lambda_r[i]
                temp_a = A[i]
                lambda_r[i] = lambda_r[i+1]
                A[i] = A[i+1]
                lambda_r[i+1] = temp_l
                A[i+1] = temp_a
                exchanged = True
        n_b -= 1

    if debug:
        for i in range(0, Y.shape[0]):
            print lambda_r[i]
            print A[i]
    return lambda_r, A

def p(t, p_0, p_inf, lambda_r, A):
    """Returns the probality of a channel being in a certain 
    state at time t. The state is characterised by p_inf (the
    probablity at equilibrium) and p_0 (the initial probability).
    The rates lambda_r and the amplitude terms A are the eigen-
    values and the eigenvectors of -Q, respectively. Eq. 24,26,27."""
    p_ret = np.empty((p_0.shape[0], len(t)))
    for j in range(0, p_0.shape[0]):
        sum_p = 0
        for i in range(1, len(A)):
            w_ij = 0
            for r in range(0, A[i].shape[1]):
                w_ij += p_0[r] * A[i][r,j]
            sum_p += w_ij*np.exp(-t*lambda_r[i])
 
        p_ret[j] = p_inf[j] + sum_p

    return p_ret
