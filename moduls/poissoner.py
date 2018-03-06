import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import splu,factorized

def laplace_splu(shape):
    xbins, ybins = shape # getting the shape of the global calculation domain
    # definig the diagonals of 2D Laplace matrix operator
    diag0 = (-4.0)*np.ones((xbins*ybins),dtype = 'd')	# central diagonal
    diag1 = np.ones((xbins*ybins),dtype = 'd')
    diag1[(ybins-1)::ybins]=0.0			#  Dirichlet boundary condition
    diag2 = np.ones((xbins*ybins),dtype = 'd')
    diag2[:-(ybins):ybins]=0.0			#  Dirichlet boundary condition
    diag3 = np.ones((xbins*ybins),dtype = 'd')
    diag3[0:(ybins-1)]=0.0				#  Dirichlet boundary condition
    diag4 = np.ones((xbins*ybins),dtype = 'd')
    diag4[-(ybins):]=0.0				#  Dirichlet boundary condition
    data = np.array([diag0,diag1,diag2,diag3,diag4])
    diags = np.array([0,1,-1, ybins, -ybins])

    A = spdiags(data, diags, xbins*ybins,xbins*ybins,format='csc')  # sparse Laplace matrix
    sp = splu(A)
    return sp.solve # factorization of inversed Laplace matrix (long)

#   sp = splu(A, permc_spec='NATURAL', options=dict(Equil=False, \
#   SymmetricMode = True, IterRefine=None))
#       return factorized(A) #factorization of inversed Laplace matrix (long)


def poiss_splu_solve(solve,b, dx, dy, charge = -1.):
    shape = b.shape		# saving the shape of RHS (shape of global calculation domain)
    potent = solve((-1.)*charge*b.flatten()*dx*dy)	# solve Poisson equation with factorized operator (fast)
    potent.shape = shape				# reshaping the flattened solution
    return potent
