import numpy as np
import numpy.linalg as LA

#Initialise a 3x3 matrix with random matrix elements

#Solve eigenvalue problem
A = np.random.random((3,3))
eigenValues, eigenVectors = LA.eig(A)

#Sort eigenvalues and eigenvectors using the numpy.argsort function
idx = eigenValues.argsort()   
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:,idx]

#print eigenvalues
print(eigenValues)
