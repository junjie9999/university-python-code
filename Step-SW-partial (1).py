from pylab import *
from scipy.integrate import quad
from numpy import linalg as LA

#Define parameters
L=2
Nbasis=120
k=1000
#DelV=1
#DelV=1/2*(k*(x)**2)
#Matrix elements of the kinetic energy operator
def KE(i,j):
    if i==j:
       ke = (1/2)*((i*pi/L)**2)#*((hbar)**2/2*m)# kinetic hamiltonian
    else:
       ke=0
    return ke

#Basis functions
def basis(i,x):
    if i % 2 == 0:
       fn = sin(i*pi*x/L)
    else:
       fn = cos(i*pi*x/L)
    return fn

#Potential energy
def steppot(x):
    DelV=1/2*(k*(x)**2)
    if abs(x) <= L/2:
       fn = DelV
    else:
       fn = 0
    return fn

#Integrand in calculation of matrix energy of potential energy 
def integrand(x,i,j):
    return basis(i,x)*steppot(x)*basis(j,x)

#Matrix elements of the potential energy
def V(i,j):
    #call scipy function for numerical integration
    result = quad(integrand, -L/2, L/2, args=(i,j))
    return result[0]

#Initialise Hamiltonian matrix with zeros
H=np.zeros((Nbasis,Nbasis))

#Assign Hamiltonian matrix elements
for i in range(0,Nbasis):
    for j in range(0,Nbasis):
        H[i][j] = KE(i+1,j+1)+V(i+1,j+1)

#TODO: Solve eigenvalue problem using numpy and sort in order
eigenValues, eigenVectors = LA.eig(H)
idx = eigenValues.argsort()
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:,idx]
#Print eigenvalues
for state in range(0,Nbasis):
 print('n=',state+1, 'E=',eigenValues[state])
#TODO: Print eigenvalues

xdata = []
eigfunc =[]

def xrange(start, end, step):
 while start <= end: 
  yield start
  start += step
  
for state in range(0,2): #print for the first five state 
  del xdata[:]
  del eigfunc[:]
  for x in xrange(-L/2, L/2, 0.01):
    xdata.append(x)
    psi=0
    for j in range(0,Nbasis):
        psi=psi+eigenVectors[j,state]*basis(j+1,x)
    eigfunc.append(psi)
    eig=state+1
    
  plt.plot(xdata,eigfunc,label='n= %i' %eig)
  
plt.xlabel('x (Bohr)')
plt.ylabel('psi')
plt.margins(x=0.1, y=0.1)
plt.title('Eigenfunctions')
plt.grid(True)
plt.savefig("eigenfunctions.png")
plt.legend(loc=4)
plt.show()
