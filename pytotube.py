# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 09:54:35 2013

@author: ahmed
PytoTube : Elliptic 1d Steady State Heat Equation Solver

"""

from numpy import zeros, array, nan, sqrt
from numpy.linalg.linalg import solve
from numpy.lib.shape_base import tile
from numpy.core.fromnumeric import shape

from scimath.units.length import feet, meter, cm
from scimath.units.api import has_units, UnitArray

from pylab import plot, subplot, grid, legend, title


class post1d():

    @has_units
    def ft2m(self, a):
        """ Add two arrays in ft and convert them to m.

        """
        return a * feet/m

class Laplace():
    """
    1d Laplace equation solver

    Parameters
    ----------
    L : Float :
      The length of rod.
    N : Integer
      Number of elements (Uniform)

    Returns
    -------

    Results
    -------

    Examples
    --------
    >>> MyProblem = Laplace(N = 5, L = 1.5 )

    """
    # --- DEFAULT METHODS ---
    def __del__(self):
        pass

    def __init__(self, L, N):
        self.deltaX = float(L / float(N))
        self.N = int(N + 1)
        self.M = self.N
        self.L = float(L)
        self.X = []
        # TODO !!! move this part to __call__ method

        self.A = array(zeros([self.M, self.M]))



        #self.filtera = array(self.A +1, dtype=bool)
        self.filtera = self.A == self.A         # COOL ! :)
        self.DOF = array(tile(True, (self.M)), dtype=bool)
        for i in range(self.M):         # assuming DX=Constant
            self.X.append(round(self.deltaX*i, 15))   # To avoid  stupid 0.199

        self.B = array(zeros([self.M, 1]))
        self.T = array(tile(nan, (self.M,)))

        # TODO for solve method
        self.bcConduction = array(tile(False, (self.M, )), dtype=bool)
        self.bcConvection = array(tile(False, (self.M, )), dtype=bool)
        self.bcT = array(tile(True, (self.M, )), dtype=bool)
        #self.bcRadiation = array(tile(False, (self.M, 1)), dtype=bool)

        internalNodes = array([1., -2., 1.])
#
        for i in range(1, int(self.M-1)):
#            self.A[i, range(i-1, i+2)] = internalNodes
            self.B[i] *= self.deltaX**2
# %------------------------------------------------

    def __call__(self):
#        self.print_matrice(self.A, 'A')
#        self.print_matrice(self.B, 'B')
#        self.print_matrice(self.bcConduction, 'Conduction')
#        self.print_matrice(self.bcConvection, 'Convection')
#        self.print_matrice(self.bcT, 'Temperature')
        # Construct matrix A and B according to self.bcXXX matrices
        pass

    # --- USER METHODS ---
    def solve(self):

        # If Convection NOT Implemented
        # If Conduction (Neumann) [T(i) - T(i-1)] /DX B *= DX
        # If Temperature
        # TODO IMPORTANT
        # DO NOT GIVE 1-self.M DIRECTLY [1 -2 1]
        # Instead check If Convection / Conduction if internal nodes are BC
        # Now only
        internalNodes = array([1., -2., 1.])

        for i in range(1, int(self.M-1)):
            self.A[i, range(i-1, i+2)] = internalNodes
            self.B[i] *= self.deltaX**2
        # %------------------------------------------------


        A = self.A[self.filtera]
        a = int(sqrt(A.shape[0]))
        A.shape = ([a, a])
        B = self.B[self.bcT]

        self.T[self.DOF] = solve(A, B)


    def plot(self):
        self.ax = subplot(1,1,1)
        self.label += "$\Delta_{x}$=" + str(round(self.deltaX,3)) + "m"
        self.p, = self.ax.plot(self.X, self.T, '-o',label=self.label)
        grid(True)
        handles, labels = self.ax.get_legend_handles_labels()
        # reverse the order
        self.ax.legend(handles[::-1], labels[::-1])
        # or sort them by labels
        import operator
        hl = sorted(zip(handles, labels),
        key=operator.itemgetter(1))
        handles2, labels2 = zip(*hl)
        self.ax.legend(handles2, labels2,bbox_to_anchor=(0., 1), loc=2, borderaxespad=0.)
        pass

    def print_matrice(self, val, label, n=13):
        print u'\x11' * n + ' ' + str(label) + ' ' + '\x10'*n
        for i in range(shape(val)[0]):
            print str(i) + u' \x10 ' + str(val[i, :])
        print u'\x16' * 30
        print

    def x2index(self, x):
        index = self.X.index(x)
        return index

    # --- LAPLACE BOUNDARY CONDITIONS ---
    def dirichlet(self, index, Scalar):

        self.B += -Scalar * self.A[:, index:index+1].reshape(self.B.shape)

#        self.B += Scalar * vstack(self.A[:, index:index+1])
        # TODO Check if the index is compatible with multiple values (lists)
        # Remove the index of known(s) scalar(s) from B
        self.DOF[index] = False

        # Divide A and B
        # TODO !!!
        # Sum the row of A with i5ndex-1 row
        self.A[index-1] += self.A[index]
        self.B[index-1] += self.B[index]

    def neumann(self, index, Value):
        """
        dT/dx = Value
        """
        self.B[index] = self.deltaX * Value
        if index == self.M - 1:
            self.A[index, self.M-3:] = array([0, -1, 1])

        else:
            self.A[index, index:index+3] = array([1, -1, 0])
            # Fixed for -/+ wrong direction of heat flow


class Heat1d(Laplace,post1d):
    # --- HEAT 1D CLASS ---
    def __del__(self):
        pass

# TODO FIX THIS Heat1d.__init__()  problem
#    def __init__(self):
#        print 'hellllll'

# TODO !!! Support list with Bcs,
# example self.dirichlet(x = [0, 0.1, 0.2] , T = [100,200,300])
    # --- HEAT BOUNDARY CONDITIONS ---

    def Conduction(self, x, k):
        """
        dT/dx x | X  =  K1
        """

        # Call Neumann and assign to self.bcConduction[index] => True
        index = self.x2index(x)
        self.bcConduction[index] = True

        self.neumann(index, k)

    def Convection(self, x, k, h, Ta):

        """
        dT/dx = -h/k * (T_i - T_ambient)
        h = W/(m^2K)
        """
        # TODO !!! NOT Implemented
        # Call Neuman with two coefficent [C1 C0] C1*Ti + C0
        # and assign to self.bcConvection[index] => True
        index = self.x2index(x)
        self.label = "h=" + str(h) + " $W /(m^2 \cdot K) \;$"
        self.bcConvection[index] = True

        self.B[index] = -h * self.deltaX * Ta
        if index == self.M - 1:
            self.A[index, self.M-3:] = array([0, k, -k-h*self.deltaX])

        else:
            self.A[index, index:index+3] = array([-k-h*self.deltaX, k, 0])
            # Fixed for -/+ wrong direction of heat flow

    def Temperature(self, x, T):
        """
        T(x) = T
        """
        index = self.x2index(x)
        self.dirichlet(index, T)
        self.filtera[:, index] = False
        self.filtera[index, :] = False
        self.bcT[index] = False
        self.T[index] = T

        pass

    def source(self, k, Q, index=0):
        """
        Q = \dot{q}/k
        dot_q = W/m^3
        k = W/(mK)
        """
        self.B[1:self.M-1] += -Q/k
        pass



if __name__ == '__main__':
#    # TODO !!! There is a problem with internal node BCs
#    # Because __init__ gives [1 -2 1] every internal node
#    # At now, only x= 0 and x=L Bcs are supported
#    # How to give internal nodes as Bs
#    # Give all BCs, except self.filtera rows, give [1 -2 1]

    mesh = [10]
    for h  in [0, 10, 20, 50,100]:
        for i in mesh:
            w = Heat1d(1.0, i)
            w.source(75, 1e5)
            w.Conduction(0.0, -100.)
            w.Temperature(.0, 400)
    #        w.Temperature(0.0, 200)
            w.Convection(1.0, 75, h, 200)
            w.solve()
            w.plot()
            title('For N : ' + str(mesh) )
def heat():
    cla();clf();
    mesh = [10]
    for h  in [0, 10, 20, 50,100]:
        for i in mesh:
            w = Heat1d(L=1.0, N=i)
            w.source(k=75, Q=1e5)
            w.Conduction(x=0.0, k= -10.)
            w.Temperature(x=0.0, T=700.)
            w.Convection(x=1.0, k=75, h=h, Ta=200.)
            w.solve()
            w.plot()
            title('For N : ' + str(mesh) )
