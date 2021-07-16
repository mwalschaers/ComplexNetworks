#!users/jussieu/walschaers/.conda/envs/mattia_graph_env/bin python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 10:20:40 2018

Evaluating intensity correlations for photon-subtracted states

@author: Mattia Walschaers
"""

import igraph as graph
import sys
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import datetime
import math

def matrix_J(m):
    One = np.eye(m)
    Zero = np.zeros((m,m))
    return np.block([[Zero, -1.*One],[One, Zero]])


"""Evaluate correlations <a_i a_j>, <a_i^{\dag} a_j>, etc, based on the covariance matrix"""
def correlation_2(f1,f2):
    global V
    Jf1 = np.dot(matrix_J(int(np.size(f1[0])/2.)),f1[0])
    Jf2 = np.dot(matrix_J(int(np.size(f2[0])/2.)),f2[0])
    if f1[1]+f2[1]==0:
        return (np.dot(f1[0],np.dot(V,f2[0]))-np.dot(Jf1,np.dot(V,Jf2))+1j*(np.dot(f1[0],np.dot(V,Jf2))+np.dot(Jf1,np.dot(V,f2[0]))))/4.
    elif f1[1]+f2[1]==1:
        return (np.dot(f1[0],np.dot(V,f2[0]))+np.dot(Jf1,np.dot(V,Jf2))+1j*(np.dot(f1[0],np.dot(V,Jf2))-np.dot(Jf1,np.dot(V,f2[0]))))/4. - np.dot(f1[0],f2[0])/2. -1j*np.dot(f1[0],Jf2)/2.
    else:
        return (np.dot(f1[0],np.dot(V,f2[0]))-np.dot(Jf1,np.dot(V,Jf2))-1j*(np.dot(f1[0],np.dot(V,Jf2))+np.dot(Jf1,np.dot(V,f2[0]))))/4.

"""Creation of cluset state covariance matrix"""
def createCluster(V0,g):
    a = g.get_adjacency()
    A = np.array(a.data)
    m = g.vcount()
    One = np.eye(m)
    Zero = np.zeros((m,m))
    CZ = np.block([[One,Zero],[A,One]])
    V = np.dot(CZ,np.dot(V0,np.transpose(CZ)))
    return V

def createSqueezedVacuum(squeezing_dB):
    diag=np.power(10,0.1*squeezing_dB)
    return np.diag(np.concatenate((diag,np.power(diag,-1.))))

def addNoise(V, nu):
    l = int(np.size(V[0]))
    return V + np.diag(nu*np.random.random_sample((l)))


def factorial2(n):

    if n <= 0 or n == 1:
        return 1

    def multiply(coef):
        val = 1
        for i in range(n // 2):
            val *= (2 * i + coef)
        return val

    if n % 2 == 1:
        return n * multiply(1)

    return multiply(2)

##########################################
"""Group all different classes of terms for number-correlation matrix"""
##########################################
def const(x,n,m):
    nn=min(n,m)
    mm=max(n,m)
    if x == 0:
        return factorial2(n-1)*factorial2(m-1)
    else:
        return math.factorial(nn)*math.factorial(mm)/(math.factorial(nn-x)*math.factorial(mm-x)*math.factorial(x))*factorial2(n-x-1)*factorial2(m-x-1)


def singleModeCor(m,n,g):
    ada = correlation_2([g,1],[g,0])
    aa = correlation_2([g,0],[g,0])
    adad = correlation_2([g,1],[g,1])
    r = min(m,n)
    som=0

    if (m==0) and (n==0):
        return 1

    elif(m+n)% 2 != 0:
        return 0

    elif (m%2 == 0) and (n%2==0):
        for x in range(0,r+2,2):
           cx = const(x,n,m)
           som = som + cx*(ada**x)*(aa**((n-x)/2))*(adad**((m-x)/2))

    else:
        for x in range(1,r+2,2):
            cx = const(x,n,m)
            som = som + cx*(ada**x)*(aa**((n-x)/2))*(adad**((m-x)/2))

    return som


def class1(f1,f2,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    return ada1*a2a + ada2*a1a

def class21(f1,f2,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    return ada1*ada2

def class22(f1,f2,g):
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    return a1a*a2a

def class3(f1,f2,f3,f4,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    ada3 = correlation_2([g,1],f3)
    ada4 = correlation_2([g,1],f4)
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    a3a = correlation_2(f3,[g,0])
    a4a = correlation_2(f4,[g,0])
    return (ada1*ada2*a3a*a4a + ada1*a2a*ada3*a4a
            +a1a*ada2*ada3*a4a + ada1*a2a*a3a*ada4
            +a1a*ada2*a3a*ada4 + a1a*a2a*ada3*ada4)

def class41(f1,f2,f3,f4,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    ada3 = correlation_2([g,1],f3)
    ada4 = correlation_2([g,1],f4)
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    a3a = correlation_2(f3,[g,0])
    a4a = correlation_2(f4,[g,0])
    return (ada1*ada2*ada3*a4a + ada1*ada2*a3a*ada4
            + ada1*a2a*ada3*ada4 + a1a*ada2*ada3*ada4)

def class42(f1,f2,f3,f4,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    ada3 = correlation_2([g,1],f3)
    ada4 = correlation_2([g,1],f4)
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    a3a = correlation_2(f3,[g,0])
    a4a = correlation_2(f4,[g,0])
    return (a1a*a2a*a3a*ada4 + a1a*a2a*ada3*a4a
            +a1a*ada2*a3a*a4a + ada1*a2a*a3a*a4a)

def class51(f1,f2,f3,f4,g):
    ada1 = correlation_2([g,1],f1)
    ada2 = correlation_2([g,1],f2)
    ada3 = correlation_2([g,1],f3)
    ada4 = correlation_2([g,1],f4)
    return ada1*ada2*ada3*ada4

def class52(f1,f2,f3,f4,g):
    a1a = correlation_2(f1,[g,0])
    a2a = correlation_2(f2,[g,0])
    a3a = correlation_2(f3,[g,0])
    a4a = correlation_2(f4,[g,0])
    return (a1a*a2a*a3a*a4a)

def class0(f1,f2,f3,f4):
    f12 = correlation_2(f1,f2)
    f13 = correlation_2(f1,f3)
    f14 = correlation_2(f1,f4)
    f23 = correlation_2(f2,f3)
    f24 = correlation_2(f2,f4)
    f34 = correlation_2(f3,f4)
    return (f12*f34 + f13*f24 + f14*f23)

def singleModeSubCorrelator(f,g,n):
    f1d = f[0]
    f2d = f[1]
    f2 = f[2]
    f1 = f[3]
    A0 = singleModeCor(n,n,g)
    som = A0*class0(f1d,f2d,f2,f1)

    som1 = A0*correlation_2(f1d,f1)
    som2 = A0*correlation_2(f2d,f2)

    if n > 0 :
        A1 = singleModeCor(n-1,n-1,g)
        som = som + n**2 * A1 * (
                         class1(f1d,f2d,g)*correlation_2(f2,f1)
                        +class1(f2,f1,g)*correlation_2(f1d,f2d)
                        +class1(f1d,f2,g)*correlation_2(f2d,f1)
                        +class1(f2d,f1,g)*correlation_2(f1d,f2)
                        +class1(f1d,f1,g)*correlation_2(f2d,f2)
                        +class1(f2d,f2,g)*correlation_2(f1d,f1))

        som1 = som1 + n**2 * A1 * class1(f1d,f1,g)
        som2 = som2 + n**2 * A1 * class1(f2d,f2,g)

    if n > 1:
        A21 = singleModeCor(n-2,n,g)
        som = som + n*(n-1)*A21*(class21(f1d,f2d,g)*correlation_2(f1,f2)+class21(f1d,f2,g)*correlation_2(f2d,f1)
                        +class21(f1d,f1,g)*correlation_2(f2d,f2)+class21(f1,f2,g)*correlation_2(f1d,f2d)
                        +class21(f2d,f1,g)*correlation_2(f1d,f2)+class21(f2d,f2,g)*correlation_2(f1d,f1))
        som1 = som1 + n*(n-1)*A21*class21(f1d,f1,g)
        som2 = som2 + n*(n-1)*A21*class21(f2d,f2,g)
        A22 = singleModeCor(n,n-2,g)
        som = som + n*(n-1)*A22*(class22(f1d,f2d,g)*correlation_2(f1,f2)+class22(f1d,f2,g)*correlation_2(f2d,f1)
                        +class22(f1d,f1,g)*correlation_2(f2d,f2)+class22(f1,f2,g)*correlation_2(f1d,f2d)
                        +class22(f2d,f1,g)*correlation_2(f1d,f2)+class22(f2d,f2,g)*correlation_2(f1d,f1))
        som1 = som1 + n*(n-1)*A22*class22(f1d,f1,g)
        som2 = som2 + n*(n-1)*A22*class22(f2d,f2,g)

        A3 = singleModeCor(n-2,n-2,g)
        som = som + n**2 * (n-1)**2 * A3 * class3(f1d,f2d,f2,f1,g)

    if n > 2:
        A41 = singleModeCor(n-3,n-1,g)
        som = som + n**2 * (n-1) * (n-2) * A41 *class41(f1d,f2d,f2,f1,g)

        A42 = singleModeCor(n-1,n-3,g)
        som = som + n**2 * (n-1) * (n-2) * A42 *class42(f1d,f2d,f2,f1,g)

    if n > 3:
        A51 = singleModeCor(n-4,n,g)
        som = som + n * (n-1) * (n-2) * (n-3) * A51 *class51(f1d,f2d,f2,f1,g)

        A52 = singleModeCor(n,n-4,g)
        som = som + n * (n-1) * (n-2) * (n-3) * A52 *class52(f1d,f2d,f2,f1,g)


    return (som/A0 - som1*som2/(A0**2))


def singleModeSubCorrelatorDiag(f,g,n):
    f1d = f[0]
    f2d = f[1]
    f2 = f[2]
    f1 = f[3]
    A0 = singleModeCor(n,n,g)
    som = A0*class0(f1d,f2d,f2,f1)

    som1 = A0*correlation_2(f1d,f1)
    som2 = A0*correlation_2(f2d,f2)

    if n > 0 :
        A1 = singleModeCor(n-1,n-1,g)
        som = som + n**2 * A1 * (
                         class1(f1d,f2d,g)*correlation_2(f2,f1)
                        +class1(f2,f1,g)*correlation_2(f1d,f2d)
                        +class1(f1d,f2,g)*correlation_2(f2d,f1)
                        +class1(f2d,f1,g)*correlation_2(f1d,f2)
                        +class1(f1d,f1,g)*correlation_2(f2d,f2)
                        +class1(f2d,f2,g)*correlation_2(f1d,f1))

        som1 = som1 + n**2 * A1 * class1(f1d,f1,g)
        som2 = som2 + n**2 * A1 * class1(f2d,f2,g)

    if n > 1:
        A21 = singleModeCor(n-2,n,g)
        som = som + n*(n-1)*A21*(class21(f1d,f2d,g)*correlation_2(f1,f2)+class21(f1d,f2,g)*correlation_2(f2d,f1)
                        +class21(f1d,f1,g)*correlation_2(f2d,f2)+class21(f1,f2,g)*correlation_2(f1d,f2d)
                        +class21(f2d,f1,g)*correlation_2(f1d,f2)+class21(f2d,f2,g)*correlation_2(f1d,f1))
        som1 = som1 + n*(n-1)*A21*class21(f1d,f1,g)
        som2 = som2 + n*(n-1)*A21*class21(f2d,f2,g)
        A22 = singleModeCor(n,n-2,g)
        som = som + n*(n-1)*A22*(class22(f1d,f2d,g)*correlation_2(f1,f2)+class22(f1d,f2,g)*correlation_2(f2d,f1)
                        +class22(f1d,f1,g)*correlation_2(f2d,f2)+class22(f1,f2,g)*correlation_2(f1d,f2d)
                        +class22(f2d,f1,g)*correlation_2(f1d,f2)+class22(f2d,f2,g)*correlation_2(f1d,f1))
        som1 = som1 + n*(n-1)*A22*class22(f1d,f1,g)
        som2 = som2 + n*(n-1)*A22*class22(f2d,f2,g)

        A3 = singleModeCor(n-2,n-2,g)
        som = som + n**2 * (n-1)**2 * A3 * class3(f1d,f2d,f2,f1,g)

    if n > 2:
        A41 = singleModeCor(n-3,n-1,g)
        som = som + n**2 * (n-1) * (n-2) * A41 *class41(f1d,f2d,f2,f1,g)

        A42 = singleModeCor(n-1,n-3,g)
        som = som + n**2 * (n-1) * (n-2) * A42 *class42(f1d,f2d,f2,f1,g)

    if n > 3:
        A51 = singleModeCor(n-4,n,g)
        som = som + n * (n-1) * (n-2) * (n-3) * A51 *class51(f1d,f2d,f2,f1,g)

        A52 = singleModeCor(n,n-4,g)
        som = som + n * (n-1) * (n-2) * (n-3) * A52 *class52(f1d,f2d,f2,f1,g)


    return (som/A0 - som1*som2/(A0**2) + som1/A0)

###########################################

##########################################
""" Make number-correlation matrix """
##########################################
def makeCorrelationMatrixSingleModeSub(n,g):
    global V
    m = np.size(V[0])
    One = np.eye(m)
    Cmatrix= np.empty([int(m/2.),int(m/2.)])


    for i in range(int(m/2)):
        for j in range(int(i+1)):
            f = [[One[i],1],[One[j],1],[One[j],0],[One[i],0]]
            if i!=j :
                cx = singleModeSubCorrelator(f,g,n)
                Cmatrix[i,j] = np.real(cx)
                Cmatrix[j,i] = np.real(cx)
            else:
                cx = singleModeSubCorrelatorDiag(f,g,n)
                Cmatrix[i,i] = np.real(cx)

    return Cmatrix



##########################################
""" Make quadrature covariance matrix after subtraction """
##########################################
def getCovarianceMatrix(n,g):

    global V
    m = np.size(V[0])
    One = np.eye(m)
    qq = np.empty([m,m])

    for i in range(int(m/2)):
        for j in range(int(i+1)):
            A0 = singleModeCor(n,n,g)
            aidajd = A0*correlation_2([One[i],1],[One[j],1])
            aiaj = A0*correlation_2([One[i],0],[One[j],0])
            aidaj = A0*correlation_2([One[i],1],[One[j],0])
            ajdai = A0*correlation_2([One[j],1],[One[i],0])
            if n > 0:
                A1 = singleModeCor(n-1,n-1,g)
                aidajd = aidajd + n**2 * A1 * class1([One[i],1],[One[j],1],g)
                aiaj = aiaj + n**2 * A1 * class1([One[i],0],[One[j],0],g)
                aidaj = aidaj + n**2 * A1 * class1([One[i],1],[One[j],0],g)
                ajdai = ajdai + n**2 * A1 * class1([One[j],1],[One[i],0],g)
            if n > 1:
                A21 = singleModeCor(n-2,n,g)
                aidajd = aidajd + n*(n-1)*A21*class21([One[i],1],[One[j],1],g)
                aiaj = aiaj + n*(n-1)*A21*class21([One[i],0],[One[j],0],g)
                aidaj = aidaj + n*(n-1)*A21*class21([One[i],1],[One[j],0],g)
                ajdai = ajdai + n*(n-1)*A21*class21([One[j],1],[One[i],0],g)

                A22 = singleModeCor(n,n-2,g)
                aidajd = aidajd + n*(n-1)*A22*class22([One[i],1],[One[j],1],g)
                aiaj = aiaj + n*(n-1)*A22*class22([One[i],0],[One[j],0],g)
                aidaj = aidaj + n*(n-1)*A22*class22([One[i],1],[One[j],0],g)
                ajdai = ajdai + n*(n-1)*A22*class22([One[j],1],[One[i],0],g)

            qq[i,j] = np.real((aidaj + ajdai + (aidajd + aiaj)))
            qq[i+int(m/2),j+int(m/2)] = np.real(aidaj + ajdai - (aidajd + aiaj))
            qq[i,j+int(m/2)] = np.real(1j*(aidajd - aiaj - aidaj + ajdai))
            qq[i+int(m/2),j] = np.real(1j*(aidajd - aiaj + aidaj - ajdai))
            if i != j:
                qq[j,i] = np.real(aidaj + ajdai + (aidajd + aiaj))
                qq[j+int(m/2),i+int(m/2)] = np.real(aidaj + ajdai - (aidajd + aiaj))
                qq[j,i+int(m/2)] = np.real(1j*( aidajd - aiaj + aidaj - ajdai))
                qq[j+int(m/2),i] = np.real(1j*( aidajd - aiaj - aidaj + ajdai))

    return (qq)/np.real(A0)+One



#############################################################################
""" Using all of the above """
#############################################################################
ID=int(sys.argv[1])
n=int(sys.argv[2])
m=int(sys.argv[3])
SQ = float(sys.argv[4])
repetitions = int(sys.argv[5])
rand = int(sys.argv[6])

filename = 'READ_ME_ID'+ str(ID) + '.txt'
f = open(filename, 'w')

f.write('This file was created at ' + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
f.write('\n')
f.write('We created a Barabasi-Albert graph state, with the following properties:')
f.write('Size of the network: ' + str(n) +'\n')
f.write('Number connection for every newly added vertex: ' + str(m)+'\n')
f.write('Squeezing: ' + str(SQ)+'\n')
f.write('Number of different networks generated: ' + str(repetitions)+'\n')
f.write('The photons were subtracted in a random vertex (0=no, 1=yes)? ' + str(rand)+'\n')
f.write('\n')



"""Generate all the different realisations and subtract the relevant number of photons"""
for rep in range(repetitions):

    connected = False
    while not connected:
        """Here one has to change the graph type to generate Watts-Stragatz etc"""
        g = graph.Graph.Barabasi(n=n,m=m)
        connected = g.is_connected()

    squeeze=SQ*np.ones(n)
    V0 = createSqueezedVacuum(squeeze)

    ad = g.get_adjacency()
    adj = np.array(ad.data)

    filenameA = 'AdjacencyMatrix_Realisation'+ str(rep) + '_ID' + str(ID) + '.txt'
    np.savetxt(filenameA, adj, delimiter=',')

    deGraph = np.sum(adj, axis=0).tolist()
    maxDegVert = np.argmax(deGraph)

    V = createCluster(V0,g)
    M = np.size(V[0])
    One = np.eye(int(M/2))
    Zero = np.zeros([int(M/2),int(M/2)])
    if rand == 1:
        randomVertex = np.random.randint(n)
        subtractionMode = np.array(One[randomVertex].tolist()+Zero[0].tolist())
    else:
        subtractionMode = np.array(One[maxDegVert].tolist()+Zero[0].tolist())


    for nSubtractions in range(11):


        correlationMatrix = makeCorrelationMatrixSingleModeSub(nSubtractions,subtractionMode)
        covarianceMatrix = getCovarianceMatrix(nSubtractions,subtractionMode)
        filenameC = 'CorrelationMatrix_sub' + str(nSubtractions) + '_Realisation'+ str(rep) + '_ID' + str(ID) + '.csv'
        np.savetxt(filenameC, correlationMatrix, delimiter=',')
        filenameV = 'CovarianceMatrix_sub' + str(nSubtractions) + '_Realisation'+ str(rep) + '_ID' + str(ID) + '.csv'
        np.savetxt(filenameV, covarianceMatrix, delimiter=',')


f.write('The execution of the code finished at: ' + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
f.close()
