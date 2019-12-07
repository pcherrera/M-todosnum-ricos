import numpy as np
def phi(i,x):

    if i == 0:
        return (1-x[0]-x[1])*(1-2*x[0]-2*x[1])
    elif i == 1:
        return x[0]*(2*x[0]-1)
    elif i == 2:
        return x[1]*(2*x[1]-1)
    elif i == 3:
        return 4*x[0]*x[1]
    elif i == 4:
        return 4*x[1]*(1-x[0]-x[1])
    elif i == 5:
        return 4*x[0]*(1-x[0]-x[1])

def F_inv(Tri,x):
    p0 = Tri[0,:]
    p1 = Tri[1,:]
    p2 = Tri[2,:]
    BK = np.array([p1-p0,p2-p0]).T
    inv = np.linalg.inv(BK)
    return inv@(x-p0)

def FK(Points,x):
    N = len(x)
    p0 = Points[0,:]
    p1 = Points[1,:]
    p2 = Points[2,:]

    v = np.zeros([N,2])
    for i in range(0,N):
        v[i,:]= (np.array([p1-p0,p2-p0]).T)@x[i,:]+p0
    return v


def f(x):
    return x[0]**6

def int_FK(Tri,k):
    p0 = Tri[0,:]
    p1 = Tri[1,:]
    p2 = Tri[2,:]
    BK = np.array([p1-p0,p2-p0]).T
    det = np.linalg.det(BK)
    x = np.array([[0.24928674517091 ,   0.24928674517091],
                   [0.24928674517091 ,   0.50142650965818],
                   [0.50142650965818 ,   0.24928674517091],
                   [0.06308901449150 ,   0.06308901449150],
                   [0.06308901449150 ,   0.87382197101700],
                   [0.87382197101700 ,   0.06308901449150],
                   [0.31035245103378 ,   0.63650249912140],
                   [0.63650249912140 ,   0.05314504984482],
                   [0.05314504984482 ,   0.31035245103378],
                   [0.63650249912140 ,   0.31035245103378],
                   [0.31035245103378 ,   0.05314504984482],
                   [0.05314504984482 ,   0.63650249912140]])

    w = np.array([0.11678627572638,
                   0.11678627572638,
                   0.11678627572638,
                   0.05084490637021,
                   0.05084490637021,
                   0.05084490637021, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837])

    v = FK(Tri,x)
    r = 0
    for i in range(len(w)):
        r+= w[i]*phi(k,v[i])*np.cos(v[i][0])

    return abs(det)*0.5*r

ref = np.array([[0,0],[0,1],[1,0]])


