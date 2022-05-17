import numpy as np


def computeQ(n,r,t1,t2):
    """
    This function is to compute Q matrix in the QP problem.

    Parameters
    ----------
    n : int
        Polynormial order.
    r : int
        derivertive order:{1-4}
        1:minimum vel 
        2:minimum acc 
        3:minimum jerk 
        4:minimum snap
    t1 : float 
        start timestamp for polynormial.
    t2 : float
        end timestamp for polynormial.

    Returns
    -------
    Q : (n+1)*(n+1) np.array 
        Q matrix in solve_qp function.

    """
    shape = (n-r) * 2 + 1
    T = np.zeros(shape) 
    for i in range(shape):
        T[i] = t2**(i+1)-t1**(i+1)
    
    Q = np.zeros((n+1,n+1))
    
    for i in range(r,n+1):
        # print('i=',i)
        for j in range(i,n+1):
            # print('j=',j)
            k1 = i - r 
            k2 = j - r 
            k = k1 + k2 + 1 
            a = np.prod(np.arange(k1+1, k1+r+1, 1))
            b = np.prod(np.arange(k2+1, k2+r+1, 1))    
            Q[i,j] = a*b/k*T[k-1]
            Q[j,i] = Q[i,j]
            
    return Q

def calc_tvec(t,n_order,r):
    """
    This function is to compute pva constraints.

    Parameters
    ----------
    t : float
        Initial time for different section
    n_order : int
        Polynormial order.
    r : int
        derivertive order:{1-4}
        1:minimum vel 
        2:minimum acc 
        3:minimum jerk 
        4:minimum snap

    Returns
    -------
    tvec : 1 * (n_order+1) np.array
        pva constraints.

    """
    tvec = np.zeros((1,n_order+1))

    for i in range(r,n_order+1):
        #print("i:",i)
        a = np.prod(np.arange(i+1-r, i+1, 1))
        #print("a=",a)
        b = i - r
        #print("b=",b)
        tvec[0,i] = a * (t**b)
    return tvec

def re_shape_vaj(vaj):
    """

    Parameters
    ----------
    vaj : 1 * n np.array

    Returns
    -------
    vaj : (n,) np.array
        DESCRIPTION.

    """
    vaj = np.reshape(vaj,(len(vaj[0]),),order='F') 
    return vaj

def poly_val(poly,t,r):
    val = 0
    n = len(poly)-1   
    
    if r <= 0:
        for i in range(n+1):
            val = val + poly[i] * t**i 
    else:
        for i in range(r,n+1):
            # print(i)
            a = poly[i] * np.prod(np.arange(i-r+1, i+1, 1)) * t**(i-r)
            val = val + a 
            # print(val)
    return val
