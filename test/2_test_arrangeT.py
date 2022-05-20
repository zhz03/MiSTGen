
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 12:21:24 2022

@author: zhaoliang
"""

import numpy as np
from scipy.linalg import block_diag
# from qpsolvers import solve_qp
# from quadprog import solve_qp
import cvxopt
import matplotlib.pyplot as plt


def quadprog(H, f, L=None, k=None, Aeq=None, beq=None, lb=None, ub=None):
    """
    Input: Numpy arrays, the format follows MATLAB quadprog function: https://www.mathworks.com/help/optim/ug/quadprog.html
    Output: Numpy array of the solution
    """
    n_var = H.shape[1]

    P = cvxopt.matrix(H, tc='d')
    q = cvxopt.matrix(f, tc='d')

    if L is not None or k is not None:
        assert(k is not None and L is not None)
        if lb is not None:
            L = np.vstack([L, -np.eye(n_var)])
            k = np.vstack([k, -lb])

        if ub is not None:
            L = np.vstack([L, np.eye(n_var)])
            k = np.vstack([k, ub])

        L = cvxopt.matrix(L, tc='d')
        k = cvxopt.matrix(k, tc='d')

    if Aeq is not None or beq is not None:
        assert(Aeq is not None and beq is not None)
        Aeq = cvxopt.matrix(Aeq, tc='d')
        beq = cvxopt.matrix(beq, tc='d')

    sol = cvxopt.solvers.qp(P, q, L, k, Aeq, beq)

    return np.array(sol['x'])

def arrangeT(waypts,T):
    """
    This function is to calculate the initial 
        time series for each waypoint section.

    Parameters
    ----------
    waypts : 2d np.array
        2d way points.
    T : float
        Total execution time.

    Returns
    -------
    ts : np.array[float]
        Time to reach each waypoint.

    """
    x = waypts[:,1:] - waypts[:,0:-1] 
    
    dist = sum(x**2)**0.5
    k = T/sum(dist)
    new = np.cumsum(dist*k)
    ts = np.insert(new,0,0)
    
    return ts

def init_T(waypts,ave_speed):
    """
    This functino is to initialize total estimated T
    The idea in this function is to assign T by 
        using total path length / average speed 

    Parameters
    ----------
    waypts :  2 * n np.array
        A series of 2d waypoints.
    ave_speed : float
        A estimated average speed.

    Returns
    -------
    T : float 
        Initial total estimated time.

    Example
    -------
    waypts = np.array([[1,2,3],[1,2,3]])
    ave_speed = 2 # m/s
    
    """
    x = waypts[:,1:] - waypts[:,0:-1] 
    dist = sum(x**2)**0.5
    total_dist = sum(dist)
    T = total_dist/ave_speed
    return T

def computeQ(n,r,t1,t2):
    """
    

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

def minimum_snap_single(waypts,ts,n_order,n_deri,v0,a0,ve,ae):
    """
    This function is to calculate the soluation of the QP program
        given the velocity, acceleration constraints

    Parameters
    ----------
    waypts : TYPE
        DESCRIPTION.
    ts : TYPE
        DESCRIPTION.
    n_order : int
        Polynormial order.
    n_deri : int
        derivertive order:{1-4}
        1:minimum vel 
        2:minimum acc 
        3:minimum jerk 
        4:minimum snap
    v0 : float 
        Initial set velocity.
    a0 : float
        Initial set acceleration.
    ve : float
        End set velocity.
    ae : 
        End set acceleration.

    Raises
    ------
    Exception
        Q_all matrix must be a n*n  matrix

    Returns
    -------
    polys : 6 * (len(waypts)-1) np.array(float)
        Solution of the qp problem during the waypoint 
            and conforming to the 6 constraints.

    """
    p0 = waypts[0]
    pe = waypts[-1]
    n_poly = len(waypts)-1
    n_coef = n_order+1
    
    # n_deri = 3
    # compute Q
    Q_all = np.empty(0)
    
    for i in range(n_poly):
        # print(i)
        q_new = computeQ(n_order,n_deri,ts[i],ts[i+1])
        # print(len(q_new))
        Q_all = block_diag(Q_all,q_new)
    Q_all = np.delete(Q_all,0,axis=0)

    if len(Q_all) != len(Q_all[0]):
        raise Exception('Q_all is not a square matrix')
        
    b_all = np.zeros((len(Q_all),1))

    Aeq = np.zeros((4*n_poly+2,n_coef*n_poly))
    beq = np.zeros((4*n_poly+2,1))

    # start/terminal pva constraints  (6 equations)
    A1 = calc_tvec(ts[0],n_order,0)
    A2 = calc_tvec(ts[0],n_order,1)
    A3 = calc_tvec(ts[0],n_order,2)
    Aeq[0,0:n_coef] = A1
    Aeq[1,0:n_coef] = A2
    Aeq[2,0:n_coef] = A3
    
    Ae1 = calc_tvec(ts[-1],n_order,0)
    Ae2 = calc_tvec(ts[-1],n_order,1)
    Ae3 = calc_tvec(ts[-1],n_order,2)
    
    start_i = n_coef*(n_poly-1)
    end_i = n_coef*n_poly
    Aeq[3,start_i:end_i] = Ae1
    Aeq[4,start_i:end_i] = Ae2
    Aeq[5,start_i:end_i] = Ae3
    

    pva = [p0,v0,a0,pe,ve,ae]
    beq[0:6,0] = pva

    # mid p constraints    (n_ploy-1 equations)
    neq = 6
    for i in range(n_poly-1):
        # print(i)
        neq = neq + 1
        col1 = n_coef*(i+1)+1
        col2 = n_coef*(i+1+1)
        #print('col1',col1)
        #print('col2',col2)
        Aeq[neq-1,col1-1:col2] = calc_tvec(ts[i+1],n_order,0)
        beq[neq-1] = waypts[i+1]

    # continuous constraints  ((n_poly-1)*3 equations)
    
    for i in range(n_poly-1):
        tvec_p = calc_tvec(ts[i+1],n_order,0)
        tvec_v = calc_tvec(ts[i+1],n_order,1)
        tvec_a = calc_tvec(ts[i+1],n_order,2)
        # print('neq:',neq)
        neq = neq + 1
        st = n_coef*(i)
        end = n_coef*(i+1+1)
        # print('st:',st)
        # print("end",end)
        Aeq[neq-1,st:end] = np.concatenate((tvec_p,-tvec_p),axis=1)
        neq = neq + 1
        Aeq[neq-1,st:end] = np.concatenate((tvec_v,-tvec_v),axis=1)
        neq = neq + 1
        Aeq[neq-1,st:end] = np.concatenate((tvec_a,-tvec_a),axis=1)
        

    Aieq = None
    bieq = None
    x = quadprog(Q_all,b_all,Aieq,bieq,Aeq,beq)
    
    polys = np.reshape(x,(n_coef,n_poly),order='F')  
    return polys 

def polys_vals(polys,ts,tt,r):
    idx = 1
    N = len(tt)
    vals = np.zeros((1,N))
    
    for i in range(N):
        # print(i)
        t = tt[i]
        if t < ts[idx-1]:
            vals[0][i] = 0
        else: 
            while idx < len(ts) and t > ts[idx] + 0.0001:
                idx = idx + 1
                # print(idx)
            a = polys[:,idx-1]
            vals[0][i] = poly_val(a,t,r)
    
    return vals 

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

def minimum_snap_2d(waypts,ts,n_order,n_deri,v0,a0,ve,ae):
    """
    This function is to output polys_x and polys_y given 2d waypoints

    Parameters
    ----------
    waypts : 2 * n np.array
        A series of 2d waypoints.
    ts : np.array(float)
        A series of timestamp to reach to each waypoint.
    n_order : int
        Polynormial order.
    n_deri : int
        derivertive order:{1-4}
        1:minimum vel 
        2:minimum acc 
        3:minimum jerk 
        4:minimum snap
    v0 : 1 * 2 np.array
        Velocity constrains at initial waypoint.
    a0 : 1 * 2 np.array
        Acceleration constrains at initial waypoint.
    ve : 1 * 2 np.array
        Velocity constrains at end waypoint.
    ae : 1 * 2 np.array
        Acceleration constrains at end waypoint.

    Returns
    -------
    polys_x,polys_y : 6 * (len(waypts)-1) np.array(float)
        Solution of the qp problem during the waypoint 
            and conforming to the 6 constraints.

    """
    polys_x= minimum_snap_single(waypts[0],ts,n_order,n_deri,v0[0],a0[0],ve[0],ae[0])
    polys_y= minimum_snap_single(waypts[1],ts,n_order,n_deri,v0[1],a0[1],ve[1],ae[1])
    return polys_x,polys_y

def mist_2d_gen(waypts,v0,a0,ve,ae,T=None,ave_speed=None,interval=0.01,n_order=5,n_deri=3):
    """
    This function is to generate minimum snap trajectory 
        and its corresponding timestamp

    Parameters
    ----------
    waypts : TYPE
        DESCRIPTION.
    v0 : TYPE
        DESCRIPTION.
    a0 : TYPE
        DESCRIPTION.
    ve : TYPE
        DESCRIPTION.
    ae : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    ave_speed : TYPE
        DESCRIPTION.
    interval : TYPE
        DESCRIPTION.
    n_order : TYPE
        DESCRIPTION.
    n_deri : TYPE
        DESCRIPTION.

    Returns
    -------
    xxs,yys : n * 1 np.array[float]
        X,Y coordinate of generated minimum snap trajectory(MiST).
    tts : n * 1 np.array[float]
        Corresponding timestamp
    polys_x,polys_y : 6 * (len(waypts)-1) np.array(float)
        Solution of the qp problem during the waypoint 
            and conforming to the 6 constraints.

    """
    if T == None and ave_speed != None:
        T = init_T(waypts,ave_speed)
    
    ts = arrangeT(waypts,T)
    
    polys_x,polys_y = minimum_snap_2d(waypts,ts,n_order,n_deri,v0,a0,ve,ae)
    
    xxs = np.array([])
    yys = np.array([])
    tts = np.array([])
    
    for i in range(len(polys_x[0])):
        tt = np.arange(ts[i],ts[i+1],interval)
        xx = polys_vals(polys_x, ts, tt, 0)
        yy = polys_vals(polys_y, ts, tt, 0)
        # plt.plot(xx[0],yy[0])
        xxs = np.append(xxs,xx)
        yys = np.append(yys,yy)
        tts = np.append(tts,tt)
        
    return xxs,yys,tts,polys_x,polys_y,ts

def mist_2d_vaj_gen(xxs,yys,tts,polys_x,polys_y,ts):
    """
        This function is to generate corresponding MiST
            velocity, acceleration and jerk

    Parameters
    ----------
    xxs,yys : n * 1 np.array[float]
        X,Y coordinate of generated minimum snap trajectory(MiST).
    tts : n * 1 np.array[float]
        Corresponding timestamp.
    polys_x,polys_y : 6 * (len(waypts)-1) np.array(float)
        Solution of the qp problem during the waypoint 
            and conforming to the 6 constraints.
    ts : np.array[float]
        Time to reach each waypoint.

    Returns
    -------
    vxx,vyy : n * 1 np.array[float]
        Corresponding MiST velocity.
    axx,ayy : n * 1 np.array[float]
        Corresponding MiST accleration.
    jxx,jyy : n * 1 np.array[float]
        Corresponding MiST jerk.

    """    
    tt = tts
    vxx = polys_vals(polys_x,ts,tt,1)
    axx = polys_vals(polys_x,ts,tt,2)
    jxx = polys_vals(polys_x,ts,tt,3)
    vyy = polys_vals(polys_y,ts,tt,1)
    ayy = polys_vals(polys_y,ts,tt,2)
    jyy = polys_vals(polys_y,ts,tt,3)

    vxx = re_shape_vaj(vxx)
    axx = re_shape_vaj(axx)
    jxx = re_shape_vaj(jxx)
    vyy = re_shape_vaj(vyy)
    ayy = re_shape_vaj(ayy)
    jyy = re_shape_vaj(jyy)
    
    return vxx,axx,jxx,vyy,ayy,jyy

def mist_2d_vis(waypts,xxs,yys,tts,avj_xxyy,show_wp=True,show_mist_xy=True,show_avj=True):
    """
    This function is to visualize waypts, MiST x,y and MiST avj.

    Parameters
    ----------
    waypts : TYPE
        DESCRIPTION.
    xxs,yys : n * 1 np.array[float]
        X,Y coordinate of generated minimum snap trajectory(MiST).
    tts : n * 1 np.array[float]
        Corresponding timestamp.
    avj_xxyy : tuple
        [vxx,axx,jxx,vyy,ayy,jyy]
    show_wp : bool, optional
        Whether to show waypts on a image. The default is True.
    show_mist_xy : bool, optional
        Whether to show MiST x,y on a image. The default is True.
    show_avj : bool, optional
        Whether to show MiST avj on a image. The default is True.

    Returns
    -------
    None.

    """
    if show_wp == True or show_mist_xy==True or show_avj==True:
        plt.figure()
    else:
        print("Warning: No figure has been plotted!")
        print("To fix: Set the correct show figure parameters")
    
    if show_wp == True:
        plt.plot(waypts[0,:],waypts[1,:],'*r')
        plt.plot(waypts[0,:],waypts[1,:],'b--')
    
    if show_mist_xy==True:
        plt.plot(xxs,yys)
    
    if show_avj == True:
        
        vxx,axx,jxx,vyy,ayy,jyy = avj_xxyy[0],avj_xxyy[1],avj_xxyy[2],avj_xxyy[3],avj_xxyy[4],avj_xxyy[5]
        
        # plot 4 row, 2 column subplots
        row = 4
        col = 2
        fig, ax = plt.subplots(row, col, sharex='col', sharey='row')
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,wspace=None,hspace=1.0)
        tt = tts
        ax[0][0].plot(tt,xxs)
        ax[0][0].set_title('x position')
        ax[0][1].plot(tt,yys)
        ax[0][1].set_title('y position')
        ax[1][0].plot(tt,vxx)
        ax[1][0].set_title('x velocity')
        ax[1][1].plot(tt,vyy,label='y velocity')
        ax[1][1].set_title('y velocity')
        ax[2][0].plot(tt,axx)
        ax[2][0].set_title('x acceleration')
        ax[2][1].plot(tt,ayy)
        ax[2][1].set_title('y acceleration')
        ax[3][0].plot(tt,jxx)
        ax[3][0].set_title('x jerk')
        ax[3][1].plot(tt,jyy,label='y jerk')
        ax[3][1].set_title('y jerk')
        ax[3][0].set_xlabel('time')
        ax[3][1].set_xlabel('time')
        
        for i in range(row):
            for j in range(col):
                ax[i][j].grid(True)        

def main_test():
    ax = [0.0, 1.0, 1.0,4.0, 5.0,8.0]
    ay = [0.0, 2.0, 4.0,8.0, 2.0,3.0]
    
    ax = [0.0, 1.0, 1.5,2.0,4.0, 5.0]
    ay = [0.0, 2.0, 0.5,-1.0,8.0, 2.0]
    
    ax = [0.0, 1.0, 2.0,4.0, 5.0]
    ay = [0.0, 2.0, -1.0,8.0, 2.0]
    
    #ax = [0.0,10.0]
    #ay = [0.0,10]
    waypts_ori = np.array([ax,ay])
    
    T = 12
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    # waypts = np.array([ax])
    waypts = waypts_ori
    ts = arrangeT(waypts,T)
    
    
    n_order = 5
    n_deri = 3
    
    polys_x= minimum_snap_single(waypts[0],ts,n_order,n_deri,v0[0],a0[0],ve[0],ae[0])
    polys_y= minimum_snap_single(waypts[1],ts,n_order,n_deri,v0[1],a0[1],ve[1],ae[1])
    
    interval = 0.01
    
    xxs = np.array([])
    yys = np.array([])
    tts = np.array([])
    
    plt.figure()
    plt.plot(waypts[0,:],waypts[1,:],'*r')
    plt.plot(waypts[0,:],waypts[1,:],'b--')
    for i in range(len(polys_x[0])):
        tt = np.arange(ts[i],ts[i+1],interval)
        xx = polys_vals(polys_x, ts, tt, 0)
        yy = polys_vals(polys_y, ts, tt, 0)
        plt.plot(xx[0],yy[0])
        xxs = np.append(xxs,xx)
        yys = np.append(yys,yy)
        tts = np.append(tts,tt)
    
    tt = tts
    vxx = polys_vals(polys_x,ts,tt,1)
    axx = polys_vals(polys_x,ts,tt,2)
    jxx = polys_vals(polys_x,ts,tt,3)
    vyy = polys_vals(polys_y,ts,tt,1)
    ayy = polys_vals(polys_y,ts,tt,2)
    jyy = polys_vals(polys_y,ts,tt,3)
    

    vxx = re_shape_vaj(vxx)
    axx = re_shape_vaj(axx)
    jxx = re_shape_vaj(jxx)
    vyy = re_shape_vaj(vyy)
    ayy = re_shape_vaj(ayy)
    jyy = re_shape_vaj(jyy)

    """
    # plot 4 row, 2 column subplots
    row = 4
    col = 2
    fig, ax = plt.subplots(row, col, sharex='col', sharey='row')
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9,wspace=None,hspace=1.0)
    
    ax[0][0].plot(tt,xxs)
    ax[0][0].set_title('x position')
    ax[0][1].plot(tt,yys)
    ax[0][1].set_title('y position')
    ax[1][0].plot(tt,vxx)
    ax[1][0].set_title('x velocity')
    ax[1][1].plot(tt,vyy,label='y velocity')
    ax[1][1].set_title('y velocity')
    ax[2][0].plot(tt,axx)
    ax[2][0].set_title('x acceleration')
    ax[2][1].plot(tt,ayy)
    ax[2][1].set_title('y acceleration')
    ax[3][0].plot(tt,jxx)
    ax[3][0].set_title('x jerk')
    ax[3][1].plot(tt,jyy,label='y jerk')
    ax[3][1].set_title('y jerk')
    ax[3][0].set_xlabel('time')
    ax[3][1].set_xlabel('time')
    
    for i in range(row):
        for j in range(col):
            ax[i][j].grid(True)
    """    
def test1():
    ax = [0.0, 1.0, 1.0,4.0, 5.0,8.0]
    ay = [0.0, 2.0, 4.0,8.0, 2.0,3.0]
    waypts_ori = np.array([ax,ay])
    
    T = 5
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    xxs,yys,tts,polys_x,polys_y,ts = mist_2d_gen(waypts_ori,v0,a0,ve,ae,T,None,interval=0.01,n_order=5,n_deri=3)
    vxx,axx,jxx,vyy,ayy,jyy = mist_2d_vaj_gen(xxs,yys,tts,polys_x,polys_y,ts)
        
if __name__ == '__main__':
    ax = [0.0, 1.0, 1.0,4.0, 5.0,8.0]
    ay = [0.0, 2.0, 4.0,8.0, 2.0,3.0]
    waypts_ori = np.array([ax,ay])
    
    T = 5
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    xxs,yys,tts,polys_x,polys_y,ts = mist_2d_gen(waypts_ori,v0,a0,ve,ae,T,None,interval=0.01,n_order=5,n_deri=3)
    avj_xxyy = mist_2d_vaj_gen(xxs,yys,tts,polys_x,polys_y,ts)
    mist_2d_vis(waypts_ori,xxs,yys,tts,avj_xxyy,show_wp=True,show_mist_xy=True,show_avj=True)



    
