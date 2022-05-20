# -*- coding: utf-8 -*-
"""
@author: zhaoliang (zhz03@g.ucla.edu)
"""

import numpy as np
import math
from utils.T_functions import arrangeT
from utils.T_functions import init_T
from utils.cvxopt_qp import quadprog
from utils.other_utils import computeQ
from utils.other_utils import calc_tvec
from utils.other_utils import re_shape_vaj
from utils.other_utils import poly_val
from scipy.linalg import block_diag
import matplotlib.pyplot as plt

class mist_generator():
    def __init__(self,ts=np.array([]),interval=0.01,n_order=5,n_deri=3):

        
        self.interval = interval if interval == None else interval
        self.n_order = n_order if n_order == None else n_order
        self.n_deri = n_deri if n_deri == None else n_deri
        self.polys_x = None
        self.polys_y = None
        self.ts = np.array([]) if ts.size == 0 else ts
        
    def minimum_snap_single(self,waypts,ts,v0,a0,ve,ae):
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
        n_coef = self.n_order+1
        
        # n_deri = 3
        # compute Q
        Q_all = np.empty(0)
        
        for i in range(n_poly):
            # print(i)
            q_new = computeQ(self.n_order,self.n_deri,ts[i],ts[i+1])
            # print(len(q_new))
            Q_all = block_diag(Q_all,q_new)
        Q_all = np.delete(Q_all,0,axis=0)
    
        if len(Q_all) != len(Q_all[0]):
            raise Exception('Q_all is not a square matrix')
            
        b_all = np.zeros((len(Q_all),1))
    
        Aeq = np.zeros((4*n_poly+2,n_coef*n_poly))
        beq = np.zeros((4*n_poly+2,1))
    
        # start/terminal pva constraints  (6 equations)
        A1 = calc_tvec(ts[0],self.n_order,0)
        A2 = calc_tvec(ts[0],self.n_order,1)
        A3 = calc_tvec(ts[0],self.n_order,2)
        Aeq[0,0:n_coef] = A1
        Aeq[1,0:n_coef] = A2
        Aeq[2,0:n_coef] = A3
        
        Ae1 = calc_tvec(ts[-1],self.n_order,0)
        Ae2 = calc_tvec(ts[-1],self.n_order,1)
        Ae3 = calc_tvec(ts[-1],self.n_order,2)
        
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
            Aeq[neq-1,col1-1:col2] = calc_tvec(ts[i+1],self.n_order,0)
            beq[neq-1] = waypts[i+1]
    
        # continuous constraints  ((n_poly-1)*3 equations)
        
        for i in range(n_poly-1):
            tvec_p = calc_tvec(ts[i+1],self.n_order,0)
            tvec_v = calc_tvec(ts[i+1],self.n_order,1)
            tvec_a = calc_tvec(ts[i+1],self.n_order,2)
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
    
    def polys_vals(self,polys,ts,tt,r):
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
    
    def calc_yaw(self,vy,vx):
        """
        Calculate yaw radius and yaw angle

        Parameters
        ----------
        vy,vx : n * 1 np.array[float]
            y velocity and x velocity.

        Returns
        -------
        yaw_rad,yaw_deg : n * 1 np.array[float]
            yaw radius and yaw angle

        """
        lxy = len(vx)
        yaw_rad = np.zeros((lxy,))
        yaw_deg = np.zeros((lxy,))
        for i in range(lxy):
            yaw_i = math.atan2(vy[i],vx[i])
            yaw_rad[i] = yaw_i
            yaw_deg[i] = np.rad2deg(yaw_i)
        # correct the first one
        yaw_deg[0] = yaw_deg[1]
        yaw_rad[0] = yaw_rad[1]
            
        return yaw_rad,yaw_deg
    
    def minimum_snap_2d(self,waypts,ts,v0,a0,ve,ae):
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
        polys_x = self.minimum_snap_single(waypts[0],ts,v0[0],a0[0],ve[0],ae[0])
        polys_y = self.minimum_snap_single(waypts[1],ts,v0[1],a0[1],ve[1],ae[1])
        return polys_x,polys_y
    
    def mist_2d_gen(self,waypts,v0,a0,ve,ae,T=None,ave_speed=None):
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
        """
        if waypts == None:
            print('Warning: waypoints are not defined!')
            print('Example default waypoints will be assigned.')
        """
        if waypts.size != 0:
            waypts = waypts
        else:
            waypts = np.array([[0.0, 1.0, 1.0,4.0, 5.0],[0.0, 2.0, 4.0,8.0, 2.0]])
        v0 = v0 if v0.size != 0 else np.array([0,0])
        a0 = a0 if a0.size != 0 else np.array([0,0])
        ve = ve if ve.size != 0 else np.array([0,0])
        ae = ae if ae.size != 0 else np.array([0,0])
        
        if T == None and ave_speed != None:
            T = init_T(waypts,ave_speed)
        elif T == None and ave_speed == None:
            print('Warning: Both Total time {T} and ave_speed are not defined!')
            print('Example default T will be assigned.')
            T = 10
        else:
            T = T            
        
        if self.ts.size == 0:
            self.ts = arrangeT(waypts,T)
        
        self.polys_x,self.polys_y = self.minimum_snap_2d(waypts,self.ts,v0,a0,ve,ae)
        
        xxs = np.array([])
        yys = np.array([])
        tts = np.array([])
        
        for i in range(len(self.polys_x[0])):
            tt = np.arange(self.ts[i],self.ts[i+1],self.interval)
            xx = self.polys_vals(self.polys_x, self.ts, tt, 0)
            yy = self.polys_vals(self.polys_y, self.ts, tt, 0)
            # plt.plot(xx[0],yy[0])
            xxs = np.append(xxs,xx)
            yys = np.append(yys,yy)
            tts = np.append(tts,tt)
            
        return xxs,yys,tts
    
    def mist_2d_vaj_gen(self,xxs,yys,tts):
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
        vxx = self.polys_vals(self.polys_x,self.ts,tt,1)
        axx = self.polys_vals(self.polys_x,self.ts,tt,2)
        jxx = self.polys_vals(self.polys_x,self.ts,tt,3)
        vyy = self.polys_vals(self.polys_y,self.ts,tt,1)
        ayy = self.polys_vals(self.polys_y,self.ts,tt,2)
        jyy = self.polys_vals(self.polys_y,self.ts,tt,3)

        vxx = re_shape_vaj(vxx)
        axx = re_shape_vaj(axx)
        jxx = re_shape_vaj(jxx)
        vyy = re_shape_vaj(vyy)
        ayy = re_shape_vaj(ayy)
        jyy = re_shape_vaj(jyy)
        
        return vxx,axx,jxx,vyy,ayy,jyy
    
    def mist_2d_vis(self,waypts,xxs,yys,tts,vaj_xxyy,show_wp=True,show_mist_xy=True,show_avj=True,same_plot=False):
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
        vaj_xxyy : tuple
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
            if same_plot == False:
                plt.figure()
        else:
            print("Warning: No figure has been plotted!")
            print("To fix: Set the correct show figure parameters")
        
        if show_wp == True:
            plt.plot(waypts[0,:],waypts[1,:],'*r',label='waypoints')
            plt.plot(waypts[0,:],waypts[1,:],'b--')
            plt.legend()
        
        if show_mist_xy==True:
            plt.plot(xxs,yys,label='MiST')
            plt.legend()
        
        if show_avj == True:
            
            vxx,axx,jxx,vyy,ayy,jyy = vaj_xxyy[0],vaj_xxyy[1],vaj_xxyy[2],vaj_xxyy[3],vaj_xxyy[4],vaj_xxyy[5]
            
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
            ax[3][0].set_xlabel('time(s)')
            ax[3][1].set_xlabel('time(s)')
            
            for i in range(row):
                for j in range(col):
                    ax[i][j].grid(True)    
    def mist_2d_gen_all(self,waypts_ori,v0,a0,ve,ae,T,show_wp=True,show_mist_xy=True,show_avj=True,same_plot=False):
        """
        This function is to put MiST x,y,v,a,j output and visualization all together.
        
        Parameters
        ----------
        waypts_ori : TYPE
            DESCRIPTION.
        v0 : np.array[float] 
            Initial set velocity.
        a0 : np.array[float] 
            Initial set acceleration.
        ve : np.array[float] 
            End set velocity.
        ae : np.array[float] 
            End set acceleration.
        T : float
            Total set execution time.
        show_wp : TYPE, optional
            DESCRIPTION. The default is True.
        show_mist_xy : TYPE, optional
            DESCRIPTION. The default is True.
        show_avj : TYPE, optional
            DESCRIPTION. The default is True.
        same_plot : TYPE, optional
            DESCRIPTION. The default is False.
        
        Returns
        -------
        xxs : TYPE
            DESCRIPTION.
        yys : TYPE
            DESCRIPTION.
        tts : TYPE
            DESCRIPTION.
        vaj_xy : TYPE
            DESCRIPTION.
        self.ts : np.array[float]
            Time segment.
        
        """
        xxs,yys,tts= self.mist_2d_gen(waypts_ori,v0,a0,ve,ae,T)
        vaj_xy = self.mist_2d_vaj_gen(xxs,yys,tts)
        self.mist_2d_vis(waypts_ori,xxs,yys,tts,vaj_xy,show_wp,show_mist_xy,show_avj,same_plot)
        
        return xxs,yys,tts,vaj_xy,self.ts   
    
if __name__ == '__main__':  
    pass
