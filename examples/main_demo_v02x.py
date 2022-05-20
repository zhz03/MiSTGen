import numpy as np
import sys
sys.path.append("..")
from mistgen.mist import mist_generator 
from mistgen.utils.T_functions import arrangeT

def main_demo_v020():
    ax = [0.0, 1.0,1.0,0.0]
    ay = [0.0, 0.0,2.0,2.0]
    
    waypts_ori = np.array([ax,ay])
    
    T = 10
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    ts = arrangeT(waypts_ori,T)
    ts[1] = ts[1]
    ts[2] = ts[2]    

    myMistGen = mist_generator(interval=0.01,n_order=5,n_deri=3)
    xxs,yys,tts = myMistGen.mist_2d_gen(waypts_ori,v0,a0,ve,ae,T)
    vaj_xy = myMistGen.mist_2d_vaj_gen(xxs,yys,tts)
    myMistGen.mist_2d_vis(waypts_ori,xxs,yys,tts,vaj_xy,True,True,False,False)
    
    ts[1] = ts[1] + 0.5
    ts[2] = ts[2] - 0.5
    myMistGen = mist_generator(ts,interval=0.01,n_order=5,n_deri=3)
    xxs,yys,tts = myMistGen.mist_2d_gen(waypts_ori,v0,a0,ve,ae,T)
    vaj_xy = myMistGen.mist_2d_vaj_gen(xxs,yys,tts)
    myMistGen.mist_2d_vis(waypts_ori,xxs,yys,tts,vaj_xy,True,True,False,True) 

def main_demo_v021():
    ax = [0.0, 1.0,1.0,0.0]
    ay = [0.0, 0.0,2.0,2.0]
    
    waypts_ori = np.array([ax,ay])
    
    T = 10
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    ts = arrangeT(waypts_ori,T)
    ts[1] = ts[1]
    ts[2] = ts[2] 
    
    myMistGen = mist_generator()
    mistTraj1 = myMistGen.mist_2d_gen_all(waypts_ori,v0,a0,ve,ae,T,show_wp=True,show_mist_xy=True,show_avj=False,same_plot=False)
    
    ts = mistTraj1[-1]
    ts[1] = ts[1] + 0.5
    ts[2] = ts[2] - 0.5
    
    myMistGen = mist_generator(ts)
    mistTraj2 = myMistGen.mist_2d_gen_all(waypts_ori,v0,a0,ve,ae,T,show_wp=True,show_mist_xy=True,show_avj=False,same_plot=True)
    
    ts = mistTraj2[-1]
    ts[1] = ts[1] + 0.3
    ts[2] = ts[2] - 0.3
    
    myMistGen = mist_generator(ts)
    mistTraj3 = myMistGen.mist_2d_gen_all(waypts_ori,v0,a0,ve,ae,T,show_wp=True,show_mist_xy=True,show_avj=False,same_plot=True)
    
    ts = mistTraj3[-1]
    ts[1] = ts[1] + 0.2
    ts[2] = ts[2] - 0.2
    
    myMistGen = mist_generator(ts)
    myMistGen.mist_2d_gen_all(waypts_ori,v0,a0,ve,ae,T,show_wp=True,show_mist_xy=True,show_avj=False,same_plot=True)

def main_demo_v024():
    ax = [0.0, 1.0,1.0,0.0]
    ay = [0.0, 0.0,2.0,2.0]
     
    waypts_ori = np.array([ax,ay])
    T = 5
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    ts = arrangeT(waypts_ori,T)
   
    myMistGen = mist_generator(ts,interval=0.01,n_order=5,n_deri=3)
    xxs,yys,tts = myMistGen.mist_2d_gen(waypts_ori,v0,a0,ve,ae,T)
    vaj_xxyy = myMistGen.mist_2d_vaj_gen(xxs,yys,tts)
    myMistGen.mist_2d_vis(waypts_ori,xxs,yys,tts,vaj_xxyy,True,True,True,False) 

    vxx,_,_,vyy,_,_ = vaj_xxyy[0],vaj_xxyy[1],vaj_xxyy[2],vaj_xxyy[3],vaj_xxyy[4],vaj_xxyy[5]
    
    yaw_rad,yaw_deg = myMistGen.calc_yaw(vyy,vxx)  
    
if __name__ == '__main__':
    main_demo_v020()