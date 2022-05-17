import numpy as np
import sys
sys.path.append("..")
from mistgen.mist import mist_generator 

def main_demo_v010():
    # v 0.1.0 test
    ax = [0.0, 1.0, 1.0,4.0, 5.0,8.0]
    ay = [0.0, 2.0, 4.0,8.0, 2.0,3.0]
    # add a classic case in paper
    # "minimum snap trajectory generation and control for quadrotors"
    ax = [0.0, 5.0,5.0,0.0]
    ay = [0.0, 0.0,6.0,6.0]
    
    waypts_ori = np.array([ax,ay])
    
    T = 10
    v0 = np.array([0,0])
    a0 = np.array([0,0])
    ve = np.array([0,0])
    ae = np.array([0,0])
    
    myMistGen = mist_generator()
    xxs,yys,tts = myMistGen.mist_2d_gen(waypts_ori,v0,a0,ve,ae,T)
    vaj_xy = myMistGen.mist_2d_vaj_gen(xxs,yys,tts)
    myMistGen.mist_2d_vis(waypts_ori,xxs,yys,tts,vaj_xy,True,True,True)

if __name__ == '__main__':
    pass
    # main_demo_v010()
