

if __name__ == '__main__':  
    ax = [0.0, 1.0,1.0,0.0]
    ay = [0.0, 0.0,2.0,2.0]
    
    #ax = [0.0, -4.0]
    #ay = [0.0, -4.0]    
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

    vxx,axx,jxx,vyy,ayy,jyy = vaj_xxyy[0],vaj_xxyy[1],vaj_xxyy[2],vaj_xxyy[3],vaj_xxyy[4],vaj_xxyy[5]
    
    yaw_rad,yaw_deg = myMistGen.calc_yaw(vyy,vxx)
