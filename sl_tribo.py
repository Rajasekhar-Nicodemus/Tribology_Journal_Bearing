import streamlit as st
import matplotlib.pyplot as plt
import numpy as np



app_mode = st.sidebar.selectbox('Select Page',['Journal Bearing']) 
st.sidebar.markdown("![Alt Text](https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Hydrodynamic_lubrication_attitude_angle.svg/330px-Hydrodynamic_lubrication_attitude_angle.svg.png)")
st.sidebar.markdown("Image sorce: Wikipedia")

def pressurefem(ex,ey,rad,wid,cl,visc,g_a,g_c,U):
    
    h=np.zeros([g_c.shape[0],g_c.shape[1]])
    dh=np.zeros([g_c.shape[0],g_c.shape[1]])
    for i in range(0,g_c.shape[0]):
        for j in range(0,g_c.shape[1]):
            h[i,j]=hfilm(cl,ex,ey,g_c[i,j])
            dh[i,j]=dhfilm(cl,ex,ey,g_c[i,j])/rad
  

    p=np.zeros([g_c.shape[0],g_c.shape[1]])
    
    delx=(g_c[1,0]-g_c[0,0])*rad
    delz=g_a[0,1]-g_a[0,0]
    
    delx2=1/(delx*delx)
    delz2=1/(delz*delz)
    
    sum2=delz2+delx2
    
    for z in range(0,400):
        for i in range(0,g_c.shape[0]):
          for j in range(0,g_c.shape[1]):
              
              if j==0 or j==g_c.shape[1]-1:
                  
                  p[i,j]=0
              else:     
                 
                  A=(delx2-1.5*dh[i,j]/(delx*h[i,j]))/(2*sum2)
                  B=(delz2)/(2*sum2)
                  C=(delx2+1.5*dh[i,j]/(delx*h[i,j]))/(2*sum2)
                  D=(delz2)/(2*sum2)
                  E=-(3*visc*U*dh[i,j])/(sum2*h[i,j]**3)
                  
                  if i==0:
                      k=g_c.shape[0]-1                      
                      p[i,j]=A*p[int(k),j]+B*p[i,j-1]+C*p[i+1,j]+D*p[i,j+1]+E
                  elif i==g_c.shape[0]-1:  
                      p[i,j]=A*p[i-1,j]+B*p[i,j-1]+C*p[0,j]+D*p[i,j+1]+E
                  else:
                      p[i,j]=A*p[i-1,j]+B*p[i,j-1]+C*p[i+1,j]+D*p[i,j+1]+E
                      
              if p[i,j]<0.0:
                   p[i,j]=0.0
              

    hmin=np.min(np.min(h))
    pmax=np.max(np.max(p))    
    return(hmin,pmax,h,p)
   
def hfilm(c,ex,ey,th):
    return(c-ex*np.cos(th)-ey*np.sin(th))

def dhfilm(c,ex,ey,th):
    return(ex*np.sin(th)-ey*np.cos(th))

def loadintegral(p,delx,delz,g_c):
     px=np.zeros([g_c.shape[0],g_c.shape[1]])
     py=np.zeros([g_c.shape[0],g_c.shape[1]])
     
     sum_px=0
     sum_py=0
     
     for i in range(0,g_c.shape[0]):
        for j in range(0,g_c.shape[1]):
            px[i,j]=p[i,j]*np.cos(g_c[i,j])
            py[i,j]=p[i,j]*np.sin(g_c[i,j])
            
            
     for i in range(0,g_c.shape[0]):
        for j in range(0,g_c.shape[1]):
               if j==0 or j==g_c.shape[1]-1 or i==0 or  i==g_c.shape[0]-1: 
                 sum_px=sum_px+2*px[i,j]   
                 sum_py=sum_py+2*py[i,j]
               else:  
                 sum_px=sum_px+4*px[i,j]   
                 sum_py=sum_py+4*py[i,j]
                 
     sum_px=sum_px-px[0,0]-px[g_c.shape[0]-1,g_c.shape[1]-1]-px[0,g_c.shape[1]-1]-px[g_c.shape[0]-1,0]          
     sum_py=sum_py-py[0,0]-py[g_c.shape[0]-1,g_c.shape[1]-1]-py[0,g_c.shape[1]-1]-py[g_c.shape[0]-1,0]
     
     sum_px=sum_px*delx*delz*0.25
     sum_py=sum_py*delx*delz*0.25
     return(sum_px,sum_py)


if app_mode=='Journal Bearing':
    st.header("Simple  Journal Bearing APP: By E.Rajasekhar Nicodemus (rajasekhar.nicodemus@gmail.com)")
    st.header("Bearing Geomteric parameters")
    rad=st.number_input("Bearing Radius (mm)",0.01,100.0)
    rad=rad/(1000)
    wid=st.number_input("Bearing Width (mm)",0.01,100.0)
    wid=wid/1000
    cl=st.number_input("Bearing Clearance (microns)",0.1,100.0)
    cl=cl/(1000)
    cl=cl/(1000)
    
    st.header("Oil properties")
    visc=st.number_input("Dynamic viscoity (cp)",1.00,500.0)
    visc=visc*0.001
    
    #st.write(rad,wid,cl,visc)
    
    st.header("Bearing Operating conditions")
    speed=st.number_input("Journal Speed (rpm)",100.0,10000.0)
    speed=speed*0.10472
    wx=st.number_input("Bearing Load in X-Horizontal(N)",0.0,10000.0)
    wy=st.number_input("Bearing Load in Y-Vertical (N)",0.0,10000.0)
    #st.write(wx,wy)
    
    st.header("Mesh details")
    na=st.number_input("Nodes in axial direction",5,30)
    nc=st.number_input("Nodes in circumferential direction",60,250)
    lim=st.number_input("Number of load iterations",40,500)
    
    U=speed*rad
    #st.write(speed,rad,U)
    ga=np.linspace(0,1,int(na))
    gc=np.linspace(0,1,int(nc))
    gc=gc[0:len(gc)-1]
    ga=ga*wid
    gc=2*np.pi*gc
    g_a, g_c = np.meshgrid(ga, gc)
     
    
    g_cp=g_c*rad
    delx=(g_c[1,0]-g_c[0,0])*rad
    delz=g_a[0,1]-g_a[0,0]
    
    #ex=0.8025*cl
    #ey=0.4152*cl
    
    
    if st.button("Compute Bearing Performance"):
      cou=0
      ch=0
      ex=0
      ey=0
      status= st.empty()
      
      while ch==0:

        cou=cou+1
        
        (hmin,pmax,h,p)=pressurefem(ex,ey,rad,wid,cl,visc,g_a,g_c,U)
        (fx,fy)=loadintegral(p,delx,delz,g_c)
        
        (hmin,pmax,h,p)=pressurefem(ex+1e-11,ey,rad,wid,cl,visc,g_a,g_c,U)
        (fxx,fyx)=loadintegral(p,delx,delz,g_c)
        
        (hmin,pmax,h,p)=pressurefem(ex,ey+1e-11,rad,wid,cl,visc,g_a,g_c,U)
        (fxy,fyy)=loadintegral(p,delx,delz,g_c)
        
        Kxx=(fxx-fx)/(1e-11)
        Kyy=(fyy-fy)/(1e-11)
        Kxy=(fyx-fy)/(1e-11)
        Kyx=(fxy-fx)/(1e-11)
        
        #st.write(Kxx,Kyy,Kxy,Kyx,fx,fy)
        #st.write(cou,fx,fy)
        c1=wx-fx
        c2=wy-fy
        
        den=Kxx*Kyy-Kxy*Kyx
        
        delex=(c1*Kyy-c2*Kyx)/den
        deley=(c2*Kxx-c1*Kxy)/den
        
   
        #st.write(ex,ey)
        if(cou<=20):
           ex=ex+0.01*delex
           ey=ey+0.01*deley
        else:   
           ex=ex+0.2*delex
           ey=ey+0.2*deley
           
        #st.write(ex/cl,ey/cl)
        status.write("Load Iteration :"+str(cou)+" Error in load-X (N) : "+str(c1)+" Error in Load-Y(N) : "+str(c2))
        
        exc=ex/cl
        eyc=ey/cl
        eff_ecc=np.sqrt(exc*exc+eyc*eyc) 
        
        if(cou==lim):
            ch=1
            st.warning("Max iteraions exhausted. Check the bearing paramters or Increase max load iterations")
            
        if(eff_ecc>=0.99):
            ch=1
            st.error("Bearing will go into metal to metal contact.Please check bearing parametrs")
    
        metric_x=1
        metric_y=1
        
        if wx==0 and np.abs(c1)<=10:
            metric_x=0
            
        if wy==0 and np.abs(c2)<=10:
            metric_y=0    
            
        if wx>0 and np.abs(c1/wx)<=0.01:
            metric_x=0
            
        if wy>0 and np.abs(c2/wy)<=0.01:
            metric_y=0  
            
        if metric_x==0 and metric_y==0:
            ch=1
            st.success("Load itreations converged")    
   
           
            
                          
      (hmin,pmax,h,p)=pressurefem(ex,ey,rad,wid,cl,visc,g_a,g_c,U)
      (fx,fy)=loadintegral(p,delx,delz,g_c)
      #st.write(ex,ey)
      #st.write(fx,fy)
      
      st.header("Bearing Results") 
      
      st.write("Bearing eccentricty ratio in x direction: "+ str(ex/cl))
      st.write("Bearing eccentricty ratio in y direction: "+ str(ey/cl))
      
      st.write("Load carried in x direction: "+ str(fx))
      st.write("Load carried in y direction: "+ str(fy))
     
      
      fig_1=plt.figure(figsize=(10, 4))
      ax=fig_1.add_subplot(111)
      cp=ax.contourf(g_cp, g_a, h,cmap='PuBuGn')
      ax.set_title('Gap Thickness Countour Plot')
      ax.set_ylabel('Axial (mm)')
      ax.set_xlabel('Circumferential (mm)')
      #plt.scatter(g_cp,g_a,marker='x')
      fig_1.colorbar(cp)
      st.pyplot(fig_1)
      st.write("Minimum film thickness:"+str(hmin))
      
      
      fig_2=plt.figure(figsize=(10, 4))
      ax2=fig_2.add_subplot(111)
      cp2=ax2.contourf(g_cp, g_a, p,cmap='gist_rainbow')
      ax2.set_title('Pressure Countour Plot')
      ax2.set_ylabel('Axial (mm)')
      ax2.set_xlabel('Circumferential (mm)')
      #plt.scatter(g_cp,g_a,marker='x')
      fig_2.colorbar(cp2)
      st.pyplot(fig_2)
      st.write("Maximum film pressure:"+str(pmax))

  



    

    
    

