from tkinter import *
win=Tk()
#function to finnd rk  for gauss jacobi and gasuss sadel
def gjf(eq,x,y,i):
    if i==1:
        res=(eq[3]-(eq[1]*x)-(eq[2]*y))/eq[0]
        return round(res,3)
    elif i==2:
        res=(eq[3]-(eq[0]*x)-(eq[2]*y))/eq[1]
        return round(res,3)
    elif i==3:
        res=(eq[3]-(eq[0]*x)-(eq[1]*y))/eq[2]
        return round(res,3)
#function to do partial piviting
def maxm(x,y,z,i):
    if abs(x[i])>abs(y[i]) and abs(x[i])>abs(z[i]):
        return x
    if abs(y[i])>abs(x[i]) and abs(y[i])>abs(z[i]):
        return y
    if abs(z[i])>abs(x[i]) and abs(z[i])>abs(y[i]):
        return z



#residual function for residual method
def rk(eq,x,y,z):
    re=eq[3]-(eq[0]*x)-(eq[1]*y)-(eq[2]*z)
    return round(re,4)
#largeest residual
def largres(x,y,z):
    if abs(x)>abs(y) and abs(x)>abs(z):
        return 1
    if abs(y)>abs(x) and abs(y)>abs(z):
        return 2
    if abs(z)>abs(x) and abs(z)>abs(y):
        return 3
#increment
def incre(eq,x,i):
    if i==1:
        res=x/eq[0]
    if i==2:
        res=x/eq[1]
    if i==3:
        res=x/eq[3]
    return round(res,4)



    
    
def calc():
    a=[eval(A_5_1.get()),eval(B_6_1.get()),eval(C_7_1.get()),eval(D_8_1.get())]
    b=[eval(A_10_1.get()),eval(B_11_1.get()),eval(C_12_1.get()),eval(D_13_1.get())]
    c=[eval(A_15_1.get()),eval(B_16_1.get()),eval(C_17_1.get()),eval(D_18_1.get())]
    #applying partial piviting   
    e1=maxm(a,b,c,0)
    e2=maxm(a,b,c,1)
    e3=maxm(a,b,c,2)
    x,y,z=0,0,0
    xk,yk,zk=0,0,0
    pp=Label(win,text="The Equations after Partial piviting are:",font=('Arial',12,'bold','underline'))
    pp.grid(row=22,column=0,columnspan=3,sticky="n")
    Eq1=Label(win,text=" "+str(e1[0])+"x+"+str(e1[1])+'y+'+str(e1[2])+"z="+str(e1[3]),font=('Arial',11,'bold'))
    Eq1.grid(row=23,column=0,columnspan=3,sticky='s')

    Eq2=Label(win,text=" "+str(e2[0])+"x+"+str(e2[1])+'y+'+str(e2[2])+"z="+str(e2[3]),font=('Arial',11,'bold'))
    Eq2.grid(row=24,column=0,columnspan=3,sticky='s')

    Eq3=Label(win,text=" "+str(e3[0])+"x+"+str(e3[1])+'y+'+str(e3[2])+"z="+str(e3[3]),font=('Arial',11,'bold'))
    Eq3.grid(row=25,column=0,columnspan=3,sticky='s')
    if(gjc):
        x,y,z=0,0,0
        xk,yk,zk=0,0,0
        for i in range(12):
            xk=gjf(e1,y,z,1)
            yk=gjf(e2,x,z,2)
            zk=gjf(e3,x,y,3)
            lst.append([i,x,y,z,xk,yk,zk])
            x=xk
            y=yk
            z=zk
    if(gse):
        x,y,z=0,0,0
        xk,yk,zk=0,0,0
        for i in range(5):
            xk=gjf(e1,y,z,1)
            yk=gjf(e2,xk,z,2)
            zk=gjf(e3,xk,yk,3)
            lst_s.append([i,x,y,z,xk,yk,zk])
            x=xk
            y=yk
            z=zk

    if(rel):
        x,y,z=0,0,0
        xk,yk,zk=0,0,0
        for i in range(6):
            rx=rk(e1,x,y,z)
            ry=rk(e2,x,y,z)
            rz=rk(e3,x,y,z)
            kl=largres(rx,ry,rz)
            inct=0
            if kl==1:
                inct=incre(e1,rx,kl)
            elif kl==2:
                inct=incre(e2,ry,kl)
            elif kl==3:
                inct=incre(e3,rz,kl)
            lst_r.append([i,x,y,z,rx,ry,rz,inct])
            
            if kl==1:
                x+=inct
            elif kl==2:
                y+=inct
            elif kl==3:
                z+=inct
        
    t = Table(win)
    


class Table:
     
    def __init__(self,root):
        global mat
        
        if(gjc.get()):
            lgjc=Label(root,text="Gauss Jacobi Table",font=('Arial',12,'bold','underline'))
            lgjc.grid(row=mat,column=5,columnspan=7,sticky="n")
            mat=mat+1
            for i in range(mat,12+mat,1):
                for j in range(7):
                    self.e= Entry(root, width=15,font=('Arial',10,'bold'))
                        
                    #self.e.grid_columnconfigure(1, weight=1)
                         
                    self.e.grid(row=i, column=j+5,ipadx=1)
                    self.e.insert(END, lst[i-mat][j])
            mat=mat+12
            
        if(gse.get()):
            lgse=Label(root,text="Gauss Seidal Table",font=('Arial',11,'bold','underline'))
            lgse.grid(row=mat,column=5,columnspan=7,sticky="n")
            mat=mat+1
            for i in range(mat,5+mat,1):
                for j in range(7):
                     
                    self.e = Entry(root, width=15,font=('Arial',10,'bold'))
                    
                    #self.e.grid_columnconfigure(1, weight=1)
                     
                    self.e.grid(row=i, column=j+5,ipadx=1)
                    self.e.insert(END, lst_s[i-mat][j])
            mat=mat+5
        if(rel.get()):
            lrel=Label(root,text="Residual Table",font=('Arial',12,'bold','underline'))
            lrel.grid(row=mat,column=5,columnspan=7,sticky="n")
            mat=mat+1
            for i in range(mat,6+mat,1):
                for j in range(8):
                     
                    self.e = Entry(root, width=15,font=('Arial',10,'bold'))
                    
                    #self.e.grid_columnconfigure(1, weight=1)
                     
                    self.e.grid(row=i, column=j+5,ipadx=1)
                    self.e.insert(END, lst_r[i-mat][j])
            mat=mat+6
            
 
# take the data
lst=[['k','x','y','z',"x(k+1)",'y(k+1)','z(k+1)']]
lst_s=[['k','x','y','z',"x(k+1)",'y(k+1)','z(k+1)']]
lst_r=[['k','x','y','z',"Rx",'Ry','Rz','increment']]

  
# find total number of rows and
# columns in list



head=Label(win,text='EQUATION CALCULATOR',font=('Arial',15,'bold','underline'))
head.grid(row=0,columnspan=13,sticky="n",ipadx=50)

k=0

gjc=IntVar()
gj=Checkbutton(win,text="Gauss Jacobi",variable=gjc,onvalue=1,offvalue=0)
gj.grid(row=1,column=0,sticky='w')

gse=IntVar()
gs=Checkbutton(win,text="Gauss Seidal",variable=gse,onvalue=1,offvalue=0)
gs.grid(row=2,column=0,sticky='w')

rel=IntVar()
rx=Checkbutton(win,text="Residual",variable=rel,onvalue=1,offvalue=0)
rx.grid(row=3,column=0,sticky='w')
j=4

for i in range(1,4,1):
    exec(f"e_{i}_{j}_{k}=Label(win,text='Enter coeficients of equation {i} in the form Ax+By+Cz=D')")
    exec(f"e_{i}_{j}_{k}.grid(row={j},columnspan=2,sticky='w',ipadx=10)")
    j=j+1;
    exec(f"e_{j}_{k}=Label(win,text='\tA{i}=')")
    exec(f"e_{j}_{k}.grid(row={j},column={k},sticky='e')")
    k+=1
    exec(f"A_{j}_{k}=Entry(win,border=5)")
    exec(f"A_{j}_{k}.grid(row={j},ipadx=1,columnspan=1,column={k},sticky='w')")
    k=0
    j+=1
    exec(f"e_{j}_{k}=Label(win,text='\tB{i}=')")
    exec(f"e_{j}_{k}.grid(row={j},column={k},sticky='e')")
    k+=1
    exec(f"B_{j}_{k}=Entry(win,border=5)")
    exec(f"B_{j}_{k}.grid(row={j},columnspan=1,column={k},sticky='w')")
    k=0
    j+=1
    exec(f"e_{j}_{k}=Label(win,text='\tC{i}=')")
    exec(f"e_{j}_{k}.grid(row={j},column={k},sticky='e')")
    k+=1
    exec(f"C_{j}_{k}=Entry(win,border=5)")
    exec(f"C_{j}_{k}.grid(row={j},columnspan=1,column={k},sticky='w')")
    k=0
    j+=1
    exec(f"e_{j}_{k}=Label(win,text='\tD{i}=')")
    exec(f"e_{j}_{k}.grid(row={j},column={k},sticky='e')")
    k+=1
    exec(f"D_{j}_{k}=Entry(win,border=5)")
    exec(f"D_{j}_{k}.grid(row={j},columnspan=1,column={k},sticky='w')")
    k=0
    j+=1


mat=2
        
en=Button(win,text='Enter',command=lambda:calc())
en.grid(row=j,columnspan=3,ipadx=5,sticky='ns')
j+=1



        

    
        



