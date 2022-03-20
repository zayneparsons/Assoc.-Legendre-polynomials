'''
Created on Mar 9, 2022

@author: zayneparsons
'''
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib import cm
from matplotlib.collections import LineCollection
import math
import sympy as sp
from sympy import latex
from fractions import Fraction as frac 
import numpy as np
def FAC(v):              #Almost positive numpy has a factorial function, but I used this anyways
    H=v
    if v>0:
        for j in range(1,v,1):
            v=v*(H-j)
        return v
    elif v==0:
        return 1
###############################################################################
import tkinter as tk
data= tk.Tk()
data.title('Associated Legendre Polynomials')  #Setup Tkinter window and entry windows
canvas1 = tk.Canvas(data, width = 1000, height = 155)
canvas1.pack()
label1 = tk.Label(data, text='Enter l:')
canvas1.create_window(350, 50, window=label1)
label2 = tk.Label(data, text='Enter m')
canvas1.create_window(350, 75, window=label2)
entry1 = tk.Entry (data) 
canvas1.create_window(500, 50, window=entry1)
entry2 = tk.Entry (data) 
canvas1.create_window(500, 75, window=entry2)
labelleg = tk.Label(data, text='Red indicates when the polynomials are positive (plots are of the absolute value of them); blue indicates when they are negative.')
canvas1.create_window(500, 150, window=labelleg)
Colors=[0]*int((2*np.pi-0)/0.001+1)
X=[None]*int((2*np.pi-0)/0.001+1)            #Create an initial plot of zeros so we can update it after user enters numbers
for i in range(0, int((2*np.pi-0)/0.001)):
    X[i]=0  
Y=[None]*int((2*np.pi-0)/0.001+1)
for i in range(0, int((2*np.pi-0)/0.001)):
    Y[i]=0
exp="$P_{l}^{m}(x)=(-1)^m(1-x^2)^\\frac{m}{2}\\left(\\frac{d}{dx}\\right)^mP_{l}(x)$"       #Adding LaTEX expressions
exp2="$P_{l}^{m}(x)=(-1)^\\frac{3m}{2}\\frac{1}{2^l}\\sum_{k=m}^{l}\\frac{l!(l+m)!}{k!(l+m-k)!(k-m)!(l-k)!}(x-1)^{k-\\frac{m}{2}}(x+1)^{l-k+\\frac{m}{2}}$"
exp3="$P_{l}^{m}(cos(\\theta))=(-1)^m \\frac{1}{2^l} sin^m(\\theta) \\sum_{i=0}^{(k-m)} \\sum_{j=0}^{(l-k)} \\sum_{k=m}^{l} (-1)^i \\frac{l!(l+m)!}{i!j!k!(l+m-k)!(k-m-i)!(l-k-j)!} cos^{(l-m-i-j)}(\\theta)$"
fig = Figure(figsize=(150, 150), dpi=100)           #Adding the initial plot to the window
axis = fig.add_subplot(111, projection='polar')
axis.set_title(exp, fontsize=14)
axis.plot(X, Y, color='black')#
axis.set_thetalim(0,2*np.pi)
canvas2 = FigureCanvasTkAgg(fig, master=data)
canvas2.draw()
toolbar = NavigationToolbar2Tk(canvas2, data)
toolbar.update()
canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.X, expand=1)
fig3 = Figure(figsize=(14, 5), dpi=100)           #Adding the LaTEX expressions
axis3 = fig3.add_subplot(111)
axis3.text(.5, 1, exp2, color="black", fontsize=13, horizontalalignment="center", verticalalignment="top")
axis3.text(.5, .5, exp3, color="black", fontsize=13, horizontalalignment="center", verticalalignment="top")
canvas3 = FigureCanvasTkAgg(fig3, master=data)
axis3.get_xaxis().set_visible(False)
axis3.get_yaxis().set_visible(False)
canvas3.draw()
canvas3.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.X, expand=1)
def AssLeg(v,t, l, m):
        if abs(m)<=l:                                #The Associated Legendre Polynomial. Again There's probably a pre-packaged function but here it is in its raw form
            if m>=0:
                P=0
                for k in range(m, l+1):
                    for i in range(0,k-m+1):
                        for j in range(0,l-k+1):
                            P=P+(-1)**i*frac((FAC(l)*FAC(l+m)*FAC(k-m)*FAC(l-k))/(FAC(i)*FAC(j)*FAC(k)*FAC(l+m-k)*FAC(k-m)*FAC(l-k)*FAC(k-m-i)*FAC(l-k-j)))*v**(l-m-i-j)     
                return (-1)**m*frac(1/2**l)*t**(m)*P
            else:
                M=-m
                P=0
                for k in range(M, l+1):
                    for i in range(0,k-M+1):
                        for j in range(0,l-k+1):
                            P=P+(-1)**i*frac((FAC(l)*FAC(l+M)*FAC(k-M)*FAC(l-k))/(FAC(i)*FAC(j)*FAC(k)*FAC(l+M-k)*FAC(k-M)*FAC(l-k)*FAC(k-M-i)*FAC(l-k-j)))*v**(l-M-i-j)             
                return (-1)**M*frac((1/2**l))*t**(M)*P*FAC(l-M)/FAC(l+M) 
        else:
            return 0
def getdata():   
    global l
    global m
    global U
    global u
    global o  
    global tmptext               #populating values to be used in other functions called when buttons are pressed
    G=entry1.get()
    Gp=entry2.get()
    try:
        l=int(G)
        m=int(Gp)
    except:
        l=None
        m=None
    Endp1=0
    U=float(Endp1)
    Endp2=2*np.pi
    u=float(Endp2)
    step=0.001
    o=float(step)
    ###############################################################################
    ###############################################################################
    x = sp.Symbol('cos(x)')
    y = sp.Symbol('sin(x)')
    if l!=None and m!=None:
        exp=AssLeg(x,y, l, m)
        tmptext = "$"+latex(exp)+"$"  
    ###############################################################################
###############################################################################   
def reg():   
    getdata()
    X=[None]*int((u-U)/o+1)                     #making plots of just the polynomials (or the absolute value of them)
    p=[None]*(int((u-U)/o)+1)
    Colors=[None]*(int((u-U)/o)+1)
    if l!=None and m!=None and abs(m)<=l:
        for i in range(0, int((u-U)/o)+1):
            X[i]=U+i*o
            if abs(m)<=l:
                p[i]=((2*l+1)/(2)*FAC(l-m)/FAC(l+m))**(1/2)*abs(AssLeg(math.cos(X[i]), math.sin(X[i]), l, m))
                Colors[i]=np.sign(AssLeg(math.cos(X[i]), math.sin(X[i]), l, m))  
            else:
                p[i]=None
        segmentsl= []
        for i in range(int((u-U)/o+1)-1):
                segmentsl += [[[X[i],p[i]], [X[i+1],p[i+1]]]]
        segments=np.array(segmentsl)
        if abs(m)<=l:
            axis.cla()       #clearing previous plot
            lc = LineCollection(segments, cmap='coolwarm',array=Colors,linewidth=3)  #creating a set of color coordinated line segments according to the sign of the polynomials
#i.e. if the polynomial is negative, the plot will appear blue. If the absolute value were not plotted, the plot would trace over the same graph twice...
            axis.add_collection(lc)
            O=max(p)                                 #Updating the plot
            axis.set_rlim(0, O)
            axis.set_title(tmptext)
            canvas2.draw()
    else:
        labelack = tk.Label(data, text='Either l or m is not an integer, OR, abs(m)>l')
        canvas1.create_window(750, 50, window=labelack)
        labelack.after(3000, labelack.destroy)       
def mag(): 
    getdata()
    X=[None]*int((u-U)/o+1)                                 #making plots of the magnitude squared
    p=[None]*(int((u-U)/o)+1)
    if l!=None and m!=None and abs(m)<=l:
        for i in range(0, int((u-U)/o)+1):
            X[i]=U+i*o
            if abs(m)<=l:
                p[i]=(((2*l+1)/(2)*FAC(l-m)/FAC(l+m))**(1/2)*AssLeg(math.cos(X[i]), math.sin(X[i]), l, m))**2
            else:
                p[i]=None
        if abs(m)<=l:
            axis.cla() 
            O=max(p)
            axis.plot(X,p, color='red')                                 #Updating the plot
            axis.set_rlim(0, O)
            axis.set_title(tmptext)
            canvas2.draw()
    else:
        labelack = tk.Label(data, text='Either l or m is not an integer, OR, abs(m)>l')
        canvas1.create_window(900, 50, window=labelack)
        labelack.after(3000, labelack.destroy) 
def _quit():
    data.quit()     
    data.destroy()  
button1 = tk.Button(text='Associated Legendre Poly', command=reg)
canvas1.create_window(500, 20, window=button1)
button2=tk.Button(text='Quit/Exit', command=_quit)
canvas1.create_window(500, 125, window=button2)
button3=tk.Button(text='Magnitude squared', command=mag)
canvas1.create_window(500, 101, window=button3)
data.mainloop()
