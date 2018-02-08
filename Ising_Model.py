# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 17:49:48 2017

@author: Colm
"""

#ISING MODEL

import numpy as np
import matplotlib.pylab as plt
import math

#setting up a matrix of -1s and 1s
def matgen(a,b): #inputs should be: a*b = size of system
    global s
    s =[[np.random.choice(range(-1,2,2)) for x in range(a+1)]for y in range(b+1)] 




def mainloop(a,b,T,it):
    #initial parameters
    i=np.arange(0,a)
    j=np.arange(0,b)
    Kb=1.38e-23
    t=T/1.38e-23
    B=1./(Kb*t)
    E2=0.
    S2=0.
    avE=0.
    avS=0.
    c=0
    h=0
    global Slist1,avElist
    #setting up empty matrices for computed variables to be iterably stored in
    Slist1=[0 for x in range(it)]
    avElist1=[0 for x in range(it)]
    z=np.arange(0,it)
    for element in z:
        E2=0
        S2=0
        for g in i:    #iterating for each element in a row 
            for element in j: #iterating for each element in a column
                
                #setting periodic boundary conditions
                if g>=a-1:
                    s[g+1][element]=s[0][element]         
                if g<=0:
                    s[g-1][element]=s[a][element]
                if element>=b-1:
                    s[g][element+1]=s[g][0]                
                if element<=0:
                    s[g][element-1]=s[g][b]
                
                #Calculating the energy at each lattice point
                E=np.float(2*s[g][element]*(s[g-1][element]+s[g+1][element]+s[g][element-1]+s[g][element+1])) -h*s[g][element]
                E2=(E2+E) 
                global E2 
                
                #metropolis algorithm
                if E<0.:
                    s[g][element]=s[g][element] * -1
                if E>0.:
                    if np.exp(-B*(E))>np.random.random():
                        s[g][element]=s[g][element] * -1
                global S2        
                S2=S2 + np.float(s[g][element])
        
        #calculating variables and assigning them to a matrix position.
        global avE
        avE=E2/(a*b)
        global avS
        avS=S2/(a*b)
        Slist1[c]=avS
        if Slist1[c]<0:
           matgen(a,b)
        avElist1[c]=avE
        c=c+1                    
        



#printing the matrix s.
def matprint():
    plt.matshow(s)
    plt.title('System size:[a*b]')
    plt.xlabel('a')
    plt.ylabel('b')





#defining a function that iterates the mainloop over a range of temperatures and saves the variables.
def ising(a,b,temp1,temp2,its,it): #inputs should be:a*b=system size, temp1->temp2=temperature range,
                             # its=np. of steps in temp range, it=number of iterations that the system will run through
    global y
    y=np.arange(temp1,temp2,(temp2-temp1)/its)
    c=0
    for i,element in enumerate(y): # enumerating the array y so that each element is referable.
        if i==0:
                #initialising matrices for storing variables.
                global avSlist,avElist,chilist,Cvlist
                avSlist=[0. for x in range(its)]
                Elist=[0. for x in range(its)]
                avElist=[0. for x in range(its)]
                Slist=[0. for x in range(its)]
                chilist=[0. for x in range(its)]
                Cvlist=[0. for x in range(its)]
        
        #defining the temperature before the first figure for the temperature so that 
        # Cv & chi can be calculated for the first element in the temperature array. 
        if i<=0:
            y[i-1]=y[i]-(temp2-temp1)/its
        #generating the matrix and running the main steps.
        matgen(a,b)
        mainloop(a,b,element,it)
        #saving the final average energy,magnetisation etc when the system has 
        #reached equalibrium.
        avSlist[c]=avS
        avElist[c]=avE
        Elist[c]=E2
        Slist[c]=S2
        global Cv
        Cv=(Elist[c]-Elist[c-1])/(element-y[i-1])
        #defining Cv & chi for temperatures previous to the first element in temp
        if i<=0:
            Cv=0
        global chi
        chi=(Slist[c]-Slist[c-1])/(element-y[i-1])
        if i<=0:
            chi=0
        chilist[c]=np.abs(chi)
        Cvlist[c]=np.abs(Cv)
        c=c+1
    #plotting each variable Vs temp
    plots(temp1,temp2,its)




def plots(temp1,temp2,its):
    temp=np.arange(temp1,temp2,(temp2-temp1)/its)
    #absolute value of the average magnetisation(avS) and average energy(avE)
    w=np.abs(avSlist)
    w1=np.abs(avElist)
    plt.scatter(temp,w)
    plt.plot(temp,w)
    plt.title('Temp Vs Magnetisation')
    plt.xlabel('Temp [J/Kb]')
    plt.ylabel('<|Magnetic spin|> [1/n]')
    plt.show()
    plt.scatter(temp,w1)
    plt.plot(temp,w1)
    plt.title('Temp Vs Energy')
    plt.xlabel('Temp [J/Kb]')
    plt.ylabel('<|Energy|> [1/n]')
    plt.figure()
    plt.scatter(temp,Cvlist)
    plt.plot(temp,Cvlist)
    plt.title('Temp Vs Cv(Specific Heat Capacity)')
    plt.xlabel('Temp [J/Kb]')
    plt.ylabel('<Specific Heat Capacity> ')
    plt.figure()
    plt.scatter(temp,chilist)
    plt.plot(temp,chilist)
    plt.title('Temp Vs Magnetic Suceptability')
    plt.xlabel('Temp [J/Kb]')
    plt.ylabel('<|magnetic spin|> (chi)')
    plt.figure()
    plt.show()




#plotting how the amount of iterations to reach equalibrium depends on system size.
def syssize(xy1,xy2,step): #inputs should be:(start,stop,step) for the range of system sizes.
    a=np.arange(xy1,xy2,step)
    global avSlist
    avSlist=[0. for x in range((xy1-xy2)/step)]
    c=0
    for h in a:
        matgen(h,h)
        # The temp is set to 1 and thus should reach equalibrium at approx. avS=1.0
        mainloop(h,h,1,1000)
        avSlist[c]=avS
        c=c+1
    #Plots the system size vs the average Magnetisation per site
    syssizeplot(xy1,xy2,step)




def syssizeplot(xy1,xy2,step):
    xy=np.arange(xy1,xy2,step)
    plt.plot(xy,avSlist)
    plt.title('Effect of System Size on the Results')   
    plt.xlabel('System Size[n*n]')
    plt.ylabel('Average Magenetisation per site[1/n]' )
    plt.show()


