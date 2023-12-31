#!/usr/bin/env python
# coding: utf-8

# In[1]:


#problem 1


# In[10]:


import math
def ludcmp(a, n, indx, d): #define the ludcmp func
    d[0] = 1.0
    vv=indx.copy()
    for i in range(0,n):
        big = 0.0
        for j in range(0,n):
            temp = math.fabs(a[i][j])
            if (temp > big):
                big = temp
        vv[i]=1.0/big
    for j in range(0,n): #the whole progress was considered by the crout algorythm
        big = 0.0
        for i in range(0,n):           
            if(i<j):
                l = i
            else:
                l = j
            sum = a[i][j]    
            for k in range(0,l):
                sum -= a[i][k]*a[k][j] #sum of alpha ik*betha kj in crout algorythm
            
                a[i][j] = sum
           
            if(i>=j):
                dum = vv[i]*math.fabs(sum)
                if (dum >= big):
                    big = dum
                    imax = i
                    
        # pivoting part, swap row j with row imax, a[j] is a whole row
        if (j != imax):
            dum = a[imax]
            a[imax] = a[j]
            a[j] = dum
            d[0] = - d[0]
            vv[imax] = vv[j]
            
        # divide by the beta diagonal value
        indx[j] = imax
        dum = 1.0/a[j][j]
        for i in range(j+1,n):
            a[i][j] *= dum

# define lubksb() take the row swapped LU decomposed

def lubksb(a, n, indx, b):
    ii = -1
    # forward
    for i in range(0,n):
        ip = indx[i]
        sum = b[ip]
        b[ip] = b[i]
        if(ii != -1):
            for j in range(ii,i):
                sum -= a[i][j]*b[j]
        elif(sum != 0):
            ii = i
        b[i] = sum
        
    #  backward
    for i in range(n-1,-1,-1):
        sum = b[i]
        for j in range(i+1,n):
            sum -= a[i][j]*b[j]
        b[i] = sum/a[i][i]



def solution(a,n,b):#define the function to solve the equation
    indx=[]
    for e in range(n):
        indx.append(e) #fill the indx list
    d =[]
    d.append(1)#let d=[1], which is a list with element 1
    ludcmp(a,n,indx,d)  #use the ludcmp func as the problem mentioned
    x = b.copy()
    lubksb(a,n,indx,x) #use the ludcmp func as the problem mentioned
    return x

unknown_matrix=[[1.083565459610028+34,-1],[2,3]]
n=2
value=[3.3481894150417832,5]
solution(unknown_matrix,n,value)


# In[13]:


#problem 2


# In[70]:


import numpy #set a list for len(b)=k and for e in b: len(e)=k (here k=L+1 since L=2, it is a 3*3 mattei)
def bize(n):
    s=((n**2-2),(n**2-2))
    b=numpy.zeros(s,int)
    return b  


# In[71]:


#b=bize(k) #set a list with len=n**2-2, and every element is a list with n**2-2's 0
#b
#the thought for solution is firstly set a set of unknown data for (k**2-2) which delete the V11 and Vkk cause V11=0
#Vkk=1, and then set (k**2-2) equations to solve them
#after receive every node Voltage, minus the connected voltages with V11 or Vkk and take their sum, since r=1, then can get their current.
#finally the resistance, as well
#make sure this approach must guarantee that the equality for current to V11 and current from Vkk
#V21-V11+V12-V11=Vkk-Vkk-1+Vkk-Vk-1k

k=int(input())+1
print('the enter number turns to matrix'+ str(k))


# In[72]:


b=bize(k)
if k==2:
    b=[[2,0],[0,2]]
if k==3:
    b[0][0]=3
    b[0][1]=-1
    b[0][k]=-1
    b[-1][-1]=3
    b[-1][-2]=-1
    b[-1][-1-k]=-1
    b[k-2][k-1]=3
    b[k-2][k]=-1
    b[k-2][k-1+k]=-1
    b[-k+1][-k-1]=-1
    b[-k+1][-k]=3
    b[-k+1][-2*k]=-1
    b[2*k-4][k-2]=2
    b[2*k-4][k-3]=-1
    b[2*k-4][k-2+k]=-1
    b[-2*k+3][-k+1]=2
    b[-2*k+3][-k+1-k]=-1
    b[-2*k+3][-k+2]=-1
    b[1+k+k-4][1+k-1]=4
    b[1+k+k-4][1-1]=-1
    b[1+k+k-4][1+k-2]=-1
    b[1+k+k-4][1+k]=-1
    b[1+k+k-4][1+k+k-1]=-1
if k>3:
    b=bize(k)
    b[0][0]=3
    b[0][1]=-1
    b[0][k]=-1 #3V21-V11-V31-V22=0 here V11=0
    def ni(b,i):#set the first set of function 3Vi1-V(i-1)1-V(i+1)1-Vi2=0 and the first solution x was [V21,31,41,51,61,71,81,91]
        b[i][i-1]=-1
        b[i][i]=3  
        b[i][i+1]=-1
        b[i][i+k]=-1
        return b
    def chanumber(b,i):
    
        b=ni(b,i)
        if i<k-3:
            return chanumber(b,i+1)
        elif i>=k-3:
            return b  #so the formula will cost 7equations

    b[-1][-1]=3
    b[-1][-2]=-1
    b[-1][-1-k]=-1#3Vk-1k-Vk-2k-Vk-1k-1=1
    def nk(b,i):  #here is the last column,  3Vik-V(i-1)k-V(i+1)k-Vi(k-1)=0 the last element for solution x was[........V1k,2k,3k,4k,5k,6k,7k,8k....k-1k]
        b[-i-1][-i-2]=-1  
        b[-i-1][-i-1]=3   #so the formula will cost from b[-7to-1]
        b[-i-1][-i]=-1
        b[-i-1][-i-1-k]=-1
        return b
    def nechanumber(b,i):
        b=nk(b,i)
        if i<k-3:
            return nechanumber(b,i+1)
        elif i>=k-3:
            return b #so the formula will also cost 7equations
    
#repeat the previous process by lines
#since the solutions are in the order[V21,31,41,51,61,71,81.....k1,12,22,32,42............13,23,33,43.............k-1k]
#the first column 3V1i-V1i+1-V1i-1-V2i=0(k>i>1)
#since 3V12-V11-V13-V22=0 V11=0
    b[k-2][k-1]=3
    b[k-2][k]=-1
    b[k-2][k-1+k]=-1
    def revbi(b,i):
        b[i+k-2][k+(i-1)*k-1]=-1
        b[i+k-2][k+(i+0)*k]=-1
        b[i+k-2][k+(i+0)*k-1]=3
        b[i+k-2][k+(i+1)*k-1]=-1
        return b
    def rechanumber(b,i):
        b=revbi(b,i)
        if i<k-3:
            return rechanumber(b,i+1)
        if i>=k-3:
            return b
#next step return the similar process for k column around all lines
#3Vkj-Vkj-1-Vkj+1-Vk-1j=0
#think of 3Vkk-1-Vkk-2-Vkk-Vk-1k-1=0 since Vkk=1, 3Vkk-1-Vkk-2-Vk-1k-1=1 canresult
    b[-k+1][-k-1]=-1
    b[-k+1][-k]=3
    b[-k+1][-2*k]=-1
    def revnk(b,i):
        b[-i-k+1][-(i+1)*k-1]=-1
        b[-i-k+1][-(i+1)*k]=3
        b[-i-k+1][-(i+2)*k]=-1
        b[-i-k+1][-i*k]=-1
        return b
    def revnkchanumber(b,i):
        b=revnk(b,i)
        if i<k-3:
            return revnkchanumber(b,i+1)
        elif i>=k-3:
            return b
#for the next step find the top 2 points except the V11(V11=0),and V99(V99=1)
#2Vk1-Vk-11-Vk2=0 also including. 2V1k-V1k-1-V2k=0
    b[2*k-4][k-2]=2
    b[2*k-4][k-3]=-1
    b[2*k-4][k-2+k]=-1
    b[-2*k+3][-k+1]=2
    b[-2*k+3][-k+1-k]=-1
    b[-2*k+3][-k+2]=-1

#final equations for nodes with 4 connections which was 4Vij-Vij-1-Vi-1j-Vij+1-Vi+1j=0

#j with i in range[2,k-1] wile like the previous functions i j start with 2

    def fochkijf(b,i,j):
        b[i+k+k-4][j+k-1]=4
        b[i+k+k-4][j-1]=-1
        b[i+k+k-4][j+k-2]=-1
        b[i+k+k-4][j+k]=-1
        b[i+k+k-4][j+k+k-1]=-1
        return b
    def fchanumber(b,i,j):
        b=fochkijf(b,i,j)
        if i<(k-2)**2:
            if i%(k-2)!=0:
                return fchanumber(b,i+1,j+1)
            elif i%(k-2)==0:
                return fchanumber(b,i+1,j+3)
        else:
            if i>=(k-2)**2:
                return b
#print(b)



# In[73]:


if k==2:
    final_unknown=b
    final_matrix=numpy.array(final_unknown)
    x_value=[1,1]
if k==3:
    final_matrix=numpy.array(b)
    final_unknown=b.tolist()
    x_value=[0]*(3**2-2) #set a b value matr while the -1 and -k+1 is 1, others are 0 Ax=b
    x_value[-k+1]=1
    x_value[-1]=1
if k>3:
    set_1=chanumber(b,1)
    set_2=nechanumber(set_1,1)
    set_3=rechanumber(set_2,1)
    set_4=revnkchanumber(set_3,1)

    final_matrix=fchanumber(set_4,1,1)# set an array for Vi+1i to Vk-1k total (k**2-2)

    final_unknown=final_matrix.tolist()#transform the array to list

#print(final_unknown)
    x_value=[0]*(k**2-2) #set a b value matr while the -1 and -k+1 is 1, others are 0 Ax=b
    x_value[-k+1]=1
    x_value[-1]=1


# In[74]:


import math
import time
def ludcmp(a, n, indx, d):
    d[0] = 1.0
    # looking for the largest a in each row and store it in vv as inverse
    # We need a new list same size as indx, for this we use .copy() 
    vv=indx.copy()
    for i in range(0,n):
        big = 0.0
        for j in range(0,n):
            temp = math.fabs(a[i][j])
            if (temp > big):
                big = temp
        vv[i]=1.0/big
    #
    # run Crout's algorithm
    for j in range(0,n):
        # top half & bottom part are combined
        # but the upper limit l for k sum is different
        big = 0.0
        for i in range(0,n):           
            if(i<j):
                l = i
            else:
                l = j
            sum = a[i][j]    
            for k in range(0,l):
                sum -= a[i][k]*a[k][j]
            # end for k
            a[i][j] = sum
            # for bottom half, we keep track which row is larger
            if(i>=j):
                dum = vv[i]*math.fabs(sum)
                if (dum >= big):
                    big = dum
                    imax = i
            # end if (i>= ...)
        # end for i
        # pivoting part, swap row j with row imax, a[j] is a whole row
        if (j != imax):
            dum = a[imax]
            a[imax] = a[j]
            a[j] = dum
            d[0] = - d[0]
            vv[imax] = vv[j]
        # end if (j != ...)
        # divide by the beta diagonal value
        indx[j] = imax
        dum = 1.0/a[j][j]
        for i in range(j+1,n):
            a[i][j] *= dum
        # end for i 
    # end for j 
# end of def ludcmp

# We do backward substitution in lubksb() take the row swapped LU decomposed
# a, size n, and swapping indx, and b vector as input.  The output is
# in b after calling.
def lubksb(a, n, indx, b):
    ii = -1
    # forward
    for i in range(0,n):
        ip = indx[i]
        sum = b[ip]
        b[ip] = b[i]
        if(ii != -1):
            for j in range(ii,i):
                sum -= a[i][j]*b[j]
        elif(sum != 0):
            ii = i
        b[i] = sum
    # bote alpha_{ii} is 1 above
    #  backward
    for i in range(n-1,-1,-1):
        sum = b[i]
        for j in range(i+1,n):
            sum -= a[i][j]*b[j]
        b[i] = sum/a[i][i]
# end lubksb()

            
# unfortunately a is destroyed (become swapped LU)
def linearsolver(a,n,b):
    indx = list(range(n))
    d =[1]
    ludcmp(a,n,indx,d)
    x = b.copy()
    lubksb(a,n,indx,x)
    print("x=",x)
    return x
# end linearsolver

n=k**2-2
unknownV=linearsolver(final_unknown,n,x_value) #set the values of voltages as a list for unknownV

#since the resistor network is k*k square(L+1=k, so we just need to focus on the V11 and Vkk, and the current is also
#current=V21-V11+V12-V11=Vkk-Vkk-1+Vk-Vk-1k) which means we need to make sure the equation is verified, andthe value was current
current_1=unknownV[0]+unknownV[k-1]
current_2=2-unknownV[-1]-unknownV[-k]

current_difference=current_1-current_2
#we set a range that if the difference is far less than 0.0001, we can say that they are equal because the result has error while calcualting
if current_1-current_2<0.0001: #justify the equation
    print('equal')

time.process_time()#note the time


# In[75]:


import scipy,time
import numpy as np


# In[76]:


from scipy import linalg
if k==2:
    h=[1,1]
if k>2:
    h=[0]*(k**2-2)#set a b value matr while the -1 and -k+1 is 1, others are 0 Ax[-k+1]=1  
    h[1-k]=1
    h[-1]=1

final_maarrey=np.array(final_matrix)
scipy.linalg.solve(final_maarrey,h)
print(scipy.linalg.solve(final_matrix,h)) #here we can not use the list final_unknown, because of the setting for scipy


#since the resistor network is k*k square(L+1=k, so we just need to focus on the V11 and Vkk, and the current is also
#current=V21-V11+V12-V11=Vkk-Vkk-1+Vk-Vk-1k) which means we need to make sure the equation is verified, andthe value was current
current_1=unknownV[0]+unknownV[k-1]
current_2=2-unknownV[-1]-unknownV[-k]

current_difference=current_1-current_2
#we set a range that if the difference is far less than 0.0001, we can say that they are equal because the result has error while calcualting
if current_1-current_2<0.0001: #justify the equation
    print('equal')
    

time.process_time()#we must use the matrix set which is final_matrix


# In[77]:


current_total=2-unknownV[-1]-unknownV[-k]
print('current is'+' '+ str(current_total))
resistance_total=1/current_total
print('resistance is'+' '+ str(resistance_total))
print('end of lab')
#V11=[0]
#totalV=V11+unknownV#here since we just solve the V from(Vi+1i to Vk-1k), 
#totalV.append(1)
#totalV #final total voltage with V11=0 in the first element and Vkk=1 in the last


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




