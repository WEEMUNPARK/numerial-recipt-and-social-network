#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
import matplotlib.pyplot as plt

def polint(xa, ya, n, x, y, dy):
    # initialize c and d to be equal to ya

    c = ya.copy()
    d = ya.copy()
    # then find an index ns which is closest to the value x in xa()
    ns = 0
    dif = math.fabs(x - xa[0])
    for i in range(1, n):
        dift = math.fabs(x - xa[i])
        if (dift < dif):
            ns = i
            dif = dift
    y[0] = ya[ns]
    ns -= 1

    # do double loop over column m and row i on the triangular Tableau.
    # we don't try to catch the dividing by 0 error, let Python itself do it
    for m in range(1, n):
        for i in range(0, n - m):
            ho = xa[i] - x
            hp = xa[i + m] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            den = w / den
            d[i] = hp * den
            c[i] = ho * den
        # end for i
        # tricking coding here in C
        if (2 * (ns + 1) < (n - m)):
            dy[0] = c[ns + 1]
        else:
            dy[0] = d[ns]
            ns -= 1
        # end if else
        y[0] += dy[0]
    # end for m
    for e in y:
        return e
    


# end def polint(...)

# a test run
xa = [-1, 1, 2, 4]
ya = [1.25, 2, 3, 0]
y = [0]
dy = [0]
n = 4

def y_range(x,y,xa,ya,n,dy):
    a=polint(x,y,xa,ya,n,dy)
    return a
x1=np.linspace(-1,4,100)
y1=[]
for i in x1:
    y1.append(y_range(xa,ya,n,i,y,dy))
plt.plot(x1,y1)
plt.show()


# In[2]:


Problem 2


# In[3]:


import numpy

def romberg( f, a, b, n ):

    s=(n+1,n+1)
    r = numpy.zeros(s,float)
    h = b - a
    r[0][0] = 0.5 * h * ( f( a ) + f( b ) )

    powerOf2 = 1
    for i in range( 1, n + 1 ):

        # Compute the halved stepsize and use this to sum the function at
        # all the new points (in between the points already computed)

        h = 0.5 * h

        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in range( 1, powerOf2, 2 ):
            sum = sum + f( a + k * h )

        # Compute the composite trapezoid rule for the next level of
        # subdivision.  Use Richardson extrapolation to refine these values
        # into a more accurate form.

        r[i][0] = 0.5 * r[i-1][0] + sum * h

        powerOf4 = 1
        for j in range( 1, i + 1 ):
            powerOf4 = 4 * powerOf4
            r[i][j] = r[i][j-1] + ( r[i][j-1] - r[i-1][j-1] ) / ( powerOf4 - 1 )
            
    return r[-1][-1]


# In[4]:


f=lambda x: (x**4)*math.log(x+(x**2+1)**0.5)   #here print out the result list which was the result from n=1 to 29
a=0
b=2
n=29
result_list=[]
for n in range(0,30):
    x=romberg(f,a,b,n)
    result_list.append(x)
result_list


# In[ ]:


result_list.pop(0) #since n=0 is not needed we delete it then by the way
result_list


# In[ ]:



import math
diff=[]
for e in result_list:
    diff.append(e-(8 -40*math.sqrt(5) + 480*math.asinh(2))/75)
#this step was to turn the error between the accurate values and the experimental values into absolute ones for better estimating

for i in range(0,len(diff)):
    if diff[i]<0:
        diff[i]=diff[i]*-1

diff


# In[ ]:



#this was to show the minimum error of the errors
def pickerjo(diff):
    if len(diff)>1:
        if diff[1]>=diff[0]:
            diff.pop(1)
        elif diff[1]<diff[0]:
            diff.pop(0)
        return pickerjo(diff)
    if len(diff)==1:
        return diff


# In[ ]:


pickerjo(diff)
result_a=diff[0]


# In[ ]:


print(result_a+(8 -40*math.sqrt(5) + 480*math.asinh(2))/75) #print out the most accurate one we got in


# In[ ]:





# In[ ]:


Problem 3


# In[7]:


import math
import random
import numpy as np
import matplotlib.pyplot as plt
x=np.zeros(6)
def func(r,x):    #set the wave function's *opp wave function
    abs_r1 = math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])+math.sqrt((x[3]-r)*(x[3]-r)+x[4]*x[4]+x[5]*x[5])
    abs_r2 = math.sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+math.sqrt((x[0]-r)*(x[0]-r)+x[1]*x[1]+x[2]*x[2])
    func = math.pow(math.e, -(1/a0)*abs_r1)+math.pow(math.e, -(1/a0)*abs_r2)
    
    return(func)


# In[8]:


import random

def set_new_x(x):
    x1=[0,0,0,0,0,0]
    for i in range(len(x1)):
        x1[i]= x[i]+(2*random.random()-1)*a0  #set the range from (-a0,a0)
    return x1

 # since the Xi was six elements, and we need N times six elements, to find the most stable ones
def sample_made(r):
    x = np.zeros(6)
    N = 15000
    samples = []            
    while len(samples) < N: 
        
        x_new= set_new_x(x)
        
        P_oriirr =(func(r,x)**2)
        P_monver=(func(r,x_new)**2)
        ratio=min(P_oriirr,P_monver/P_oriirr)    #since it is a probability distribution with wavefunction**2, so since Markov chain P(n)=Pn**n, we can see that 
        u=random.random()
        
        if u<=ratio:               
            
            x = x_new
            samples.append(x)
    
    return samples


# In[ ]:


def get_E_Average(r):

    Kenergy= 0                     
    Potential= 0

    for i in range(0,len(sample)):   #set the matrix for list N-dimensions list filled with x0,...x5
        x = np.zeros(6)
        for j in range(len(x)):
            x[j] = sample[i][j]                        
                                               

        
        #the Laplacian part of Knetic energy were
        r_A1 = math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])   
        r_B1 = math.sqrt((x[0]-r)**2+x[1]**2+x[2]**2) 
        r_A2 = math.sqrt(x[3]**2+x[4]**2+x[5]**2)  
        r_B2 = math.sqrt((x[3]-r)**2+x[4]**2+x[5]**2)
        
        K_particle1_1=(np.exp((-1/a0)*(r_A1+r_B2))/a0)*((x[0]**2+x[1]**2+x[2]**2)/a0-3*r_A1+(x[0]**2+x[1]**2+x[2]**2)/r_A1)/r_A1**2
        K_particle1_2=(np.exp((-1/a0)*(r_A2+r_B1))/a0)*(((x[0]-r)**2+x[1]**2+x[2]**2)/a0-3*r_B1+((x[0]-r)**2+x[1]**2+x[2]**2)/r_B1)/r_B1**2
        K_1=K_particle1_1 + K_particle1_2
        
        
        #also, particle2 is similar
        K_particle2_1=(np.exp((-1/a0)*(r_A1+r_B2))/a0)*(((x[3]-r)**2+x[4]**2+x[5]**2)/a0-3*r_B2+((x[3]-r)**2+x[4]**2+x[5]**2)/r_B2)/r_B2**2
        K_particle2_2=(np.exp((-1/a0)*(r_A2+r_B1))/a0)*((x[3]**2+x[4]**2+x[5]**2)/a0-3*r_A2+(x[3]**2+x[4]**2+x[5]**2)/r_A2)/r_A2**2
        K_2=K_particle2_1 + K_particle2_2
       
       
        
        
        r_12 = math.sqrt((x[0]-x[3])*(x[0]-x[3])+(x[1]-x[4])*(x[1]-x[4])+(x[2]-x[5])*(x[2]-x[5]))
        
        Ken = (-1/2)*(K_1+K_2)        
        Poten = 1/r_12 + 1/r - 1/r_A1 - 1/r_A2 -1/r_B1 - 1/r_B2

        Kenergy += Ken/func(r,sample[i])   #sum P(x)E(x)=total E   
        Potential += Poten                 #similar to upon one      
    
    E_Average = Kenergy/(len(sample)) + Potential/(len(sample))  # here is the average energy
    
    return(E_Average)



r_range = np.linspace(0.5,8,1000)    


res = []                   
a0 = 0.8   #change the a0 to adjust the graph
for r in r_range :
    sample = []
    sample = sample_made(r)
    res.append(get_E_Average(r))

plt.plot(r_range,res)          
plt.xlabel('r')
plt.ylabel('Energy')


# In[ ]:


def frequencyfind(x):
    a=[]
    for key in x:
        if key not in a:
            a.append(key)
    return a 
def ff(a,x):
    count=0
    for i in range(len(x)):
        if x[i]==a:
            count+=1
    return count
def fin(x):
    a=frequencyfind(x)
    b=[]
    for e in a:
        b.append(ff(e,x))
    return b
def findi(x):
    a=frequencyfind(x)
    b=fin(x)
    c=dict(zip(a,b))
    return c
        
        


# In[17]:


x=['a','b','b','b','c',7,6,4,8,6,7,5,5,2,3,2,1,5,8,9,0,'o',90909090,10000000000000]
del x[0:3]
print(x)


# In[11]:


_CHO_ = 'ㄱㄲㄴㄷㄸㄹㅁㅂㅃㅅㅆㅇㅈㅉㅊㅋㅌㅍㅎ'
_JUNG_ = 'ㅏㅐㅑㅒㅓㅔㅕㅖㅗㅘㅙㅚㅛㅜㅝㅞㅟㅠㅡㅢㅣ'
_JONG_ = 'ㄱㄲㄳㄴㄵㄶㄷㄹㄺㄻㄼㄽㄾㄿㅀㅁㅂㅄㅅㅆㅇㅈㅊㅋㅌㅍㅎ'
_cho_=list(_CHO_)
_JUNG_=list(_JUNG_)
_JONG_l=list(_JONG_)
a=_cho_[0]+_JUNG_[0]
print(_JONG_l)
#class find():
    #def model(self,x):
        #for i in range(len(x)):
            #if x[0]


# In[37]:


x=['ㄱ','ㅏ','ㅇ','ㅏ']
fi(x)


# In[66]:


a='123'
a=list(a)
a

