---
jupytext:
  formats: notebooks//ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
```

# Computational Mechanics Project #01 - Heat Transfer in Forensic Science

We can use our current skillset for a macabre application. We can predict the time of death based upon the current temperature and change in temperature of a corpse. 

Forensic scientists use Newton's law of cooling to determine the time elapsed since the loss of life, 

$\frac{dT}{dt} = -K(T-T_a)$,

where $T$ is the current temperature, $T_a$ is the ambient temperature, $t$ is the elapsed time in hours, and $K$ is an empirical constant. 

Suppose the temperature of the corpse is 85$^o$F at 11:00 am. Then, 2 hours later the temperature is 74$^{o}$F. 

Assume ambient temperature is a constant 65$^{o}$F.

1. Use Python to calculate $K$ using a finite difference approximation, $\frac{dT}{dt} \approx \frac{T(t+\Delta t)-T(t)}{\Delta t}$.

```{code-cell} ipython3
dt = 2 #hours
T0 = 85 #deg F
T2 = 74 #deg F
T_a = 65 #deg F

dTdt = (T2-T0)/dt

Kback = dTdt/-(T2-T_a)
Kforw = dTdt/-(T0-T_a)

print('backward diff deriv: K = {:1.2f}' .format(Kback))
print('forward diff deriv: K = {:1.2f}' .format(Kforw))
```

2. Change your work from problem 1 to create a function that accepts the temperature at two times, ambient temperature, and the time elapsed to return $K$.

```{code-cell} ipython3
def Kvalue(T0, T2, Ta, dt):
    '''Takes initial temperature, change in time, temperature at time and 
    ambient temperature to return a K backwards and forwards value
    
    returns array of first backwards then forwards K constant'''
    dTdt = (T2-T0)/dt
    Kback = dTdt/-(T2-Ta)
    Kforw = dTdt/-(T0-Ta)
    
    return Kback, Kforw
```

```{code-cell} ipython3
dt = 2 
T0 = 85 
T2 = 74 
Ta = 65 


Kvalue(T0, T2, Ta, dt)
```

3. A first-order thermal system has the following analytical solution, 

    $T(t) =T_a+(T(0)-T_a)e^{-Kt}$

    where $T(0)$ is the temperature of the corpse at t=0 hours i.e. at the time of discovery and $T_a$ is a constant ambient temperature. 

    a. Show that an Euler integration converges to the analytical solution as the time step is decreased. Use the constant $K$ derived above and the initial temperature, T(0) = 85$^o$F. 

    b. What is the final temperature as t$\rightarrow\infty$?
    
    c. At what time was the corpse 98.6$^{o}$F? i.e. what was the time of death?

```{code-cell} ipython3
t = np.linspace(0, 2, 500)
dt = t[1]- t[0]
T_amb = 65
K = Kback

T_eul_f = np.zeros(len(t))
T_eul_b = np.zeros(len(t))

T_eul_f[0] = 85
T_eul_b[0] = 85

for i in range(1, len(t)):
    T_eul_f[i] = T_eul_f[i-1] - K*(T_eul_f[i-1] - T_amb)*dt
    T_eul_b[i] = T_eul_b[i-1] + K*(T_eul_b[i-1] - T_amb)*dt

T_eul_b[T_eul_b > 98.6] = 98.6
    
plt.plot(-t, T_eul_b,  label='Before Discovery')
plt.plot(t, T_eul_f, label='After Discovery')
plt.ylabel('Body Temperature ($^\circ F$)')
plt.xlabel('Time (hours)')
plt.legend(loc='best')
plt.grid;
```

Part b

As t $\rightarrow\infty$, $e^{\infty}\rightarrow 0$ . Due to this the temperature of the corpse reaches the ambient pressure which is 65 degrees.

```{code-cell} ipython3
time = np.where(T_eul_b==98.6)
tt=time[0][0]

time_0 = 11
time_0_minutes = 11*60
time_dead = t[tt]*60

time = (time_0_minutes-time_dead)/60
time

timehr = int(time)
timesec = time-round(time,0)
timemin = round(timesec*60,2)

print('Time of Death is',timehr,':',timemin,'AM')
```
