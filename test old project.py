from vpython import *
import statistics as s

dt = 0.001  # timestep
step = 1  # loop counter
maxstep = 1000

#  Define the star, planets and constants
M = 1000  # mass of star (G == 1)
m1 = 1  # mass of planet 1
m2 = 2
initpos1 = vector(0, 1, 0)  # initial position vector of Planet1
initpos2 = vector(0, 2, 0)
Planet1 = sphere(pos=initpos1, radius=0.05 * m1, color=color.blue)
Planet2 = sphere(pos=initpos2, radius=0.05 * m2, color=color.red)
Star = sphere(pos=vector(0, 0, 0), radius=0.1, color=color.yellow)
vel1 = -vector(25, 0, 0)  # initial velocity of planet 1
vel2 = vector(50, 0, 0)
trace1 = curve(color=color.blue)
trace2 = curve(color=color.red)
month = 50
Area = 0
startvec = vector(initpos1)
endvec = vector(0, 0, 0)
dotprod = 0
magnit = 0
numerator = 0
Keplers = list()
averagekep = 0
stdevkep = 0

while step <= maxstep:

    rate(100)  # slow down the animation

    # calculate changes in velocities
    denom1M = mag(Planet1.pos) ** 3
    dv1M = dt * Planet1.pos * M / denom1M
    denom1P = mag(Planet1.pos - Planet2.pos) ** 3
    dv1P = dt * (Planet1.pos - Planet2.pos) * m2 / denom1P

    denom2M = mag(Planet2.pos) ** 3
    dv2M = dt * Planet2.pos * M / denom1M
    denom2P = mag(Planet2.pos - Planet1.pos) ** 3
    dv2P = dt * (Planet2.pos - Planet1.pos) * m1 / denom2P
    # update velocities
    vel1 = vel1 - dv1M - dv1P
    vel2 = vel2 - dv2M - dv2P

    # update positions
    Planet1.pos = Planet1.pos + vel1 * dt
    Planet2.pos = Planet2.pos + vel2 * dt

    step = step + 1
    trace1.append(Planet1.pos)
    trace2.append(Planet2.pos)

    if step % month == 0:
        endvec = vector(Planet1.pos)
        dotprod = dot(startvec, endvec)
        magnit = mag(startvec) * mag(endvec)
        numerator = magnit * acos(dotprod / magnit)
        Area = numerator / 2
        startvec = endvec
        Keplers.append(Area)

averagekep = sum(Keplers)
stdevkep = s.stdev(Keplers)
print("Average area swept per month =", averagekep)
print("Standard deviation in average area =", stdevkep)

print("end of program")
