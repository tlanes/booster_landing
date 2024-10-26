import numpy as np
import matplotlib.pyplot as plt

# environment const
g = 9.81  # m/s^2

# time
tf = 1000  # s
dt = 0.01  # s
t_ = np.arange(0, tf, dt)
len_t_ = len(t_)

# preallocate
tuza = np.zeros(len_t_)
tuzv = np.zeros(len_t_)
tunf = np.zeros(len_t_)
vv = np.zeros(len_t_)  # vertical velocity
av = np.zeros(len_t_)  # vertical accel
mf = np.zeros(len_t_)  # fuel mass
y = np.zeros(len_t_)   # height off ground
throttle = np.zeros(len_t_)  # throttle [0, 1]
mdot = np.zeros(len_t_)      # mass flow rate
tw_rat = np.zeros(len_t_)

# mass and fuel
mw0 = 75000  # kg
mf0 = 45500  # kg
md = mw0 - mf0
mdot0 = 137  # kg/s
# mdot0 = 237  # kg/s
tunf0 = mf0/mdot0  # s
thrust0 = 490e4  # N
throttle0 = 0  # N/s

# ic
y0 = 50e3  # m
vv0 = 0  # m/s
av0 = -g

vv[0] = vv0 + av0 * dt
av[0] = av0
tunf[0] = tunf0
mf[0] = mf0
y[0] = y0
throttle[0] = throttle0
thrust = np.ones(len_t_) * thrust0
thrust_act = np.zeros(len_t_)
thrust_act[0] = 0
mdot[0] = 0
tw_rat[0] = thrust0 / (mw0 * g)
for idx, t_i in enumerate(t_):

    mw = mf[idx - 1] + md
    if idx == 0:
        continue
    else:
        # f = m*a -> a = f/m
        # a = -g + (thrust*throttle)/(mw + mdot*dt)
        # av[idx] = -g + (thrust[idx - 1]*throttle[idx - 1])/(mw + mdot[idx - 1]*dt)
        av[idx] = -g + (thrust[idx - 1] * throttle[idx - 1]) / mw
        # v = v0 + a*t
        vv[idx] = vv[idx - 1] + av[idx - 1]*dt
        # x = x0 + v0*t + 0.5*a*t^2
        y[idx] = y[idx - 1] + vv[idx-1]*dt + (1/2)*av[idx-1]*dt**2

        # calculate Time Until variables
        # Time until zero altitude : 0 = y0 + v*t + 1/2*a*t^2 ==> solve for t
        tuza[idx] = (np.sqrt(vv[idx]**2 - 2*y[idx]*av[idx]) - vv[idx])/av[idx]
        if tuza[idx] < 0:
            tuza[idx] = -(np.sqrt(vv[idx] ** 2 - 2 * y[idx] * av[idx]) + vv[idx]) / av[idx]

        # time until zero velocity : 0 = v0 + a*t ==> solve for t
        tuzv[idx] = np.abs(vv[idx]/av[idx])

        # time until zero fuel
        if mdot[idx - 1] > 0.0:  # if you're firing
            mf[idx] = mf[idx - 1] - mdot[idx - 1]*dt
            tunf[idx] = mf[idx]/mdot[idx - 1]
        else:  # if you're not firing
            mf[idx] = mf[idx - 1]
            tunf[idx] = tunf[idx - 1]

    # logic for decision-making on thrust
    if tunf[idx] > tuza[idx] and tunf[idx] > tuzv[idx]:  # if TUNF is greater than both TUZA and TUZV
        # these are unused
        dtuzv = tuzv[idx] - tuzv[idx - 1]
        dtuza = tuza[idx] - tuza[idx - 1]
        sdtuzv = np.sign(dtuzv)
        sdtuza = np.sign(dtuza)
    else:
        throttle[idx] *= 0.98
        mdot[idx] *= 0.98

    if tuzv[idx] < tuza[idx]:  # you're gonna hover
        throttle[idx] = throttle[idx - 1] - 0.03 * np.abs(tuzv[idx] - tuza[idx])
        throttle[idx] = np.max([throttle[idx], 0.0])
        mdot[idx] = mdot0 * throttle[idx]
    elif tuzv[idx] > tuza[idx]:  # you're gonna splat
        throttle[idx] = throttle[idx - 1] + 0.3 * np.abs(tuzv[idx] - tuza[idx])
        throttle[idx] = np.min([throttle[idx], 1.0])
        mdot[idx] = mdot0 * throttle[idx]

    # else:
    #     throttle[idx] *= 0
    #     mdot[idx] *= 0

    thrust_act[idx] = thrust[idx] * throttle[idx]

    if y[idx] <= 0.0:
        y[idx] = 0.0
        print(f'hit the ground at t = {t_i}:\nav_f = {av[idx]}m/s^2, vv_f = {vv[idx]}m/s')
        break

pad_xlim = 10
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
ax3.plot(t_, av, marker='.', linestyle='none')
ax3.set_xlabel('time, s')
ax3.set_ylabel('accel, m/s^2')
ax3.grid()
ax3.set_xlim([0, t_i+pad_xlim])

ax2.plot(t_, vv, marker='.', linestyle='none')
ax2.set_xlabel('time, s')
ax2.set_ylabel('velocity, m/s')
ax2.grid()
ax2.set_xlim([0, t_i+pad_xlim])

ax1.plot(t_, y, marker='.', linestyle='none')
ax1.set_xlabel('time, s')
ax1.set_ylabel('y, m')
ax1.grid()
ax1.set_xlim([0, t_i+pad_xlim])

fig, ax = plt.subplots()
ax.plot(t_, tunf, marker='.', linestyle='none', label='tunf')
ax.plot(t_, tuza, marker='.', linestyle='none', label='tuza')
ax.plot(t_, tuzv, marker='.', linestyle='none', label='tuzv')
ax.grid()
ax.set_xlabel('time, s')
ax.set_ylabel('time until [], s')
ax.legend()
ax.set_xlim([0, t_i+pad_xlim])

fig, ax = plt.subplots()
ax.plot(t_, thrust_act, marker='.', linestyle='none')
ax.grid()
ax.set_xlabel('time, s')
ax.set_ylabel('thrust_act, N')
ax.set_xlim([0, t_i+pad_xlim])

fig, ax = plt.subplots()
ax.plot(t_, throttle, marker='.', linestyle='none')
ax.grid()
ax.set_xlabel('time, s')
ax.set_ylabel('throttle, []')
ax.set_xlim([0, t_i+pad_xlim])

fig, ax = plt.subplots()
ax.plot(t_, mf, marker='.', linestyle='none')
ax.grid()
ax.set_xlabel('time, s')
ax.set_ylabel('mf, kg')
ax.set_xlim([0, t_i+pad_xlim])

plt.show()


