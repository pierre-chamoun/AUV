# %%
# libraries
import numpy as np
import matplotlib.pyplot as plt

# %%
# constant parameters of the thruster positions
a = 15
b = 7.5
c = 15
d = 15
e = 15

num_thrusters = 7
ndof = 6 

# %%
# computation of thruster forces while varying alpha
steps = 10 
delta_alpha = 90 / (steps - 1)
alpha = np.arange(0, 90, delta_alpha)
alpha = np.append(alpha, alpha[(steps - 1) - 1] + delta_alpha)
alpha = np.deg2rad(alpha)

unit_thruster_forces = np.zeros((steps ,num_thrusters, ndof))

for i in range(steps):
    TCM = np.array([[np.cos(alpha[i]), np.cos(alpha[i]), -np.cos(alpha[i]), -np.cos(alpha[i]), 0, 0, 0],
                    [np.sin(alpha[i]), -np.sin(alpha[i]), np.sin(alpha[i]), -np.sin(alpha[i]), 0, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1],
                    [0, 0, 0, 0, -c, c, 0],
                    [0, 0, 0, 0, -b, -b, a],
                    [d*np.cos(alpha[i]) + e*np.sin(alpha[i]), -d*np.cos(alpha[i]) - e*np.sin(alpha[i]), -d*np.cos(alpha[i]) - e*np.sin(alpha[i]), d*np.cos(alpha[i]) + e*np.sin(alpha[i]), 0, 0, 0]])
    
    unit_thruster_forces[i,:,:] = np.linalg.pinv(TCM) @ np.eye(6)
   

fig1, axs = plt.subplots(2,3)

subplots = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]])
generalized_forces_txt = np.array(['Fx (Surge)', 'Fy (Sway)', 'Fz (Pitch)', 'Mx (Roll)', 'My (Pitch)', 'Mz (Yaw)'])

for k in range(ndof):
    for j in range(num_thrusters):
        axs[subplots[k,0], subplots[k,1]].plot(np.rad2deg(alpha), unit_thruster_forces[:,j,k], '.-', linewidth=1.5)
        axs[subplots[k,0], subplots[k,1]].grid()
        axs[subplots[k,0], subplots[k,1]].title.set_text("Thrust vs. alpha for unit force " + generalized_forces_txt[k])
        axs[subplots[k,0], subplots[k,1]].set_xlabel("Alpha [deg]")
        axs[subplots[k,0], subplots[k,1]].set_ylabel("Thrust Force")

fig1.legend(['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7'])
fig1.suptitle('Thrust forces vs alpha (b = 7.5 cm)')
plt.show() 


# %% 
# computation of thruster forces while varying b

steps = 20
delta_b = 15 / (steps - 1)
alpha = np.deg2rad(45)

unit_thruster_forces = np.zeros((steps ,num_thrusters, ndof))
b = np.arange(0, 15, delta_b)
b = np.append(b, b[(steps - 1) - 1] + delta_b)

for i in range(steps):
    TCM = np.array([[np.cos(alpha), np.cos(alpha), -np.cos(alpha), -np.cos(alpha), 0, 0, 0],
                    [np.sin(alpha), -np.sin(alpha), np.sin(alpha), -np.sin(alpha), 0, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1],
                    [0, 0, 0, 0, -c, c, 0],
                    [0, 0, 0, 0, -b[i], -b[i], a],
                    [d*np.cos(alpha) + e*np.sin(alpha), -d*np.cos(alpha) - e*np.sin(alpha), -d*np.cos(alpha) - e*np.sin(alpha), d*np.cos(alpha) + e*np.sin(alpha), 0, 0, 0]])
    
    unit_thruster_forces[i,:,:] = np.linalg.pinv(TCM) @ np.eye(6)



fig2, axs = plt.subplots(2,3)

subplots = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]])
generalized_forces_txt = np.array(['Fx (Surge)', 'Fy (Sway)', 'Fz (Pitch)', 'Mx (Roll)', 'My (Pitch)', 'Mz (Yaw)'])

for k in range(ndof):
    for j in range(num_thrusters):
        axs[subplots[k,0], subplots[k,1]].plot(b, unit_thruster_forces[:,j,k], ".-", linewidth=1.5)
        axs[subplots[k,0], subplots[k,1]].grid()
        axs[subplots[k,0], subplots[k,1]].title.set_text("Thrust vs. b for unit force " + generalized_forces_txt[k])
        axs[subplots[k,0], subplots[k,1]].set_xlabel("b [cm]")
        axs[subplots[k,0], subplots[k,1]].set_ylabel("Thrust Force")

fig2.legend(['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7'])
fig2.suptitle('Thrust forces vs b (alpha = 45 deg)')
plt.show() 
   