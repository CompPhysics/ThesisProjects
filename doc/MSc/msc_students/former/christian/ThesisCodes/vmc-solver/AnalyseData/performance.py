import numpy as np
import matplotlib.pyplot as plt
import sys

NumProcs = [1, 2, 4, 8, 16, 32, 64]
O0 = [545.218, 276.395, 138.610, 69.4340, 35.0710, 17.7141, 8.91754]
O3 = [158.773, 80.4636, 40.0966, 20.2077, 10.1289, 5.12122, 2.57613]
O0 = np.asarray(O0)
O3 = np.asarray(O3)

plt.plot(NumProcs, O0[0]/O0, "-o")
plt.xlabel("Number of Processors")
plt.ylabel("Speed-up Factor")
plt.title("Speed-up Factor from Parallelization")
plt.figure()
plt.plot(NumProcs, O0, "-o", label="O0 (No Vectorization)")
plt.plot(NumProcs, O3, "-o", label= "O3 Vectorization")
plt.xlabel("Number of Processors")
plt.ylabel("Computation Time [s]")
plt.title("Computation Time vs. Number of Processors")
plt.legend()

plt.show()

SpeedUpO0 = O0[0]/O0
SpeedUpO3 = O0[0]/O3
print SpeedUpO0
print SpeedUpO3
