import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  cuba = np.genfromtxt("CUBACORES.txt", delimiter=",", skip_header=4)
  omp = np.genfromtxt("OMP.txt", delimiter=",", skip_header=4)
  ref_time = omp[-1,1]
  
  fig,ax = plt.subplots()
  ax.plot(cuba[:,0], cuba[:,1], color='b', label="CUBACORES", marker='o')
  ax.plot(omp[:,0], omp[:,1], color='r', label="OMP", marker='s')
  max_core = int(np.max(omp[:,0]))
  ax.plot([x+1 for x in range(max_core)],[ref_time/(x+1) for x in range(max_core)], linestyle='--', color='k', label='Ideal')
  
  ax.set_xlabel("Cores")
  ax.set_ylabel("Time (s)")
  ax.grid()
  ax.set_ylim(8,ref_time*1.1)
  ax.set_xlim([1,max_core])
  ax.set_yscale("log", base=2)
  ax.set_xscale("log", base=2)
  ax.set_yticks([int(2**x) for x in range(3,10)])
  ax.set_yticklabels([str(int(x)) for x in ax.get_yticks()])
  ax.set_xticks([int(2**x) for x in range(0,6)])
  ax.set_xticklabels([str(int(x)) for x in ax.get_xticks()])
  ax.legend()
  
  plt.tight_layout()
  plt.show()

