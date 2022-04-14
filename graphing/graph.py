import sys
from matplotlib import pyplot as plt
len_args = len(sys.argv)

if len_args == 3:
    with open(sys.argv[1], "r") as f:
        file = f.read().split(",")[:-1]
    n = []
    gflops = []
    mem_band = []
    for i in range(len(file)):
        if i % 3 == 0:
            n.append(float(file[i]))
        elif i % 3 == 1:
            gflops.append(float(file[i]))
        else:
            mem_band.append(float(file[i]))
    plt.subplot(2, 1, 1)
    plt.plot(n, gflops)
    plt.title("N vs GFLOPS")
    plt.subplot(2, 1, 2)
    plt.plot(n, mem_band)
    plt.title("N vs Mem Band.")
    plt.savefig(sys.argv[2])
