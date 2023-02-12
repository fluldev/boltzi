import numpy as np
import matplotlib.pyplot as plt 

if __name__ == "__main__":
    with open("build/temperature_39.9.boltzi_artifact", "r") as f:
        data = np.array([[*list(map(float, k.removesuffix("\n").split(","))), float(v.removesuffix("\n"))] for k,v in zip(f,f)])

    plt.plot(data[:, 0], data[:, 3])
    plt.plot([-1.1,1.1], [12, 10], c="red", ls="--")
    plt.xlabel("$z / m$")
    plt.ylabel("T / K")
    plt.show()
