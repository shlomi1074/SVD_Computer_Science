import time
import numpy as np
from SVD import svd

if __name__ == '__main__':
    start_time = time.time()
    print("Load matrix from file\n")
    A = np.loadtxt("data_8.txt", dtype=np.float64, delimiter=',')  # read data to numpy array
    load_run_time = time.time() - start_time
    start_time = time.time()
    print("Matrix loaded\n")

    print("Calculate SVD\n")
    U, S, V = svd(A)
    svd_run_time = time.time() - start_time
    print("Finished SVD\n")

    mulRes = np.matmul(np.matmul(U, S), V.T)

    # Write results to file
    f = open("output.txt", "w")
    f.write("\nA=\n" + str(A))
    f.write("\nV=\n" + str(V))
    f.write("\nS=\n" + str(S))
    f.write("\nU=\n" + str(U))
    f.write("\nUSVt=\n" + str(mulRes))
    result = np.around(mulRes, decimals=4) == np.around(A, decimals=4)
    f.write("\nSVD Results: " + str(result.all()))
    f.write("\nPerformance Results:\n")
    f.write(f"\nLoad File Duration: {str(load_run_time)}\n")
    f.write(f"\nSVD Duration: {str(svd_run_time)}\n")
    print(f"Result {result.all()}")
