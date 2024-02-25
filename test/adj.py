import numpy as np

m1 = np.matrix([[0.1, 0.2, 0.3, 0.4], 
    [0.5, 0, 0, 0], 
    [0.9, 0.10, 0, 0], 
    [0, 0.14, 0.15, 1]])
print(np.linalg.inv(m1))


