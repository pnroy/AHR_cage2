import numpy as np

sintl = [np.sin(x) for x in range(5)]

sinta = np.fromiter(sintl, np.float)
print sinta
