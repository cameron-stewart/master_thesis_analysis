import csv
import numpy as np

leng = 10000
betas = np.array([1,3,5,7,9])
nruns = betas.size
vel = np.zeros((leng,nruns))

j=0
for beta in betas:
    with open('vel' + str(beta) + '.dat', 'rb') as f:
        reader = csv.reader(f)
        i=0
        for row in reader:
            for col in row:
                if col != row[-1]:
                    if i < leng:
                        vel[i,j] = float(np.array(row)[i])
                    i+=1
    j +=1
    f.close()
    np.save('numpy.save',vel)
