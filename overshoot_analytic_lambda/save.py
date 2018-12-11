import csv
import numpy as np

leng = 12000
lambdas = np.array([1000,3000,5000,7000,9000])
nruns = lambdas.size
vel = np.zeros((leng,nruns))

j=0
for x in lambdas:
    with open('vel' + str(x) + '.dat', 'rb') as f:
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

    np.save("numpy.save",vel)
