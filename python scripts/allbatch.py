import os
import math
import numpy as np
import time
from multiprocessing import Pool

#--CK preprocessing--
#for sof in [1,2,5,10,15,20,25,30,35,40]:
#  for cov in [0.2,0.4,0.6,0.8]:
#    os.system('sbatch batchs1.sh 1 '+str(cov)+' '+str(sof))

#exit()

# create custom site investigation input files
if False:
    tests = ['dct']  #['dct', 'cts', 'CPT', 'SPT', 'DMT', 'TT']
    depths = [5,10,15,20,25,30] # [5,10,15,20,25,30,35,40]
    no_bh = [1,4,9,16,25,36,49,64,81,100]
    areas = np.array([ 5. ,  7.5, 10. , 12.5, 15. , 17.5, 20. , 22.5, 25. , 27.5, 30. ,
           32.5, 35. , 37.5, 40. , 42.5, 45. , 47.5, 50. , 52.5, 55. , 57.5,
           60. , 62.5, 65. , 67.5, 70. , 72.5, 75. , 77.5, 80. , 82.5, 85. ,
           87.5, 90. ])

    offsets = [(90-a)/2 for a in areas]

    si = open('input/si.txt','w')
    si.write(str(len(tests) * len(depths) * len(no_bh) * len(areas)) + '\n')
    si.write(str(max(no_bh)) + '\n')
    si.write('\n\n\n\n')

    x_coords = open('input/si_Xcoords.txt', 'w')
    x_coords.write('\n')
    y_coords = open('input/si_Ycoords.txt', 'w')
    y_coords.write('\n')

    j = 0
    for test in tests:
        for depth in depths:
            for a, o in zip(areas, offsets):
                for i, nbh in enumerate(no_bh):
                    spacing = a / (math.sqrt(nbh) - 1)
                    si.write(f'{nbh} {test} SD {depth}\n')
                    coords = [o + bh * spacing for bh in range(int(math.sqrt(nbh)))]
                    if len(coords) == 1:
                        coords = [45]
                    coords = np.round(coords, 1)
                    x_coords.write(' '.join(np.repeat(coords,int(math.sqrt(nbh))).astype(str)) +'\n')
                    y_coords.write(' '.join(np.tile(coords, int(math.sqrt(nbh))).astype(str)) + '\n')
                    j+=1

    si.close()
    x_coords.close()
    y_coords.close()

    
detset = np.array([1.028974  , 0.5592995 , 0.4034411 , 0.320887  , 0.2685457 ,
       0.2319366 , 0.2046802 , 0.1834855 , 0.166467  , 0.1524599 ,
       0.140704  , 0.1306791 , 0.1220162 , 0.1144464 , 0.1077667 ,
       0.1018228 , 0.09649403, 0.09168521, 0.08731996, 0.08333601,
       0.07968278])

def redcoefs(cov):
	"""
	This function returns reduction coefficients for reduction methods.
	"""
	m = -math.log(1+cov**2)/2
	v = math.sqrt(math.log(1+cov**2))
	gmean = math.exp(m)
	
	return 1.0,gmean,gmean**2,math.exp(m-v*0.675),gmean/math.exp(v)


pc = list()

pc.append([
    np.array([[1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1]]),
    1600, 80])

pc.append([
    np.array([[0, 0, 0, 0, 0, 0],
              [0, 1, 1, 1, 1, 0],
              [0, 1, 1, 1, 1, 0],
              [0, 1, 1, 1, 1, 0],
              [0, 1, 1, 1, 1, 0],
              [0, 0, 0, 0, 0, 0]]),
    900, 60])


pc.append([
    np.array([[0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 1, 1, 0, 0],
              [0, 0, 1, 1, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    100, 20])
pc.append([
    np.array([[0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0]]),
    900, 60])

pc.append([
    np.array([[1, 0, 0, 0, 0, 1],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 1]]),
    1600, 80])

pc2 = list()
pc2.append([
    np.array([[1, 1, 1, 1, 1, 0],
              [1, 1, 1, 1, 1, 0],
              [1, 1, 1, 1, 1, 0],
              [1, 1, 1, 1, 1, 0],
              [1, 1, 1, 1, 1, 0],
              [0, 0, 0, 0, 0, 0]]),
    1600, 80])
pc2.append([
    np.array([[1, 0, 1, 0, 1, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 1, 0, 1, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 1, 0, 1, 0],
              [0, 0, 0, 0, 0, 0]]),
    1600, 80])
pc2.append([
    np.array([[1, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0]]),
    1600, 80])
pc2.append([
    np.array([[1, 1, 1, 1, 0, 0],
              [1, 1, 1, 1, 0, 0],
              [1, 1, 1, 1, 0, 0],
              [1, 1, 1, 1, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    900, 60])
pc2.append([
    np.array([[1, 1, 1, 0, 0, 0],
              [1, 1, 1, 0, 0, 0],
              [1, 1, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    400, 40])
pc2.append([
    np.array([[1, 1, 0, 0, 0, 0],
              [1, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    100, 20])
pc2.append([
    np.array([[1, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    400, 40])
pc2.append([
    np.array([[1, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [1, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]]),
    900, 60])

def detcurve(smean):
    dettemp = detset * bweight / ((smean * npl * covfactor))
    des = np.interp(destol,dettemp[::-1],range(21)[::-1])
    return des


smean = 20
cov = 0.8

test = 2
red = 5
depth = 40
vsof = 10.01


#for sof in [1,2,5,10,15,20,30]:
#	for cov in [0.2,0.4,0.8]:
#		args = [1.0,cov,sof]
#		args = [str(i) for i in args]
#		args = ' '.join(args)
#		os.system('sbatch batchs1.sh '+args)

#exit()


# todo: create function that checks whether the CK and/or SI files exist, and only examine the missing cases
# check for pile 1?
# only check SI since CK is a prerequisite anyway?
# add an if-statement to the below code to simply skip adding the arguments to the list if it's present?
# have an option to ignore the file check if desired


arg_strings = []
nbh = 1
for nbh in [1]: #[1,2,3,4,5,6,8,9,12,16,25][::-1]:
    for b in range(len(pc)):
        for red in [5]:
            for cov in [0.4]: # [0.2,0.4,0.8]:
                for hsof in [50,80]: #[5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]: # [10,15,20]:# [1,2,5,10,15,20,30]:
             
                    bweight = 10 * 8 * pc[b][1]
                    npl = np.sum(pc[b][0]==1)
                    destol = 2.5 * pc[b][2] / (2 * (math.sqrt(npl) - 1))
                    covfactor = redcoefs(cov)[red-1]
        #print(covfactor)	 
                    dettemp = detset * bweight / ((smean * npl * covfactor))
                    for smean in np.linspace(200,0.1,2000):
                        if(detcurve(smean)>=5): break
                    args = [smean,cov,hsof,vsof,1,1,1,1,1] + pc[b][0].ravel().tolist() + [pc[b][1],pc[b][2]] #,test,red,depth,nbh]
                    args = [smean, cov, 1, hsof, vsof, 1, 1, 1, 1] + pc[b][0].ravel().tolist() + [pc[b][1], pc[b][2]]
#                    args = [smean,cov,sof,1,1,1,1,1] + pc[b][0].ravel().tolist() + [pc[b][1],160,test,red,depth,nbh]
                    args = [str(i) for i in args]
                    args = ' '.join(args)
                    #print('sbatch batchs1.sh '+args)
                    #os.system('sbatch batchs1.sh '+args)
                    print('soil.exe ' + args)
                    arg_strings.append('soil.exe ' + args)
                    #start = time.time()
                    #os.system('soil.exe ' + args)
                    #print(time.time() - start)

if __name__ == '__main__':
    # do parallel processing
    with Pool(5) as p:
        p.map(os.system, arg_strings)
    
    
