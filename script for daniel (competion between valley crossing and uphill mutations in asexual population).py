# packages
import numpy as np
import csv
import sys
import time
from numpy import zeros
from matplotlib.pyplot import (legend, plot, show, 
                               title, xlabel, ylabel)

# global variables
pop_size = 1.e6 # population size 
num_gen = 1000000   # number of generations
muta = .001   # A mutation rate 
mutb = 0#.00001   # B mutation rate
mutc =.00000001   # C mutation rate
sa = 0.08      # Fitness effect of a mutation A 
sb = 0.01      # Fitness effect of a mutation B
sc = sb     # Fitness effect of a mutation C
count_x = 0     # number of zero zy slices of the simulation grid starting at x = 0  
numexp = 1;
simulation_grid =  np.zeros((30,100,200)) # grid on which the simulation occurs
num_Amut = np.zeros((num_gen)) # array to save number of available A-mutations
num_Amut_hist = np.zeros((num_gen,30)) # array to save a histogram of A-mutations
increment = 1;

# get fitness of a grid point.
def getFitness(ka,kb,kc):
    return (kb*sb+kc*sc-ka*sa)
# incoming -outgoing (x,y,z)
def freq_lineage(simulation_grid,x,y,z):
    const = simulation_grid[x,y,z]*(1+getFitness(x,y,z))*(1-muta-mutb-mutc*x)
    if x>0 : const_x=simulation_grid[x-1,y,z]*(1+getFitness(x-1,y,z))*muta
    else: const_x = 0
    if y>0:  const_y=simulation_grid[x,y-1,z]*(1+getFitness(x,y-1,z))*mutb
    else: const_y=0
    if z>0 and x+1<simulation_grid.shape[0]: const_z=simulation_grid[x+1,y,z-1]*(1+getFitness(x+1,y,z-1))*(x+1)*mutc
    else: const_z=0
    return const+const_x+const_y+const_z
       
# Compute probablities of offspring given the current size of various lineages
# Note that x=ka y=kb z=kc in the implementation of this function
def probability_offspring(simulation_grid,xl,xh,yl,yh,zl,zh):
    sim_grid = np.zeros((xh+1-xl)*(yh+1-yl)*(zh+1-zl))
    i=0
    for x in range(xl,xh+1):
         for y in range(yl,yh+1):
              for z in range(zl,zh+1):
                  sim_grid[i]   =  freq_lineage(simulation_grid,x,y,z) 
                  i=i+1
    sim_grid = sim_grid/sum(sim_grid)
    return np.random.multinomial(pop_size, sim_grid) 
    
# Update the simulation grid to reflect a time step
def update_grid(simulation_grid,xl,xh,yl,yh,zl,zh):
    i=0  
    offspring_prob = probability_offspring(simulation_grid,xl,xh,yl,yh,zl,zh)
    for x in range(xl,xh+1):
        for y in range(yl,yh+1):
            for z in range(zl,zh+1):
                simulation_grid[x][y][z] = offspring_prob[i]
                i=i+1
    return    
    
# make sure that the simulation goes out of bounds by increasing the array size appropriately if need be.
def boundConditions(xmax,ymax,zmax):
    global simulation_grid,increment
    if xmax+1==simulation_grid.shape[0]:
       simulation_grid = np.pad(simulation_grid,((0,10*increment),(0,0),(0,0)),mode='constant', constant_values=0)
    if ymax+1==simulation_grid.shape[1]:
        simulation_grid = np.pad(simulation_grid,((0,0),(0,10*increment),(0,0)),mode='constant', constant_values=0)
    if zmax+1==simulation_grid.shape[2]:
        simulation_grid = np.pad(simulation_grid,((0,0),(0,0),(0,5*increment)),mode='constant', constant_values=0)
    return   


#main function which runs the code
def main():
    global mutc,simulation_grid,num_Amut
    vals = np.zeros((numexp,3))
    name = sys.argv[1]
    #mutcArr = np.array([1.e-7,10**(-6.67),10**(-6.33),10**(-6),10**(-5.67),10**(-5.33),10**(-5),10**(-4.67),10**(-4.33),10**(-4)])
    mutcArr = np.array([0,1.e-6,10**(-5),10**(-4),10**(-3)])
    sizes = np.array([[25,180,15],[25,180,15],[25,180,30],[25,165,45],[25,140,80],[24,140,80],[25,90,170],[25,70,250],[25,35,350],[25,35,450]])
    mutc = mutcArr[int(sys.argv[1])]
    for k in range(0,numexp):
        find_time = 0
        simulation_grid = np.zeros(sizes[int(sys.argv[1])])
        simulation_grid[0][0][0]=pop_size
        for i in range(0,num_gen):
            #define a time start variable to be used to compute actual run time of simulation in seconds.
            start1 = time.time()
            #obtain the coords(genomic composition) in the grid where there is a nonzero number of individuals. 
            nz = np.nonzero(simulation_grid)
            # Find the minimum coords and the minimum,maximum coords (smallest/largest number of mutations of each type(A,B,C))
            vectm = map(min,nz)
            vectx = map(max,nz)
            #computing boundary conditions
            boundConditions(vectx[0],vectx[1],vectx[2]) 
            #updating the grid 
            update_grid(simulation_grid,vectm[0],vectx[0]+1,vectm[1],vectx[1]+1,vectm[2],vectx[2]+1)
            # Computing the number of available A-mutations , histograms of lineages with x A-mutations.
            for x in range(vectm[0],vectx[0]+1):
                for y in range(vectm[1],vectx[1]+1):
                    for z in range(vectm[2],vectx[2]+1):
                        num_Amut[i] += simulation_grid[x,y,z]*x
                        num_Amut_hist[i][x]+=simulation_grid[x,y,z];
            #compute end of simulation time
            end1 = time.time()
            #compute total run time
            find_time = find_time+ end1 - start1
        #obtain the coords(genomic composition) in the grid where there is a nonzero number of individuals. 
        nz = np.nonzero(simulation_grid)
        # Find the minimum coords and the minimum,maximum coords (smallest/largest number of mutations of each type(A,B,C))
        vect = map(min,nz)
        vals[k][0] = vect[0]  # A mutations shared fixed in the population
        vals[k][1] = vect[1]  # B mutations shared fixed in the population
        vals[k][2] = vect[2]  # C mutations shared fixed in the population
	print vect[0],vect[1],vect[2]
    #create an index array (to save time in generations)
    index = zeros((num_gen))
    #Save the histogram of A-mutations in the population in a comma seperated values file.
    np.savetxt("histe6.csv", num_Amut_hist, delimiter=",")
    for i in range(0,num_gen):
        index[i] = i+1 # populate the index array (to be used as the time-axis in the plot)
    # write the number of A mutation as a function of time.
    with open(name+'6,Ub=0,Uc=0'+'.csv', 'w') as csvfile: #open a csv file
        writer = csv.writer(csvfile, delimiter = ',', lineterminator='\n')
        fieldnames = ['time','A muts'] # two columns (time, number of available A-mutations in the population)
        writer.writerow(fieldnames)  # write the field names in the first row of the csv file
        for k in range(0,num_gen):
            writer.writerow([k,num_Amut[k]]) # write the (time, number of A-mutations) in the k+1 row 
    # plotting data
    plot(index, num_Amut,'r-',label ='number of A mutations')
    legend(loc='upper center', shadow=True)
    xlabel('Time (generations)')
    ylabel('Available A mutations')
    title('Fluctuations',fontsize = 'large')
    show()
    


start = time.time()
main()
end = time.time()
t =  (end - start)/60 # time in minutes
print 'total time', t


