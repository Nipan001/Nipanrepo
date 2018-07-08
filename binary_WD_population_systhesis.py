
# coding: utf-8

# In[1]:


"""
Evolves a stellar binary and reports the mass of each star during the evolution.

Shows the type of each star as they change.
"""


#program by Rahul Kashyap and Nipan Goswami
# june,2018,ICTS


# In[2]:



from __future__ import print_function
from amuse.lab import *
import numpy
from matplotlib import pyplot #for making plots
import itertools #for making initial binaries array


# In[3]:


# main function of coding
def evolve_binary(mass_of_star1, mass_of_star2, orbital_period, eccentricity):
    code = SeBa()

    stars =  Particles(2)
    stars[0].mass = mass_of_star1
    stars[1].mass = mass_of_star2


    mu = stars.mass.sum() * constants.G
    semi_major_axis = (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)


    binaries =  Particles(1)

    binary = binaries[0]
    binary.semi_major_axis = semi_major_axis
    binary.eccentricity = eccentricity
    binary.child1 = stars[0]
    binary.child2 = stars[1]

    # we add the single stars first, as the binaries will refer to these
    code.particles.add_particles(stars)
    code.binaries.add_particles(binaries)

    from_seba_to_model = code.particles.new_channel_to(stars)
    from_seba_to_model.copy()

    from_seba_to_model_binaries = code.binaries.new_channel_to(binaries)
    from_seba_to_model_binaries.copy()

    previous_type_child1 = binary.child1.stellar_type
    previous_type_child2 = binary.child2.stellar_type

    results = []
    current_time = 0 | units.Myr
    while current_time < (1000 | units.Myr):
        code.update_time_steps()
        # The next line appears a bit weird, but saves time for this simple test.
        deltat = max(1.0*code.binaries[0].time_step, 0.1 | units.Myr)
        current_time = current_time + deltat
        code.evolve_model(current_time)
        from_seba_to_model.copy()
        from_seba_to_model_binaries.copy()

        if not binary.child1.stellar_type == previous_type_child1:
            print(binary.age, "Child 1, change of stellar type", previous_type_child1, ' -> ',binary.child1.stellar_type)
            previous_type_child1 = binary.child1.stellar_type
        if not binary.child2.stellar_type == previous_type_child2:
            print(binary.age, "Child 2, change of stellar type", previous_type_child2, ' -> ',binary.child2.stellar_type)
            previous_type_child2 = binary.child2.stellar_type
        results.append((binary.age, binary.child1.mass, binary.child1.stellar_type, binary.child2.mass, binary.child2.stellar_type))



    code.stop()
    return results


# In[5]:


n_stars=60
c=n_stars*1.3/(.5**(-1.3)-(10**(-1.3)))
c


# In[6]:


number_of_star_in_intervals=[]
df=numpy.linspace(1.0,19.0,19)
for i in df:
    c_1=(c/(-1.3))*((i/2+0.5)**(-1.3)-(i/2)**(-1.3))
    number_of_star_in_intervals.append(c_1)

print(sum(number_of_star_in_intervals))
#print(number_of_star_in_intervals)


# In[7]:


making_no_of_star_integer=[]
for i in range (0,19):
    m_1=int(round(number_of_star_in_intervals[i],0))
    making_no_of_star_integer.append(m_1)
print(making_no_of_star_integer)
print(sum(making_no_of_star_integer))


# In[8]:


arr=[]
for i in range(1,20):
    arr.append(numpy.linspace((i+.0)/2,((i+.0)/2)+0.5,making_no_of_star_integer[i-1]))
#print(arr)


# In[9]:


initial_mass_of_star=[]
for i in range(0,19):
    for j in range(making_no_of_star_integer[i]):
        initial_mass_of_star.append(arr[i][j])
#print(initial_mass_of_star)
#print(len(initial_mass_of_star))


# In[10]:



# x axis values
x = numpy.linspace(.75,9.75,19)
# corresponding y axis values
y = making_no_of_star_integer
# plotting the points 
pyplot.semilogx(x, y)
 
# naming the x axis
pyplot.xlabel('inital mass of star')
# naming the y axis
pyplot.ylabel('no of stars')
pyplot.show()


# In[11]:


#creating initial binray population
mass_of_star1 =initial_mass_of_star  
mass_of_star2 = initial_mass_of_star
orbital_period=[10,100,1000,10000]
eccentricity= numpy.linspace(.95,.25,5)

binary_tuple=[mass_of_star1,mass_of_star2,orbital_period,eccentricity]
#print(binary_tuple)


# In[12]:


#for making all the possible combination of binary stars
binary_population=[]
for element in itertools.product(*binary_tuple):
    #print(element)
    binary_population.append(element)
#print(binary_population)
#len(binary_population)


# In[ ]:


# collecting the final stage datas of stellar evolution
time=[]
mass1_final=[]
mass2_final=[]
type1_final=[]
type2_final=[]
for i, binary in enumerate(binary_population):
   print(i, "Binary mass1, mass2, period, eccentricity are:",binary)
   mass1=binary[0]; mass2=binary[1]; period=binary[2]; eccentricity=binary[3]
   table = evolve_binary( mass1 | units.MSun, mass2 | units.MSun, period | units.day, eccentricity)

   type1_final.append(table[-1][2])
   type2_final.append(table[-1][4])
   mass1_final.append(table[-1][1])
   mass2_final.append(table[-1][3])
   time.append(table[-1][0])


# In[7]:


#print(mass1_final)
#print(mass2_final)
#print(type1_final)
#print(type2_final)
#print(time)


# In[8]:


import astropy
#print(astropy.to_value(table[-1][1]))
#var=table[-1][1]
#print(var.value())
final_mass_of_star1=[]
final_mass_of_star2=[]
time_to_final_state=[]
for mass1 in mass1_final: final_mass_of_star1.append(mass1.value_in(units.MSun))
#print(final_mass_of_star1)
pyplot.hist(final_mass_of_star1,10)
pyplot.title('final_mass_1')
pyplot.xlabel('final mass of star 1')
pyplot.ylabel('no of star')
pyplot.show()
for mass2 in mass2_final: final_mass_of_star1.append(mass2.value_in(units.MSun))
#print(final_mass_of_star1)
pyplot.hist(final_mass_of_star1,10)
pyplot.title('final_mass_2')
pyplot.xlabel('final mass of star 2')
pyplot.ylabel('no of star')
pyplot.show()
for age in time: time_to_final_state.append(age.value_in(units.Myr))
#print(time_to_final_state)
pyplot.hist(time_to_final_state,10)
pyplot.title('evolution_time')
pyplot.xlabel('time')
pyplot.ylabel('no of star')
pyplot.show()

