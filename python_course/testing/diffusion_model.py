from numpy import *

def energy(particles, D):

  assert(isinstance(particles, list));

  energy = 0;
  for i in range(len(particles)):
    energy += particles[i]*(particles[i]-1);


  return double(energy) * D/2;


#print diffusion_model([1,3,7,2,9], 5.0);
