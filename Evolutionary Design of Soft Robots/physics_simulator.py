import math
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Graphics
"""

# draw and show system, wireframe
def show_bot(masses, springs):

  fig = plt.figure(figsize=(4,4))
  ax = fig.add_subplot(111, projection='3d')
  scale = 0.15
  ax.set_xlim(-1*scale, scale)
  ax.set_ylim(-1*scale, scale)
  ax.set_zlim(0, 2*scale)
  ax.grid(False)

  for m in masses:
    ax.scatter(m.p[0], m.p[1], m.p[2],s=50,c='orange')

  for s in springs:
    v1 = masses[s.m1].p
    v2 = masses[s.m2].p
    ax.plot([v1[0],v2[0]], [v1[1],v2[1]], [v1[2],v2[2]], color='blue')

  plt.show()

"""
Physics Simulator
"""

g = (0,0,-9.81)   # gravity, meters/s^2
dt = 0.0001       # timestep, seconds
T = 0             # time, seconds
kc = 100000       # ground restoration constant


class Mass:

    def __init__(self, m=1, p=[0,0,0], v=[0,0,0], a=[0,0,0]) -> None:
        self.m = np.array(m)    # mass, kg
        self.p = np.array(p)    # position, meters
        self.v = np.array(v)    # velocity, meters/sec
        self.a = np.array(a)    # acceleration, meters/sec^2


class Spring:

  def __init__(self, m1, m2, k=10000, L0=1) -> None:
    self.m1 = m1    # m1,m2 are the indices of the two masses self connects
    self.m2 = m2
    self.k = k      # spring constant, N/meter
    self.L0 = L0    # original resting length, meters


class Structure:

  def __init__(self, vertices):
    self.vertices = vertices
    self.set_masses()
    self.set_springs()
    self.time_elapsed = 0.0

  def set_masses(self):
    self.masses = []
    for vertex in self.vertices:
      m = Mass(0.1, vertex)
      self.masses.append(m)

  def set_springs(self):
    self.springs = []
    mass_indices = np.arange(len(self.masses))
    for m1,m2 in combinations(mass_indices,2):
      p1 = self.masses[m1].p
      p2 = self.masses[m2].p
      L0 = math.sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )
      s = Spring(m1,m2,10000,L0)
      self.springs.append(s)

  def get_mass_positions(self):
    positions = np.zeros((len(self.masses), 3))
    for i,mass in enumerate(self.masses):
      positions[i,:] = mass.p
    return positions
  
  """
  def get_spring_vertices(self)

    for s in springs:
    v1 = masses[s.m1].p
    v2 = masses[s.m2].p
    ax.plot([v1[0],v2[0]], [v1[1],v2[1]], [v1[2],v2[2]], color='blue')
  """
  
  # move structure
  def translate(self, dx, dy, dz):
    for mass in self.masses:
      mass.p[0] += dx
      mass.p[1] += dy
      mass.p[2] += dz
  
  # rotate structure about axis !!! DOESN'T REALLY WORK
  def rotate(self, theta, axis):
    # convert theta to radians
    theta = theta * math.pi /180
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    if axis == 0: # rotation about x axis
      for mass in self.masses:
        mass.p[1] = mass.p[1]*cos_t - mass.p[2]*sin_t
        mass.p[2] = mass.p[1]*sin_t + mass.p[2]*cos_t
    if axis == 1: # rotation about y axis
      for mass in self.masses:
        mass.p[2] = mass.p[2]*cos_t - mass.p[0]*sin_t
        mass.p[0] = mass.p[2]*sin_t + mass.p[0]*cos_t
    if axis == 2: # rotation about z axis
      for mass in self.masses:
        mass.p[0] = mass.p[0]*cos_t - mass.p[1]*sin_t
        mass.p[1] = mass.p[0]*sin_t + mass.p[1]*cos_t

  #def energy(self):

  # update for one time step
  def step(self):

    #global g
    #global dt
    #global T
    #global kc

    # tally spring force on all masses
    Fs = np.zeros((len(self.masses), 3)) # matrix of forces on masses

    for s in self.springs:
      v1 = self.masses[s.m1].p
      v2 = self.masses[s.m2].p
      L = math.sqrt( sum((v1 - v2) * (v1 - v2)) )
      #L = math.sqrt( (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2 )

      mag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
      dir1 = (v2 - v1)/L
      dir2 = -1*dir1

      Fs[s.m1,:] = Fs[s.m1,:] + mag*dir1
      Fs[s.m2,:] = Fs[s.m1,:] + mag*dir2

    for i,mass in enumerate(self.masses):

      # sum all forces on the mass (from springs, gravity, collision with the ground, etc)
      Fg = mass.m * np.array(g) # gravitational force
      Fe = np.array([0,0,0])    # external force
      Fc = np.array([0,0,0])    # restoration force from ground
      if mass.p[2] < 0:
        Fc[2] = -1*kc*mass.p[2]

      F = Fs[i,:] + Fg + Fe + Fc # total force on mass
      
      mass.a = F/mass.m
      mass.v = mass.v + mass.a * dt
      mass.p = mass.p + mass.v * dt
    
    self.time_elapsed += dt

tetra_vertices = [[0.1,0.1,0.1],[0.1,-0.1,-0.1],[-0.1,0.1,-0.1],[-0.1,-0.1,0.1]]
tetra_vertices = np.array(tetra_vertices)*0.5
tetrahedron = Structure(tetra_vertices)
tetrahedron.translate(0,0,0.15)

pos = tetrahedron.get_mass_positions()
print(pos)
show_bot(tetrahedron.masses, tetrahedron.springs)