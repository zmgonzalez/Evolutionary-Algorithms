# Physics Simulator Version 3

from vpython import *

g = vector(0,0,-9.81)   # gravity, meters/s^2
dt = 0.0001             # timestep, seconds
T = 0                   # time, seconds
kc = 100000             # ground restoration constant
ks = 1000               # default spring constant
dampening = True        # add dampening


class Mass:

    def __init__(self, m=1, p=vector(0,0,0), v=vector(0,0,0), a=vector(0,0,0)) -> None:
        self.m = m      # mass, kg
        self.p = p      # position, meters
        self.v = v      # velocity, meters/sec
        self.a = a      # acceleration, meters/sec^2
        
        self.img = sphere(pos=vector(self.p), radius=0.01, color=color.orange)
        
    def update_img(self):
        self.img.pos = self.p

class Spring:

    def __init__(self, m1, m2, masses, k=10000, L0=1) -> None:
        self.m1 = m1    # m1,m2 are the indices of the two masses self connects
        self.m2 = m2
        self.k = k      # spring constant, N/meter
        self.L0 = L0    # original resting length, meters
        #self.og_L0 = L0
        
        ax = masses[self.m2].p - masses[self.m1].p
        self.img = cylinder(pos=masses[self.m1].p, axis=ax, radius=0.005, color=color.blue)
    
    def update_img(self,masses):
        self.img.pos = masses[self.m1].p
        self.img.axis = masses[self.m2].p - masses[self.m1].p

class Robot:

    def __init__(self, vertices):
        self.vertices = vertices
        self.set_masses()
        self.set_springs()

    def set_masses(self):
        self.masses = []
        for v in self.vertices:
            m = Mass(0.1, v)
            self.masses.append(m)

    def set_springs(self):
        self.springs = []
        m1 = 0
        while m1 < len(self.masses)-1:
            m2 = m1 + 1
            while m2 < len(self.masses):
                p1 = self.masses[m1].p
                p2 = self.masses[m2].p
                L0 = sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2 )
                s = Spring(m1,m2,self.masses,ks,L0)
                self.springs.append(s)
              
                m2 = m2+1
            m1 = m1+1

    def update_img(self):
        for mass in self.masses:
            mass.update_img()
        for spring in self.springs:
            spring.update_img(self.masses)
    
    def get_mass_positions(self):
        positions = []
        for mass in self.masses:
          positions.append(mass.p)
        return positions     
  
    # move structure
    def translate(self, dx, dy, dz):
        for mass in self.masses:
            mass.p.x += dx
            mass.p.y += dy
            mass.p.z += dz
      
        self.update_img()

    def get_energy(self):
        potential = 0
        kinetic = 0

        for mass in self.masses:
            potential += abs(mass.m * g.z * mass.p.z)
            kinetic += (1/2) * mass.m * (mass.v.x**2 + mass.v.y**2 + mass.v.z**2)

        for spring in self.springs:
            v1 = self.masses[spring.m1].p
            v2 = self.masses[spring.m2].p
            L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )
            potential += (1/2) * spring.k * (spring.L0 - L)**2

        total = potential + kinetic
        return T, potential, kinetic, total

    # update for one time step
    def step(self):

        for i,mass in enumerate(self.masses):
            
            # sum all forces on the mass (from springs, gravity, collision with the ground, etc)
            Fs = vector(0,0,0)      # spring force
            for s in self.springs:
                if s.m1 == i:
                    v1 = self.masses[s.m1].p
                    v2 = self.masses[s.m2].p
                    L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )

                    fmag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
                    dir = (v2 - v1)/L
                    Fs = Fs + fmag*dir
                elif s.m2 == i:
                    v1 = self.masses[s.m1].p
                    v2 = self.masses[s.m2].p
                    L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )

                    fmag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
                    dir = (v1 - v2)/L
                    Fs = Fs + fmag*dir

            Fg = mass.m * g     # gravitational force
            Fe = vector(0,0,0)  # external force
            Fc = vector(0,0,0)  # restoration force from ground
            if mass.p.z < 0:
                Fc.z = -1*kc*mass.p.z

            F = Fs + Fg + Fe + Fc # total force on mass
            
            mass.a = F/mass.m
            mass.v = mass.v + mass.a * dt
            if dampening == True:
                mass.v = mass.v * 0.9999
            mass.p = mass.p + mass.v * dt

scene.camera.pos = vector(0,-0.6,0.5)
scene.camera.axis = vector(0,0,0) - scene.camera.pos
floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)

ball = Robot([vector(0.,0.,0.2)])

i = 0
#while i < 5000000:
while T < 5:
    rate(10000) #100000)
    ball.step()

    #if i % 5000 == 0:
    if T % 0.005:
        ball.update_img()

    if i % 500 == 0:
    #if T % 0.05 == 0:
        print(ball.get_energy())
    
    T = T + dt
    i += 1

print(i)