from vpython import *

g = vector(0,0,-9.81)   # gravity, meters/s^2
dt = 0.0001 #0.0001     # timestep, seconds
T = 0                   # time, seconds
kc = 100000             # ground restoration constant
ks = 1000                # spring constant


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
        #self.img = curve(pos=[masses[self.m1].p, masses[self.m2].p], radius=0.003, color=color.blue)
    
    def update_img(self,masses):
        self.img.pos = masses[self.m1].p
        self.img.axis = masses[self.m2].p - masses[self.m1].p
        #self.img.pos = [masses[self.m1].p, masses[self.m2].p]


class Structure:

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
        ground_vertices = [vector(0.,0.4,0.), vector(0.,-0.4,0.), vector(-0.4,0.,0.), vector(0.4,0.,0.)]

    
        #mass_indices = arange(len(self.masses))
        #for m1,m2 in combinations(mass_indices,2):
        m1 = 0
        while m1 < len(self.masses)-1:
            m2 = m1 + 1
            while m2 < len(self.masses):
                p1 = self.masses[m1].p
                p2 = self.masses[m2].p
                if p1 not in ground_vertices or p2 not in ground_vertices:
                    L0 = sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2 )
                    #if L0 > 0.1:
                        #L0 = L0 * 0.9
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
            mass.p = mass.p + mass.v * dt
            

        #self.time_elapsed += dt


scene.camera.pos = vector(0,-0.6,0.5)
scene.camera.axis = vector(0,0,0) - scene.camera.pos
floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)


"""
tetra_vertices = [vector(0.1,0.1,0.1),vector(0.1,-0.1,-0.1),vector(-0.1,0.1,-0.1),vector(-0.1,-0.1,0.1)]
tv = []
for v in tetra_vertices:
    tv.append(v*0.5)
tetra = Structure(tv)
tetra.translate(0,0,0.25)
"""
#cube_vertices = [vector(0.,0.,0.2), vector(0.,0.1,0.2), vector(0.1,0.,0.2), vector(0.1,0.1,0.2), vector(0.,0.,0.3), vector(0.,0.1,0.3), vector(0.1,0.,0.3), vector(0.1,0.1,0.3)]
#ball = Structure(cube_vertices)

ball = Structure([vector(0.,0.,0.2)])

i = 0
while i < 1000:
    rate(1000)
    ball.step()

    if i % 10 == 0:
        ball.update_img()

    if i % 100 == 0:
        print(ball.get_energy())
    
    T = T + dt
    i += 1