# Physics Simulator Version 5

from os import supports_dir_fd
from numpy import true_divide
from vpython import *

g = vector(0,0,-9.81)   # gravity, meters/s^2
dt = 0.0001             # timestep, seconds
#T = 0                   # time, seconds
kc = 100000             # ground restoration constant
ks = 1000               # default spring constant
dampening = True        # add dampening
w = 30                  # global frequency


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

    def __init__(self, m1, m2, k=10000, L0=1, a=1, b=.005, c=0) -> None:
        self.m1 = m1    # m1,m2 are the two masses self connects
        self.m2 = m2
        self.k = k      # spring constant, N/meter
        self.L0 = L0    # original resting length, meters
        self.a = a      # spring parameters
        self.b = b
        self.c = c
        
        ax = self.m2.p - self.m1.p
        self.img = cylinder(pos=self.m1.p, axis=ax, radius=0.005, color=color.blue)
    
    def update_img(self,masses):
        self.img.pos = self.m1.p
        self.img.axis = self.m2.p - self.m1.p

    def breathe(self, t):
        self.L0 = self.a + self.b*sin(w*t + self.c)

class Cube:

    def __init__(self, connecting_cube=None, connecting_side=None) -> None:
     
        self.cc = connecting_cube

        if connecting_cube == None:
            self.set_first()
        else:
            if connecting_side == 2:
                self.set_2()
            if connecting_side == 4:
                self.set_4()
        #self.set_masses()
        #self.set_springs()
        #self.initialize_img()

    def set_first(self):
        self.m123 = Mass(0.1, vector(0.,0.,0.))
        self.m124 = Mass(0.1, vector(0.,0.1,0.))
        self.m135 = Mass(0.1, vector(0.1,0.,0.))
        self.m145 = Mass(0.1, vector(0.1,0.1,0.))
        self.m236 = Mass(0.1, vector(0.,0.,0.1))
        self.m246 = Mass(0.1, vector(0.,0.1,0.1))
        self.m356 = Mass(0.1, vector(0.1,0.,0.1))
        self.m456 = Mass(0.1, vector(0.1,0.1,0.1))
        self.masses = [self.m123, self.m124, self.m135, self.m145, self.m236, self.m246, self.m356, self.m456]

        self.set_springs()

    def set_2(self):
        self.m123 = Mass(0.1, self.cc.m145.p + vector(0.1,0.,0.))
        self.m124 = Mass(0.1, self.cc.m124.p + vector(0.1,0.,0.))
        self.m236 = Mass(0.1, self.cc.m236.p + vector(0.1,0.,0.))
        self.m246 = Mass(0.1, self.cc.m246.p + vector(0.1,0.,0.))
        self.masses = [self.m123, self.m124, self.m236, self.m246]

        #self.set_springs()
        """
        self.m135 = self.cc.m123
        self.m145 = self.cc.m124
        self.m356 = self.cc.m236
        self.m456 = self.cc.m246

        for m1 in self.masses:
            for m2 in [self.m135, self.m145, self.m356, self.m456]:
                p1 = m1.p
                p2 = m2.p
                L0 = sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2 )
                a = L0
                s = Spring(m1,m2,self.masses,ks,L0,a)
                self.springs.append(s)
        """

    def set_4(self):
        self.m124 = Mass(0.1, self.cc.m124.p + vector(0.,0.1,0.))
        self.m145 = Mass(0.1, self.cc.m145.p + vector(0.,0.1,0.))
        self.m246 = Mass(0.1, self.cc.m246.p + vector(0.,0.1,0.))
        self.m456 = Mass(0.1, self.cc.m456.p + vector(0.,0.1,0.))
        self.masses = [self.m124, self.m145, self.m246, self.m456]

        self.set_springs()

        self.m123 = self.cc.m124
        self.m135 = self.cc.m145
        self.m236 = self.cc.m246
        self.m356 = self.cc.m456

        for m1 in self.masses:
            for m2 in [self.m123, self.m135, self.m236, self.m356]:
                p1 = m1.p
                p2 = m2.p
                L0 = sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2 )
                a = L0
                s = Spring(m1,m2,self.masses,ks,L0,a)
                self.springs.append(s)

    def set_masses(self):
        self.masses = []
        for i,v in enumerate(self.cube_vertices):
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
                a = L0
                s = Spring(self.masses[m1],self.masses[m2],ks,L0,a)
                #s = Spring(m1,m2,self.masses,ks,L0,a)
                self.springs.append(s)
              
                m2 = m2+1
            m1 = m1+1

    def initialize_img(self):
        P = vector(0.,0.,0.)
        max_x = -10
        min_x = 10
        max_y = -10
        min_y = 10
        max_z = -10
        min_z = 10

        for n,m in enumerate(self.masses):
            P += m.p
            if m.p.x > max_x:
                max_x = m.p.x
            if m.p.x < min_x:
                min_x = m.p.x
            if m.p.y > max_y:
                max_y = m.p.y
            if m.p.y < min_y:
                min_y = m.p.y
            if m.p.z > max_z:
                max_z = m.p.z
            if m.p.z < min_z:
                min_z = m.p.z

        P = P / n
        L = max_x - min_x
        H = max_y - min_y
        W = max_z - min_z

        self.img = box(pos=P, length=L, height=H, width=W, color=color.orange)

    def update_box(self):
        P = vector(0.,0.,0.)
        max_x = -10
        min_x = 10
        max_y = -10
        min_y = 10
        max_z = -10
        min_z = 10

        for n,m in enumerate(self.masses):
            P += m.p
            if m.p.x > max_x:
                max_x = m.p.x
            if m.p.x < min_x:
                min_x = m.p.x
            if m.p.y > max_y:
                max_y = m.p.y
            if m.p.y < min_y:
                min_y = m.p.y
            if m.p.z > max_z:
                max_z = m.p.z
            if m.p.z < min_z:
                min_z = m.p.z

        P = P / n
        L = max_x - min_x
        H = max_y - min_y
        W = max_z - min_z

        self.img.pos = P
        self.img.length = L
        self.img.height = H
        self.img.width = W

    def update_img(self):
        #self.update_box()
        for mass in self.masses:
            mass.update_img()
        for spring in self.springs:
            spring.update_img(self.masses)


class Robot:

    def __init__(self, type):
        self.cubes = []
        c0 = Cube()
        self.cubes.append(c0)
        if type == 2:
            prev = c0
            for i in range(3):
                c = Cube(prev, 4)
                self.cubes.append(c)
                prev = c
            #c5 = Cube(self.cubes[2], 2)
            #self.cubes.append(c5)
            #c6 = Cube(c5,2)
            #self.cubes.append(c6)

        self.springs = []
        for cube in self.cubes:
            for s in cube.springs:
                self.springs.append(s)

    def update_img(self):
        for cube in self.cubes:
            cube.update_img()

    # move structure
    def translate(self, dx, dy, dz):
        for cube in self.cubes:
            for mass in cube.masses:
                mass.p.x += dx
                mass.p.y += dy
                mass.p.z += dz
      
        self.update_img()

    def get_energy(self):
        potential = 0
        kinetic = 0

        for cube in self.cubes: # FIX WHEN MULTIPLE CUBES
            for mass in cube.masses:
                potential += abs(mass.m * g.z * mass.p.z)
                kinetic += (1/2) * mass.m * (mass.v.x**2 + mass.v.y**2 + mass.v.z**2)

            for spring in cube.springs:
                v1 = cube.masses[spring.m1].p
                v2 = cube.masses[spring.m2].p
                L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )
                potential += (1/2) * spring.k * (spring.L0 - L)**2

        total = potential + kinetic
        return T, potential, kinetic, total

    # update for one time step
    def step(self, t, breathe=False):

        for cube in self.cubes:
            for mass in cube.masses:
                
                # sum all forces on the mass (from springs, gravity, collision with the ground, etc)
                Fs = vector(0,0,0)      # spring force
                for s in self.springs:
                    if breathe==True:
                        s.breathe(t)

                    if s.m1 == mass:
                        v1 = s.m1.p
                        v2 = s.m2.p
                        L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )

                        fmag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
                        #fmag = 1000*(L-s.L0)
                        dir = (v2 - v1)/L
                        Fs = Fs + fmag*dir
                    elif s.m2 == mass:
                        v1 = s.m1.p
                        v2 = s.m2.p
                        L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )
                        #print(type(s))
                        #print(s.m1)
                        #print(s.k)
                        fmag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
                        #fmag = 1000*(L-s.L0)
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
                    mass.v = mass.v * 0.9995
                mass.p = mass.p + mass.v * dt

def main():
    animate = True

    scene.camera.pos = vector(0,-0.6,0.5)
    scene.camera.axis = vector(0,0,0) - scene.camera.pos
    floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)

    cube_robot = Robot(1)
    spring_types = []
    for i in range(len(cube_robot.springs)):
        draw = random()
        if draw < 0.3:     # spring type 1 - tissue
            spring_types.append(1)
        elif draw < 0.6:    # spring type 2 - bone
            spring_types.append(2)
        elif draw < 0.8:   # spring type 3 - muscle 1
            spring_types.append(3)
        else:               # spring type 4 - muscle 2
            spring_types.append(4)

    for s in cube_robot.springs:
        draw = random()
        if draw < 0.25:     # spring type 1 - tissue
            s.k = 1000
            s.b = 0
            s.c = 0
            s.img.color = color.cyan
        elif draw < 0.5:    # spring type 2 - bone
            s.k = 20000
            s.b = 0
            s.c = 0
            s.img.color = color.white
        elif draw < 0.75:   # spring type 3 - muscle 1
            s.k = 5000
            s.b = 0.02
            s.c = 0
            s.img.color = color.red
        else:               # spring type 4 - muscle 2
            s.k = 5000
            s.b = 0.02
            s.c = pi
            s.img.color = color.magenta

    cube_robot.translate(0.,0.,0.2)

    T = 0
    #i = 0
    while T < 5:
        rate(10000) # dt = 0.0001, rate is times per second. 1 sec /0.0001 = 10,000 times. so real-time
        cube_robot.step(T, False)

        if animate == True and T % 0.005:
            cube_robot.update_img()

        #if i % 500 == 0:
        #if T % 0.05 == 0:
            #print(cube_robot.get_energy())
        
        T = T + dt
        #i += 1
    
    print("Fin.")

if __name__=='__main__':
    main()

