# Zena-Marie Gonzalez (zg2326)
# MECS E4510: Evolutionary Computation and Design Algorithms
# Soft Robot Evolution - Morphology and Control System

from os import supports_dir_fd
from random import randint
from numpy import true_divide
from numpy.lib.index_tricks import fill_diagonal
from vpython import *

""" Simulation Parameters """
g = vector(0,0,-9.81)   # gravity, meters/s^2
dt = 0.0001             # timestep, seconds
#T = 0                   # time, seconds
kc = 100000             # ground restoration constant
ks = 1000               # default spring constant
dampening = True        # add dampening
w = 30                  # global frequency
mu_s = 1                # friction coefficient
mu_k = 0.8              # kinetic coefficient

""" Mutation Parameters """
add_prob = 0.4          # base probability of adding a Type_Center
del_prob = 0.4          # base probability of deleting a Type_Center
swap_prob = 0.8         # base probability of swapping the positions of 2 Type Centers in the encoding
                        # ^ order determines assignment - later in the encoding overrides earlier
type_change_prob = 0.8  # base probability of changing the type of one Type Center
rad_change_prob = 0.6   # base probability of changing the radius of a Type Center
pos_change_prob = 0.9   # base probability of changing the position of a Type Center

""" Evolution Parameters """
population_size = 20    # population size
recombine_prob = 0.25   # probability of recombining instead of mutating
lc_datapoints = 100     # number of datapoints in the learning curve

class Mass:

    def __init__(self, draw=True, m=1, p=vector(0,0,0), v=vector(0,0,0), a=vector(0,0,0)) -> None:
        self.m = m      # mass, kg
        self.p = p      # position, meters
        self.v = v      # velocity, meters/sec
        self.a = a      # acceleration, meters/sec^2
        
        if draw == True:
            self.img = sphere(pos=vector(self.p), radius=0.01, color=color.orange)
        
    def update_img(self):
        self.img.pos = self.p

class Spring:

    def __init__(self, m1, m2, draw=True, k=10000, L0=1, a=1, b=.005, c=0, col=color.blue) -> None:
        self.m1 = m1    # m1,m2 are the two masses self connects
        self.m2 = m2
        self.k = k      # spring constant, N/meter
        self.L0 = L0    # original resting length, meters
        self.a = a      # spring parameters
        self.b = b
        self.c = c
        
        if draw == True: 
            ax = self.m2.p - self.m1.p
            self.img = cylinder(pos=self.m1.p, axis=ax, radius=0.005, color=col)
    
    def update_img(self,masses):
        self.img.pos = self.m1.p
        self.img.axis = self.m2.p - self.m1.p

    def breathe(self, t):
        self.L0 = self.a + self.b*sin(w*t + self.c)

class Cube:

    def __init__(self, draw=True, connecting_cube=None, connecting_side=None) -> None:
     
        self.cc = connecting_cube
        self.draw = draw

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
        self.m123 = Mass(self.draw, 0.1, vector(0.,0.,0.))
        self.m124 = Mass(self.draw, 0.1, vector(0.,0.1,0.))
        self.m135 = Mass(self.draw, 0.1, vector(0.1,0.,0.))
        self.m145 = Mass(self.draw, 0.1, vector(0.1,0.1,0.))
        self.m236 = Mass(self.draw, 0.1, vector(0.,0.,0.1))
        self.m246 = Mass(self.draw, 0.1, vector(0.,0.1,0.1))
        self.m356 = Mass(self.draw, 0.1, vector(0.1,0.,0.1))
        self.m456 = Mass(self.draw, 0.1, vector(0.1,0.1,0.1))
        self.masses = [self.m123, self.m124, self.m135, self.m145, self.m236, self.m246, self.m356, self.m456]

        self.set_springs()

    def set_2(self):
        self.m123 = Mass(self.draw, 0.1, self.cc.m145.p + vector(0.1,0.,0.))
        self.m124 = Mass(self.draw, 0.1, self.cc.m124.p + vector(0.1,0.,0.))
        self.m236 = Mass(self.draw, 0.1, self.cc.m236.p + vector(0.1,0.,0.))
        self.m246 = Mass(self.draw, 0.1, self.cc.m246.p + vector(0.1,0.,0.))
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
                s = Spring(m1,m2,self.draw,self.masses,ks,L0,a)
                self.springs.append(s)
        """

    def set_4(self):
        self.m124 = Mass(self.draw, 0.1, self.cc.m124.p + vector(0.,0.1,0.))
        self.m145 = Mass(self.draw, 0.1, self.cc.m145.p + vector(0.,0.1,0.))
        self.m246 = Mass(self.draw, 0.1, self.cc.m246.p + vector(0.,0.1,0.))
        self.m456 = Mass(self.draw, 0.1, self.cc.m456.p + vector(0.,0.1,0.))
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
                #s = Spring(m1,m2,self.draw,20000,L0,a,0,0,color.white)
                s = Spring(m1,m2,self.draw,1000,L0,a,0,0,color.cyan)
                self.springs.append(s)

    def set_masses(self):
        self.masses = []
        for i,v in enumerate(self.cube_vertices):
            m = Mass(self.draw, 0.1, v)
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
                s = Spring(self.masses[m1],self.masses[m2],self.draw,1000,L0,a,0,0,color.cyan)
                #s = Spring(self.masses[m1],self.masses[m2],self.draw,20000,L0,a,0,0,color.white)
                #s = Spring(m1,m2,self.masses,ks,L0,a)
                self.springs.append(s)
              
                m2 = m2+1
            m1 = m1+1

class Robot:

    def __init__(self, type, draw=True):
        self.draw = draw
        self.cubes = []
        c0 = Cube(self.draw)
        self.cubes.append(c0)
        if type == 2:
            prev = c0
            for i in range(3):
                c = Cube(self.draw, prev, 4)
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
        self.masses = []
        for cube in self.cubes:
            for m in cube.masses:
                self.masses.append(m)
   
    def apply_encoding(self, encoding):

        ghost_springs = []

        for tc in encoding:
            for spring in self.springs:
                if self.distance(spring.m1.p, tc.p) < tc.r or self.distance(spring.m2.p, tc.p) < tc.r:

                    if spring in ghost_springs:
                        ghost_springs.remove(spring)

                    if tc.type == 1:        # spring type 1 - tissue
                        spring.k = 1000
                        spring.b = 0
                        spring.c = 0
                        if self.draw == True:
                            spring.img.color = color.cyan
                    elif tc.type == 2:      # spring type 2 - bone
                        spring.k = 20000
                        spring.b = 0
                        spring.c = 0
                        if self.draw == True:
                            spring.img.color = color.white
                    elif tc.type == 3:      # spring type 3 - muscle 1
                        spring.k = 5000
                        spring.b = 0.02
                        spring.c = 0
                        if self.draw == True:
                            spring.img.color = color.red
                    elif tc.type == 4:      # spring type 4 - muscle 2
                        spring.k = 5000
                        spring.b = 0.02
                        spring.c = pi
                        if self.draw == True:
                            spring.img.color = color.magenta
                    else:                   # spring type 5 - air
                        ghost_springs.append(spring)

        floating_masses = self.masses.copy()
        #ghost_springs = self.springs.copy()

        for spring in self.springs:
            if not spring in ghost_springs:
                if spring.m1 in floating_masses:
                    floating_masses.remove(spring.m1)
                if spring.m2 in floating_masses:
                    floating_masses.remove(spring.m2)

        for s in ghost_springs:
            if self.draw == True:
                s.img.visible = False
                del s.img
            self.springs.remove(s)
        for m in floating_masses:
            if self.draw == True:
                m.img.visible = False
                del m.img
            self.masses.remove(m)

    def distance(self, p1, p2):
        return mag(p1-p2)
     
    def update_img(self):
        for mass in self.masses:
            mass.update_img()
        for spring in self.springs:
            spring.update_img(self.masses)

    # move structure
    def translate(self, dx, dy, dz):
        for cube in self.cubes:
            for mass in cube.masses:
                mass.p.x += dx
                mass.p.y += dy
                mass.p.z += dz
      
        if self.draw == True:
            self.update_img()

    def get_energy(self, time):
        potential = 0
        kinetic = 0

        for mass in self.masses:
            potential += abs(mass.m * g.z * mass.p.z)
            kinetic += (1/2) * mass.m * (mass.v.x**2 + mass.v.y**2 + mass.v.z**2)

        for spring in self.springs:
            v1 = spring.m1.p
            v2 = spring.m2.p
            L = sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2 )
            potential += (1/2) * spring.k * (spring.L0 - L)**2

        total = potential + kinetic
        return time, potential, kinetic, total

    # get center of mass
    def get_cm(self):
        sum_mass_positions = vector(0.,0.,0.)
        for mass in self.masses:
            sum_mass_positions = sum_mass_positions + mass.p
        cm = sum_mass_positions / len(self.masses)
        return(cm)

    # update for one time step
    def step(self, t, breathe=False):

        for mass in self.masses:
                
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
                    fmag = s.k*(L-s.L0) # compression(L>L0) or tension(L<L0) in direction of spring
                    dir = (v1 - v2)/L
                    Fs = Fs + fmag*dir

            Fg = mass.m * g     # gravitational force
            Fe = vector(0,0,0)  # external force
            Fc = vector(0,0,0)  # restoration force from ground

            F = Fs + Fg
            if mass.p.z < 0:
                Fh = sqrt(F.x**2 + F.y**2)    # horizontal force
                Fn = F.z                      # normal force
                if Fn < 0:
                    if Fh < (-Fn)*mu_s:
                        Fe.z = -Fh
                    else:
                        Fe = -1*F / Fh  # normalize (horizontal components of) F and in opposite direction
                        Fe.z = 0        # only want horizontal forces
                        Fe = Fe * Fn * mu_k # magnitude is Fn*mu_k
                Fc.z = -1*kc*mass.p.z

            F = F + Fe + Fc # total force on mass
                
            mass.a = F/mass.m
            mass.v = mass.v + mass.a * dt
            if dampening == True:
                mass.v = mass.v * 0.990 #0.9995
            mass.p = mass.p + mass.v * dt

class Type_Center:

    def __init__(self, type, p=vector(0,0,0), r=0.05) -> None:
        self.type = type    # 1=tissue, 2=bone, 3=muscle1, 4=muscle2, 5=air
        self.p = p          # position of center
        self.r = r          # radius

        #self.img = sphere(pos=vector(self.p), radius=0.01, color=color.blue)

    def to_string(self):
        s = "enc.append(Type_Center(" + str(self.type) + ", vector(" + str(self.p.x) + "," + str(self.p.y) + ',' + str(self.p.z) + "), " + str(self.r) + "))"
        return(s)

def random_encoding(length):
    encoding = []

    for i in range(length):
        x = random() * 0.18 - 0.04
        y = random() * 0.48 - 0.04
        z = random() * 0.18 - 0.04
        type = randint(1,5)
        
        tc = Type_Center(type, vector(x,y,z))
        encoding.append(tc)

    return encoding

def evaluate(encoding):

    robot = Robot(2, False)
    robot.apply_encoding(encoding)

    cm_before = robot.get_cm()
    cm_before.z = 0

    T = 0
    i = 0
    while i < 500:
        robot.step(T, True) 
        T = T + dt
        i += 1

    cm_after = robot.get_cm()
    cm_after.z = 0

    dist_traveled = mag(cm_after - cm_before)
    speed = dist_traveled / T
    return speed

def mutate(encoding):
    if random() < add_prob:
        tc = random_encoding(1)
        encoding.append(tc[0])
        #print("tc added")
    if len(encoding) > 1 and random() < del_prob:
        pos = randint(0,len(encoding)-1)
        encoding.remove(encoding[pos])
        #print("tc at ", pos, " deleted")
    if random() < swap_prob:
        pos1 = randint(0,len(encoding)-1)
        pos2 = randint(0,len(encoding)-1)
        temp = encoding[pos1]
        encoding[pos1] = encoding[pos2]
        encoding[pos2] = temp
        #print("tcs at ", pos1, " and ", pos2, " swaped")
    if random() < type_change_prob:
        new_type = randint(1,5)
        pos = randint(0,len(encoding)-1)
        encoding[pos].type = new_type
        #print("type at ", pos, " changed to ", new_type)
    if random() < rad_change_prob:
        factor = random() + 0.5
        pos = randint(0,len(encoding)-1)
        encoding[pos].r = encoding[pos].r * factor
        #print("radius at ", pos, " changed by ", factor)
    if random() < pos_change_prob:
        a = 0.1 # ASSIGN
        b = a/2
        delta = vector(random()*a-b, random()*a-b, random()*a-b)
        pos = randint(0,len(encoding)-1)
        new_p = encoding[pos].p + delta
        if new_p.x > -0.04 and new_p.x < 0.14 and new_p.y > -0.04 and new_p.y < 0.44 and new_p.z > -0.04 and new_p.z < 0.14:
            encoding[pos].p = new_p
            #print("position at ", pos, "changed to ", new_p, " by ", delta)
        else:
            new_p = encoding[pos].p - delta
            if new_p.x > -0.04 and new_p.x < 0.22 and new_p.y > -0.04 and new_p.y < 0.42 and new_p.z > -0.04 and new_p.z < 0.22:
                encoding[pos].p = new_p
                #print("position at ", pos, "changed to ", new_p, " by ", delta)

    return encoding

def recombine(par1, par2):
    m1 = random()
    m2 = random()
    if m1 > m2:
        temp = m1
        m1 = m2
        m2 = temp

    s1 = round(m1*len(par1))
    e1 = round(m2*len(par1))
    s2 = round(m1*len(par2))
    e2 = round(m2*len(par2))

    """
    s1 = randint(0, len(par1)-1)
    e1 = randint(s1, len(par1)-1)
    s2 = randint(0, len(par2)-1)
    e2 = randint(s2, len(par2)-1)
    """

    child = par1[0:s1] + par2[s2:e2] + par1[e1:len(par1)]
    return child

def find_ind(val, list):
    for i in range(len(list)):
        if list[i] == val:
            return i
    print("ERROR: index not found")
    return 0

def random_search(evals):
    best_encoding = random_encoding(5)
    best_speed = evaluate(best_encoding)

    benchmark = round(evals / lc_datapoints)

    for e in range(evals):
        new_e = random_encoding(5)
        new_e_speed = evaluate(new_e)

        if new_e_speed > best_speed:
            best_encoding = new_e
            best_speed = new_e_speed
        
        if e % benchmark == 0:
            print(e, best_speed)

    return best_encoding, best_speed

def parallel_hillclimber(evals, n_solutions):

    #benchmark = round(evals / lc_datapoints)
    #next_bench = benchmark
    #eval_count = 0

    population = []
    scores = []
    for j in range(population_size):
        enc = random_encoding(5)
        population.append(enc)
        scores.append(evaluate(enc))
    #print("Base pop created")
    print(0, scores[0])
    eval_count = population_size
    print(eval_count, scores[0])

    generations = round(evals / population_size) - 1 # - 1 for evaluating initial population
    #cohort_size = population_size / 4
    for g in range(generations):
        child_pop = []
        child_scores = []
        for i in range(population_size):
            child = mutate(population[i])
            child_pop.append(child)
            child_scores.append(evaluate(child))

        new_pop = []
        new_scores = []

        max_par = max(scores)
        max_chl = max(child_scores)
        while len(new_pop) < population_size:
            
            if len(scores) > 1 and max_par > max_chl:
                i = find_ind(max_par, scores)
                new_pop.append(population[i])
                new_scores.append(scores[i])
                population.remove(population[i])
                scores.remove(scores[i])
                max_par = max(scores)

            else:   # if max_chl >= max_par: (priority given to younger solutions)
                i = find_ind(max_chl, child_scores)
                new_pop.append(child_pop[i])
                new_scores.append(child_scores[i])
                child_pop.remove(child_pop[i])
                child_scores.remove(child_scores[i])
                max_chl = max(child_scores)

        population = new_pop
        scores = new_scores

        #print("Generation ", g)
        eval_count += population_size
        #if eval_count > next_bench:
        print(eval_count, scores[0])
        #next_bench += benchmark

    best_encodings = population[0:n_solutions+1]
    best_speeds = scores[0:n_solutions+1] 

    return best_encodings, best_speeds

def evolutionary_search(evals, n_solutions):

    benchmark = round(evals / lc_datapoints)

    population = []
    scores = []
    for j in range(population_size):
        enc = random_encoding(5)
        population.append(enc)
        scores.append(evaluate(enc))
    print("Base pop created")
    evals = evals - population_size

    for e in range(evals):
        if random() < recombine_prob:
            ip1 = randint(0, population_size-1)
            ip2 = randint(0, population_size-1)
            child = recombine(population[ip1], population[ip2])
            child_score = evaluate(child)
            ip = ip1
            if scores[ip1] > scores[ip2]:
                ip = ip2
            if child_score > scores[ip]:
                population[ip] = child
                scores[ip] = child_score
        else:
            ip = randint(0, population_size-1) 
            child = mutate(population[ip])
            child_score = evaluate(child)
            if child_score > scores[ip]:
                population[ip] = child
                scores[ip] = child_score

        if e % benchmark == 0:
            best_speed = max(scores)
            print(e, best_speed)

    best_encodings = []
    best_speeds = [] 

    for k in range(n_solutions):
        m = max(scores)
        i = find_ind(m, scores)
        best_encodings.append(population[i])
        best_speeds.append(scores[i])
        population.remove(population[i])
        scores.remove(scores[i])

    return best_encodings, best_speeds

def main_1():

    print("Random Search")

    result = random_search(1000)
    print("Speed: ", result[1])
    encoding = result[0]
    #encoding = random_encoding(5)
    for tc in encoding:
        print(tc.to_string())

    scene.camera.pos = vector(0,-0.6,0.5)
    scene.camera.axis = vector(0,0,0) - scene.camera.pos
    floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)

    robot = Robot(2, True)
    robot.apply_encoding(encoding)

    #"""
    T = 0
    #i = 0
    while T < 5:
        rate(10000) # dt = 0.0001, rate is times per second. 1 sec /0.0001 = 10,000 times. so real-time
        robot.step(T, True)

        if T % 0.005:
            robot.update_img()

        #if i % 500 == 0:
        ## if T % 0.05 == 0:
            #print(cube_robot.get_energy(T))
        
        T = T + dt
        #i += 1
    #"""
    
    print("Fin.")

def main_2():

    result = parallel_hillclimber(1000, 4)
    #result = evolutionary_search(1000, 4)
    best_bots = result[0]
    bot_scores = result[1]

    for i in range(4):
        bot = best_bots[i]
        print("Robot ", i, ": ")
        print("Speed: ", bot_scores[i])
        for tc in bot:
            print(tc.to_string())
        print(" ")

    scene.camera.pos = vector(0,-0.6,0.5)
    scene.camera.axis = vector(0,0,0) - scene.camera.pos
    floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)

    robot = Robot(2, True)
    robot.apply_encoding(best_bots[0])

    #"""
    T = 0
    #i = 0
    while T < 5:
        rate(10000) # dt = 0.0001, rate is times per second. 1 sec /0.0001 = 10,000 times. so real-time
        robot.step(T, True)

        if T % 0.005:
            robot.update_img()

        #if i % 500 == 0:
        ## if T % 0.05 == 0:
            #print(cube_robot.get_energy(T))
        
        T = T + dt
        #i += 1
    #"""
    
    print("Fin.")

def main(): # display robot from encoding:
    enc = []
    enc.append(Type_Center(4, vector(0.1307945919550911,0.0395324609026174,0.05512654861076575), 0.029967232503743196))
    enc.append(Type_Center(1, vector(0.12730562256838113,0.18742721729184636,0.08711949318011505), 0.05))
    enc.append(Type_Center(1, vector(-0.005739830150065205,0.07421911656680072,-0.03683555672002461), 0.0610894758873259))
    enc.append(Type_Center(1, vector(0.03967784662869471,0.19499912815014622,0.09732119062865908), 0.15177656496826464))
    enc.append(Type_Center(4, vector(-0.008584886286624474,0.31513674878655873,-0.01724328330314578), 0.033430876113585896))
    enc.append(Type_Center(5, vector(0.02379319606242218,0.009566879616287112,0.014770842262386577), 0.06586270387058234))
    enc.append(Type_Center(1, vector(0.00735581038605506,-0.03808115203946532,0.03913080797404194), 0.05))

    scene.camera.pos = vector(0,-0.6,0.5)
    scene.camera.axis = vector(0,0,0) - scene.camera.pos
    floor = box(pos=vector(0,0,-0.0001), length=2, height=2, width=0.0002)

    robot = Robot(2, True)
    robot.apply_encoding(enc)

    #"""
    T = 0
    #i = 0
    while T < 5:
        rate(10000) # dt = 0.0001, rate is times per second. 1 sec /0.0001 = 10,000 times. so real-time
        robot.step(T, True)

        if T % 0.005:
            robot.update_img()

        #if i % 500 == 0:
        ## if T % 0.05 == 0:
            #print(cube_robot.get_energy(T))
        
        T = T + dt
        #i += 1
    #"""
    
    print("Fin.")

def main_test():
    enc = random_encoding(5)
    for tc in enc:
        print(tc.to_string())

if __name__=='__main__':
    main()
