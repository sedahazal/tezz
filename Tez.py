from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset
import material
import sketch
import part
import section
import mesh
from odbAccess import *
from numpy import float,cumsum,array,abs,where,cos,arange,pi,zeros,unique,append,loadtxt,ones

class Model:
    def __init__(self,model_name,width,height,time_step, duration,mesh_size,file_name,D,layer,K0,void_ratio,damping,friction_angle,dilation_angle):
        self.f = lambda L,x : [L[i:i+1] for i in x]
        self.width = width
        self.height = height
        self.model_name = model_name
        self.time_step = time_step
        self.duration = duration
        self.mesh_size = mesh_size
        self.file_name = file_name
        self.density = D
        self.layer = layer
        self.layer_heights = append([self.height], self.height -  cumsum(ones(self.layer)*self.height/self.layer))
        self.K0 = float(K0)
        self.void_ratio = float(void_ratio)
        self.damping = float(damping)
        self.friction_angle = friction_angle
        self.dilation_angle = dilation_angle
        self.inf_size = width/8
        mdb.models.changeKey(fromName="Model-1",toName=model_name)
        self.model = mdb.models[model_name]
    
    def create_part(self):
        soilprofile = self.model.ConstrainedSketch(name="__profile__",sheetSize=self.width*2)
        soilprofile.rectangle(point1=(0,0),point2=(self.width,self.height))
    
        self.soilPart = self.model.Part(name="Soil Part",dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
        self.soilPart.BaseShell(sketch=soilprofile)
    
        face_index = self.soilPart.faces.findAt((self.width/2,self.height/2,0)).index
        self.soilPart.Set(faces=self.f(self.soilPart.faces,[face_index,]),name="SoilFace")
    
        del self.model.sketches["__profile__"]
        self.model.ConstrainedSketch(gridSpacing=3.53, name='__profile__', sheetSize=self.width * 2,
                                    transform=self.soilPart.MakeSketchTransform(
                                        sketchPlane=self.soilPart.faces[0],
                                        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0.0)))
    
        self.soilPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch = self.model.sketches['__profile__'])
        self.model.sketches['__profile__'].Line(point1=(self.inf_size , self.height), point2=(self.inf_size, 0))
        self.model.sketches['__profile__'].Line(point1=(self.width - self.inf_size , self.height), point2=(self.width - self.inf_size, 0))
        self.soilPart.PartitionFaceBySketch(faces=self.soilPart.faces[0],
                                      sketch=self.model.sketches['__profile__'])
    
        del self.model.sketches["__profile__"]
    def face_partition(self):
        for y in self.layer_heights[1:-1]:
            face = self.soilPart.faces.findAt((self.width/2,y-0.1,0))
            self.model.ConstrainedSketch(gridSpacing=3.53, name='__profile__', sheetSize=self.width * 2,
                                         transform=self.soilPart.MakeSketchTransform(
                                             sketchPlane=face,
                                             sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0.0)))

            self.soilPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=self.model.sketches['__profile__'])
            self.model.sketches['__profile__'].Line(point1=(0, y), point2=(self.width, y))
            self.soilPart.PartitionFaceBySketch(faces=face,
                                                sketch=self.model.sketches['__profile__'])

            del self.model.sketches["__profile__"]
        for y in self.layer_heights[1:-1]:
            face = self.soilPart.faces.findAt((0.1,y-0.1,0))
            self.model.ConstrainedSketch(gridSpacing=3.53, name='__profile__', sheetSize=self.width * 2,
                                         transform=self.soilPart.MakeSketchTransform(
                                             sketchPlane=face,
                                             sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0.0)))

            self.soilPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=self.model.sketches['__profile__'])
            self.model.sketches['__profile__'].Line(point1=(0, y), point2=(self.inf_size, y))
            self.soilPart.PartitionFaceBySketch(faces=face,
                                                sketch=self.model.sketches['__profile__'])

            del self.model.sketches["__profile__"]
        for y in self.layer_heights[1:-1]:
            face = self.soilPart.faces.findAt((self.width-0.1,y-0.1,0))
            self.model.ConstrainedSketch(gridSpacing=3.53, name='__profile__', sheetSize=self.width * 2,
                                         transform=self.soilPart.MakeSketchTransform(
                                             sketchPlane=face,
                                             sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0.0)))

            self.soilPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=self.model.sketches['__profile__'])
            self.model.sketches['__profile__'].Line(point1=(self.width - self.inf_size, y), point2=(self.width, y))
            self.soilPart.PartitionFaceBySketch(faces=face,
                                                sketch=self.model.sketches['__profile__'])

            del self.model.sketches["__profile__"]

    def natural_frequency(self,mode,Vs,H):
        f = (2*mode-1)*Vs/(4*H)
        return f
    def rayleigh_damping(self,modes,Vs,H,damping_ratio):
        if Vs == 0:
            return 0,0
        else:
            wi = self.natural_frequency(modes[0],Vs,H)
            wj = self.natural_frequency(modes[1],Vs,H)
            alpha = damping_ratio*2*wi*wj/(wi+wj)
            beta = damping_ratio*2/(wi+wj)
            return alpha,beta


    def layer_centers(self):
        self.centers = []
        layer_heights = cumsum(ones(self.layer)*self.height/self.layer)
        for i in range(len(layer_heights)):
            if i == 0:
                center = layer_heights[i] / 2
            else:
                center = layer_heights[i - 1] + (layer_heights[i] - layer_heights[i - 1]) / 2
            self.centers.append(center)
    def create_material(self):
        stress_list=cumsum(array(self.centers)*self.density)
        for i in range(self.layer):
            stress= stress_list[i]
            mean_effective=((stress+(2*self.K0*stress))/3)
            Gmax=(3230*((2.97-self.void_ratio)**2)*(mean_effective)**0.5)/(1+self.void_ratio)
            E=(Gmax*2.5)
            Vs=(Gmax/self.density)**0.5
            H=self.height/self.layer
            name="sand{}".format(i)
            sand = self.model.Material(name)
            sand.Density(table=((self.density,),))
            sand.Elastic(table=(((E,0.25),)))
            alpha, beta = self.rayleigh_damping([1, 4], Vs, H, self.damping)
            sand.Damping(alpha=alpha, beta=beta)
            #sand.MohrCoulombPlasticity(table=((self.friction_angle,self.dilation_angle),))
            #sand.mohrCoulombPlasticity.MohrCoulombHardening(table=((0.1, 0.0),))
            #sand.mohrCoulombPlasticity.TensionCutOff(dependencies=0, table=((0.0, 0.0),), temperatureDependency=OFF)

    

    def create_section(self):
        for i in range(self.layer):
            name = "sand{}".format(i)
            face1=self.soilPart.faces.findAt((self.width/2,self.height-self.centers[i],0)).index
            face2=self.soilPart.faces.findAt((0.1,self.height-self.centers[i],0)).index
            face3=self.soilPart.faces.findAt((self.width-0.1,self.height-self.centers[i],0)).index
            self.soilPart.Set(faces=self.f(self.soilPart.faces,[face1,face2,face3]),name=name)
            self.model.HomogeneousSolidSection(name=name+"Section",material=name)
            self.soilPart.SectionAssignment(region=self.soilPart.sets[name],sectionName=name+"Section")

    def create_instance(self):
        self.model.rootAssembly.DatumCsysByDefault(CARTESIAN)
        self.model.rootAssembly.Instance(dependent=ON,name="Soil Part-1", part=self.soilPart)
        self.model.rootAssembly.regenerate()
    
    def edge_sets(self):
        x1 = (0,self.height/2,0)
        x2 = (self.width,self.height/2,0)
        x3 = (self.inf_size,self.height-0.01,0)
        x4 =(self.width - self.inf_size,self.height-0.01,0)
        y = (self.width/2,0,0)
        y2=(self.width/2,self.height,0)
    
        left_edge = self.soilPart.edges.findAt(x1).index
        left_edge2 = self.soilPart.edges.findAt(x3).index
        left_edge3 = self.soilPart.edges.findAt(x4).index
        right_edge = self.soilPart.edges.findAt(x2).index
        bottom_edge = self.soilPart.edges.findAt(y).index
        top_edge = self.soilPart.edges.findAt(y2).index
    
        self.soilPart.Set(edges=self.f(self.soilPart.edges,[left_edge,right_edge,left_edge2,left_edge3,bottom_edge,top_edge]),name = "Edges")
    
        edge1 = self.soilPart.edges.findAt((0.01,self.height,0)).index
        edge2 = self.soilPart.edges.findAt((0.01,0,0)).index
        edge3 = self.soilPart.edges.findAt((self.width-0.01,self.height,0)).index
        edge4 = self.soilPart.edges.findAt((self.width-0.01,0,0)).index
        self.soilPart.Set(edges=self.f(self.soilPart.edges, [bottom_edge,edge2,edge4 ]), name="HorizontalEdges")
        self.soilPart.Set(edges=self.f(self.soilPart.edges, [edge1,edge2,edge3,edge4]), name="singleseed")
    def boundary_conditions(self):
        sets = self.model.rootAssembly.instances["Soil Part-1"].sets
    
        self.model.DisplacementBC(createStepName = "Initial",name="BC-Y",region=sets["HorizontalEdges"],u1=0,u2=0,u3=0)
        self.model.AccelerationBC(name="Vibration",createStepName="VibrationStep", a1= 1, amplitude= "vibration",
                             region= self.model.rootAssembly.instances["Soil Part-1"].sets["HorizontalEdges"])
    
    def create_step(self):
        self.model.StaticStep(name="GravityStep",previous="Initial")
        self.model.ImplicitDynamicsStep(initialInc=self.time_step, timePeriod= self.duration, maxInc= self.time_step, maxNumInc=1000000,
                                  name="VibrationStep",previous="GravityStep")
    
    
    def gravity_load(self):
        self.model.Gravity(comp2 = -9.81,createStepName="GravityStep",name="Gravity_Load")
    
    
    def geostaticstress(self):
        s2 = -self.density*self.height*9.81
        lc = 0.33
        self.model.GeostaticStress(lateralCoeff1 = lc,name="GS",region=self.model.rootAssembly.instances['Soil Part-1'].sets["SoilFace"],
                              stressMag1 = 0,stressMag2 = s2,vCoord1=self.height,vCoord2 = 0)
    def create_faces(self):
        for i in range(self.layer):
            y = self.centers[i]
            x1 = (self.width / 4-0.01, y , 0)
            x2 = (3 * self.width / 4+0.01, y, 0)
            left_face = self.soilPart.faces.findAt(x1).index
            right_face = self.soilPart.faces.findAt(x2).index
            self.soilPart.Set(faces=self.f(self.soilPart.faces, [left_face, right_face]), name="infinite_face"+str(i))

    def set_mesh_control(self):
        for y in self.layer_heights[:-1]:
            face1 = self.soilPart.faces.findAt((self.width - 0.01,y-0.01,0))
            face2 = self.soilPart.faces.findAt((0.01,y-0.01,0))
            edge1 = self.soilPart.edges.findAt((self.width-0.01,y,0))
            edge2 = self.soilPart.edges.findAt((0.01,y,0))
            self.soilPart.setSweepPath(edge=edge1, region=face1, sense=FORWARD)
            self.soilPart.setSweepPath(edge=edge2, region=face2, sense=REVERSE)

    def generate_mesh(self):
        self.set_mesh_control()
        self.soilPart.setElementType(elemTypes=(mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
            hourglassControl=DEFAULT, distortionControl=DEFAULT), mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)),regions=(self.soilPart.faces,))
        for i in range(self.layer):
            self.soilPart.setElementType(elemTypes=(mesh.ElemType(elemCode=CPE4I, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
                                                                  hourglassControl=DEFAULT, distortionControl=DEFAULT),
                                                    mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)),
                                         regions=(self.soilPart.sets["infinite_face"+str(i)]))
        self.soilPart.seedEdgeBySize(edges=self.soilPart.sets["Edges"].edges, size=self.mesh_size,
                                     constraint=FIXED)
    
        self.soilPart.seedEdgeByNumber(edges=self.soilPart.sets["singleseed"].edges, number=1,
                                       constraint=FIXED)
        self.soilPart.generateMesh()
    
    def create_nodes(self):
        node_list = []
    
        coordinates = {"a1":(32,0,0),
                       "a2":(32,6,0),
                       "a3":(40.8,6,0),
                       "a4":(40.8,12.64,0),
                       "a5":(40.8,13.64,0),
                       "a6":(40.8,14.64,0),
                       "a7":(40.8,16.04,0),
                       "a8":(32,12.64,0),
                       "a9":(32,13.64,0),
                       "a10":(32,14.64,0),
                       "a11":(29,12.64,0),
                       "a12":(29,13.64,0),
                       "a13":(29,14.64,0),
                       "a14":(29,16.04,0),
                        "a15":(32,16.04,0)}

        for node in coordinates:
            node_list.append(self.soilPart.nodes.getClosest(coordinates[node]))
    
        self.soilPart.Set(nodes=mesh.MeshNodeArray(node_list),name="Accelometers")
    
    def create_vibration(self):
        time,accelarations = loadtxt(self.file_name,skiprows=0,unpack=True)
        data = [[time[i], accelarations[i]] for i in range (len(time))]
    
        self.model.TabularAmplitude(name="vibration",timeSpan=STEP, smooth=SOLVER_DEFAULT, data=data)

    def create_history_output(self):
        self.model.fieldOutputRequests["F-Output-1"].deactivate("GravityStep")
        self.model.fieldOutputRequests["F-Output-1"].deactivate("VibrationStep")
        self.model.HistoryOutputRequest(createStepName = "VibrationStep", frequency = 1, name = "H-Output-2", variables = ('A2','A1'),
                                   region= self.model.rootAssembly.allInstances['Soil Part-1'].sets['Accelometers'])
        del self.model.historyOutputRequests ['H-Output-1']
        del self.model.fieldOutputRequests['F-Output-1']
    
    def create_job(self):
        mdb.Job(model=self.model_name,name= self.model_name)
        mdb.jobs[self.model_name].writeInput(consistencyChecking=OFF)
    
    def change_element_type(self):
        inp = open(self.model_name + ".inp")
        data = inp.read().replace("CPE4I","CINPE4")
        new_inp = open(self.model_name + ".inp", "w")
        new_inp.write(data)
        new_inp.close()
        inp.close()

    def operator(self):
        self.create_part()
        self.layer_centers()
        self.create_material()
        self.edge_sets()
        self.face_partition()
        self.create_section()
        self.create_instance()
        self.create_step()
        self.create_faces()
        #self.geostaticstress()
        #self.gravity_load()
        self.create_vibration()
        self.boundary_conditions()
        self.generate_mesh()
        self.create_nodes()
        self.create_history_output()
        self.create_job()
        self.change_element_type()

model=Model("deneme",64,16.64,0.00313,92.84519,0.5,"C:\\Users\\gpumachine\\Desktop\\hazal tez\\self025g-70.txt",1.581,4,0.388,0.64,0.005,38,12)
model.operator()