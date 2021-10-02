"""Created on 2 Oct 2021

Interactive Collision Detection (ICoD) Demo.
@author: Alex Grimwood

run from the command line:
python icod_demo.py

See requirements.txt for dependencies

Ask the author for the necessary .stl files

Based upon previous work by kpiqu
"""
import vtkmodules.all as vtk
from vedo import *
import argparse
import os
import numpy as np

def get_program_parameters():
    ''' Set input parameters from command line '''
    description = 'Collision detection demo.'
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gantry', default='Gantry.stl', type=str, help='gantry stl file name')
    parser.add_argument('--collimator', default='Collimator.stl', type=str, help='collimator stl file name')
    parser.add_argument('--couch', default='TableTop.stl', type=str, help='couch-top stl file name')
    parser.add_argument('--body', default='Body.stl', type=str, help='patient stl file name')
    parser.add_argument('--machine_dir', default='.\stl\VarianTrueBeamSTx', type=str, help='machine stl directory')
    parser.add_argument('--patient_dir', default='.\stl\Patient', type=str, help='patient stl directory')
    args = parser.parse_args()
    return args


def collisionfilter(mesh1, mesh2, clean1=True, clean2=True):
    ''' Create a vtk collision detection filter between two mesh objects '''
    # extract polydata from mesh objects
    if clean1 is True:
        poly1 = mesh1.clean().polydata()
    else:
        poly1 = mesh1.polydata()

    if clean2 is True:
        poly2 = mesh2.clean().polydata()
    else:
        poly2 = mesh2.polydata()

    # object transformation inputs
    tform1 = vtk.vtkTransform()
    tform2 = vtk.vtkTransform()

    # collision detector
    collide = vtk.vtkCollisionDetectionFilter()
    collide.SetInputData(0, poly1)
    collide.SetTransform(0, tform1)
    collide.SetInputData(1, poly2)
    collide.SetTransform(1, tform2)
    collide.SetBoxTolerance(0.0)
    collide.SetCellTolerance(0.0)
    collide.SetNumberOfCellsPerNode(2)
    collide.SetCollisionModeToFirstContact()
    collide.GenerateScalarsOn()
    
    # collision mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(collide.GetContactsOutputPort())
    mapper.SetResolveCoincidentTopologyToPolygonOffset()
    
    # collision actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor((0,0,0))
    actor.GetProperty().SetLineWidth(3.0)
    
    return collide, actor, tform1, tform2


def main():
    # Load CLI args
    args = get_program_parameters()

    # Initialise display window
    plt = Plotter(title="Collision Demo")
        
    # Define colors
    colors = vtk.vtkNamedColors()

    # Specify stl directories
    machine_dir = os.path.normpath(args.machine_dir) 
    patient_dir = os.path.normpath(args.patient_dir)

    # Specify machine stl files
    fileGantry = os.path.join(machine_dir,args.gantry)
    fileCollimator = os.path.join(machine_dir,args.collimator)
    fileCouch = os.path.join(machine_dir,args.couch)
    
    # Specify patient stl file
    filePatient = os.path.join(patient_dir,args.body)
    
    # Specify model colours
    gantry_colour = (0.678431, 0.847059, 0.901961)
    collimator_colour = (0, 0, 0.545098)
    couchtop_colour = (0, 0, 0.545098)
    body_colour = (0.564706, 0.933333, 0.564706)
    
    # Load machine stl
    gantry = Mesh(fileGantry).c(gantry_colour)
    collimator = Mesh(fileCollimator).c(collimator_colour)
    couchtop = Mesh(fileCouch).c(couchtop_colour)
    gantry.name = "gantry"
    collimator.name = "collimator"
    couchtop.name = "couchtop"

    # Load patient stl
    body = Mesh(filePatient).c(body_colour)
    body.name = "body"

    # Adjust machine to realistic orientation
    initial_machine_rot = 90
    initial_couch_offset = 0
    gantry.rotateX(initial_machine_rot)
    collimator.rotateX(initial_machine_rot)
    couchtop.rotateX(initial_machine_rot)
    couchtop.z(initial_couch_offset) 

    # Create collision filters
    collide0, actor0, tform0_1, tform0_2 = collisionfilter(body, collimator)
    collide1, actor1, tform1_1, tform1_2 = collisionfilter(couchtop, collimator)

    # Model colour selection function
    def set_colours(model_list, colour_list):
        """ Change mesh colours """
        for model, colour in zip(model_list, colour_list):
            model_properties = model.GetProperty()
            model_properties.SetColor(colour)

    # Collision details reporting function
    keys = ["Gantry", "Snout", "BodyX", "BodyY", "BodyZ"]
    global geometry
    geometry = {key: 0 for key in keys}

    def collision_vis(collisionfilter, meshes):
        ''' Identify collision objects and their geometries '''
        if collisionfilter.GetNumberOfContacts() > 0:
            print("Collision between "+meshes[0].name+" and "+meshes[1].name+": "+str(geometry))
            set_colours(meshes,[(1,0,0),(1,0,0)])

    # Initialise gantry and snout values
    global g_theta, s_ext
    g_theta = 0
    s_ext = 0

    # Slider callback functions
    def slider_x(widget, event):
        """ Moves patient LEFT-RIGHT """
        value = widget.GetRepresentation().GetValue()
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # move patient
        body.x(value)  # set patient x position
        tform0_1.SetMatrix(body.GetMatrix())
        # record new geometry
        global geometry
        geometry["BodyX"]=round(value,1)
        # detect collisions
        collision_vis(collide0,[body,collimator])
        collision_vis(collide1,[couchtop,collimator])

    def slider_y(widget, event):
        """ Moves patient ANT-POST """
        value = widget.GetRepresentation().GetValue()
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # move patient
        body.y(value)  # set patient y position
        tform0_1.SetMatrix(body.GetMatrix())
        # record new geometry
        global geometry
        geometry["BodyY"]=round(value,1)
        # detect collisions
        collision_vis(collide0,[body,collimator])
        collision_vis(collide1,[couchtop,collimator])

    def slider_z(widget, event):
        """ Moves patient SUP-INF """
        value = widget.GetRepresentation().GetValue()
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # move patient
        body.z(value)  # set patient z position
        tform0_1.SetMatrix(body.GetMatrix())
        # record new geometry
        global geometry
        geometry["BodyZ"]=round(value,1)
        # detect collisions
        collision_vis(collide0,[body,collimator])
        collision_vis(collide1,[couchtop,collimator])

    def slider_g(widget, event):
        """ Rotates gantry and collimator  """
        value = widget.GetRepresentation().GetValue()
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # rotate gantry
        global g_theta
        g_theta = value
        T = vtk.vtkTransform()
        T.RotateZ(value)
        gantry.SetUserMatrix(T.GetMatrix())
        # transform colimator
        global s_ext
        T = vtk.vtkTransform()
        T.RotateZ(g_theta)
        T.Translate(0., s_ext, 0.)
        collimator.SetUserMatrix(T.GetMatrix())
        tform0_2.SetMatrix(collimator.GetUserMatrix())
        tform1_2.SetMatrix(collimator.GetUserMatrix())
        # record new geometry
        global geometry
        geometry["Gantry"]=round(value,1)
        # detect collisions
        collision_vis(collide0,[body,collimator])
        collision_vis(collide1,[couchtop,collimator])

    def slider_s(widget, event):
        """ Extends collimator """
        value = widget.GetRepresentation().GetValue()
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # transform colimator
        global g_theta, s_ext
        s_ext = value
        T = vtk.vtkTransform()
        T.RotateZ(g_theta)
        T.Translate(0., s_ext, 0.)
        collimator.SetUserMatrix(T.GetMatrix())
        tform0_2.SetMatrix(collimator.GetUserMatrix())
        tform1_2.SetMatrix(collimator.GetUserMatrix())
        # record new geometry
        global geometry
        geometry["Snout"]=round(value,1)
        # detect collisions
        collision_vis(collide0,[body,collimator])
        collision_vis(collide1,[couchtop,collimator])
                
    plt.addSlider3D(
        slider_x,
        pos1=[-450., 1250., -2050.],
        pos2=[-450.+1000, 1250., -2050.],
        xmin=-500,
        xmax=500,
        value=0,
        s=0.08,
        c="r",
        rotation=180,
        title="x patient (mm)",
    )
                
    plt.addSlider3D(
        slider_y,
        pos1=[600., -800., -2050.],
        pos2=[600., -800.+2000, -2050.],
        xmin=-500,
        xmax=500,
        value=0,
        s=0.08,
        c="g",
        rotation=90,
        title="y patient (mm)",
    )
    
    plt.addSlider3D(
        slider_z,
        pos1=[600., 1250., -2000.],
        pos2=[600., 1250., -2000+2000.],
        xmin=-500,
        xmax=500,
        value=0,
        s=0.08,
        c="b",
        rotation=180,
        title="z patient (mm)",
    )

    plt.addSlider2D(
        slider_g,
        xmin=-180,
        xmax=180,
        value=0,
        pos=4,
        title='gantry angle (deg)',
        titleSize=1,
        c='m',
        showValue=True,
    )

    plt.addSlider2D(
        slider_s,
        xmin=0,
        xmax=250,
        value=0,
        pos=3,
        title='snout extension (mm)',
        titleSize=1,
        c='c',
        showValue=True,
    )

    # Visualise
    plt.show([body, collimator, gantry, couchtop, actor0, actor1], axes=1, camera={"pos": (0.,0.,-12000.), "viewup": [0,-1,0], "focalPoint": (0,0,0)}, bg='black', bg2='dark grey').close()

if __name__ == '__main__':
    main()