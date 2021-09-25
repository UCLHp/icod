"""Created on Sep 23 2021

@author: Alex Grimwood

Interactive Collision Detection (ICoD) Demo.

run from the command line:
python icod_demo.py

See requirements.txt for dependencies

Ask the author for the necessary .stl files

Based upon previous work by kpiqu
"""
import numpy as np
import vtk
from vedo import Plotter, Mesh
import os

def main():
    # initialise display
    plt = Plotter(title="Collision Demo")

    # specify stl directories
    machine_dir = "stl\\VarianTrueBeamSTx"
    patient_dir = "stl\\Patient"

    # specify machine stl files
    fileGantry = os.path.join(machine_dir,"Gantry.stl")
    fileCollimator = os.path.join(machine_dir,"Collimator.stl")
    fileCouch = os.path.join(machine_dir,"TableTop.stl")
    fileTable = os.path.join(machine_dir,"PatientSupport.stl")
    # specify patient stl file
    filePatient = os.path.join(patient_dir,"BODY.stl")
    
    # specify model colours
    gantry_colour = (0.678431, 0.847059, 0.901961)
    collimator_colour = (0, 0, 0.545098)
    couchtop_colour = (0.678431, 0.847059, 0.901961)
    couch_colour = (0, 0, 0.545098)
    body_colour = (0.564706, 0.933333, 0.564706)
    
    # load machine stl
    gantry = Mesh(fileGantry).c('lb')
    collimator = Mesh(fileCollimator).c('dark blue')
    couchtop = Mesh(fileCouch).c('lb')
    couch = Mesh(fileTable).c('dark blue')
    # load patient stl
    body = Mesh(filePatient).c('light green')

    # adjust geometries to realistic relative positions
    initial_machine_rot = 90
    initial_body_rot = 0
    initial_couch_offset = 0
    gantry.rotateX(initial_machine_rot)
    collimator.rotateX(initial_machine_rot)
    couchtop.rotateX(initial_machine_rot)
    couch.rotateX(initial_machine_rot)
    body.rotateX(initial_body_rot)
    couchtop.z(initial_couch_offset)
    couch.z(initial_couch_offset)   

    # model colour selection functions
    def set_colours(model_list, colour_list):
        """ Change mesh colours """
        for model, colour in zip(model_list, colour_list):
            model_properties = model.GetProperty()
            model_properties.SetColor(colour)
    
    # collision detection functions
    global collision_count
    collision_count = 0 # collision counter
    
    keys = ["GantryAngle", "SnoutExt", "BodyY", "BodyZ"]
    global geometry
    geometry = {key: 0 for key in keys}

    def collivision(a,b, a_name, b_name):
        """ Collision detection between two objects (a and b).
        
        Args:
        a (vedo mesh object)
        b (vedo mesh object)
        
        Returns:
        Collision geometry printed to terminal
        Mesh collisions visualised in red
        """
        ap = a.clean().polydata()
        bp = b.clean().polydata()
        transform0 = vtk.vtkTransform()
        transform1 = vtk.vtkTransform()
        collide = vtk.vtkCollisionDetectionFilter()
        collide.SetInputData(0, ap)
        collide.SetTransform(0, transform0)
        collide.SetInputData(1, bp)
        collide.SetTransform(1, transform1)
        collide.SetBoxTolerance(0.0)
        collide.SetCellTolerance(0.0)
        collide.SetNumberOfCellsPerNode(2)
        collide.SetCollisionModeToFirstContact()
        collide.GenerateScalarsOn()
        collide.Update()
        if collide.GetNumberOfContacts() > 0:
            global geometry
            print("Collision between "+a_name+" and "+b_name+"! "+str(geometry))
            set_colours([a,b],[(1,0,0),(1,0,0)])

    # slider callback functions
    def slider_y(widget, event):
        """ Moves patient ANT-POST """
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # move patient
        value = widget.GetRepresentation().GetValue()
        body.y(value)  # set patient y position
        # record new geometry
        global geometry
        geometry["BodyY"]=round(value,1)
        # detect collisions
        collivision(collimator,couchtop, "collimator", "couchtop")
        collivision(collimator,body, "collimator", "body")

    def slider_z(widget, event):
        """ Moves patient SUP-INF """
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        value = widget.GetRepresentation().GetValue()
        body.z(value)  # set patient z position
        # record new geometry
        global geometry
        geometry["BodyZ"]=round(value,1)
        # detect collisions
        collivision(collimator,couchtop, "collimator", "couchtop")
        collivision(collimator,body, "collimator", "body")

    def slider_g(widget, event):
        """ Rotates gantry and collimator  """
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # rotate gantry
        value = widget.GetRepresentation().GetValue()
        T = vtk.vtkTransform()
        T.RotateZ(value)
        gantry.SetUserTransform(T)
        collimator.SetUserTransform(T)  # set gantry angle
        # record new geometry
        global geometry
        geometry["GantryAngle"]=round(value,1)
        # detect collisions
        collivision(collimator,couchtop, "collimator", "couchtop")
        collivision(collimator,body, "collimator", "body")

    def slider_s(widget, event):
        """ Extends collimator """
        # reset mesh colours
        set_colours([collimator, couchtop, body], [collimator_colour, couchtop_colour, body_colour])
        # extend snout
        value = widget.GetRepresentation().GetValue()
        collimator.SetPosition(0,value,0)
        # record new geometry
        global geometry
        geometry["SnoutExt"]=round(value,1)
        # detect collisions
        collivision(collimator,couchtop, "collimator", "couchtop")
        collivision(collimator,body, "collimator", "body")

    # create sliders
    plt.addSlider3D(
        slider_y,
        pos1=[0., -2550., -1500.],
        pos2=[0., -450., -1500.],
        xmin=-500,
        xmax=500,
        value=0,
        s=0.08,
        c="r",
        rotation=90,
        title="y position (mm)",
    )

    plt.addSlider3D(
        slider_z,
        pos1=[0., -2550., -1500.],
        pos2=[0., -2550., 500.],
        xmin=-500,
        xmax=500,
        value=0,
        s=0.08,
        c="g",
        rotation=180,
        title="z position (mm)",
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
    
    # display
    plt.show([body, gantry, collimator, couchtop, couch], axes=11, camera={"pos": (8000,0,-8000), "viewup": [-0.5,-0.5,0.5]}, bg='black', bg2='dark grey').close()

if __name__ == "__main__":
    main()
    