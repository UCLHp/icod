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

    # load machine stl
    gantry = Mesh(fileGantry).c('lb')
    collimator = Mesh(fileCollimator).c('dark blue')
    couchtop = Mesh(fileCouch).c('lb')
    couch = Mesh(fileTable).c('dark blue')

    # load patient stl
    body = Mesh(patient_dir).c('light green')

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

    # collision detection functions
    global collision_count
    collision_count = 0 # collision counter

    def collivision(a,b):
        """ Collision detection between two objects (a and b).
        
        Parameters:
        a (vedo mesh object)
        b (vedo mesh object)
        
        Returns:
        Running total of collisions printed to terminal
        """
        global collision_count
        a = a.clean().polydata()
        b = b.clean().polydata()
        transform0 = vtk.vtkTransform()
        transform1 = vtk.vtkTransform()
        collide = vtk.vtkCollisionDetectionFilter()
        collide.SetInputData(0, a)
        collide.SetTransform(0, transform0)
        collide.SetInputData(1, b)
        collide.SetTransform(1, transform1)
        collide.SetBoxTolerance(0.0)
        collide.SetCellTolerance(0.0)
        collide.SetNumberOfCellsPerNode(2)
        collide.SetCollisionModeToFirstContact()
        collide.GenerateScalarsOn()
        collide.Update()
        if collide.GetNumberOfContacts() > 0:
            collision_count +=1
            print("Collision!"+str(collision_count))

    # slider callback functions
    global gantry_angle
    gantry_angle = 0

    global snout_extension
    snout_extension = 0

    def slider_y(widget, event):
        """ Moves patient ANT-POST """
        value = widget.GetRepresentation().GetValue()
        body.y(value)  # set patient y position

    def slider_z(widget, event):
        """ Moves patient SUP-INF """
        value = widget.GetRepresentation().GetValue()
        body.z(value)  # set patient z position

    def slider_g(widget, event):
        """ Rotates gantry and collimator """
        value = widget.GetRepresentation().GetValue()
        T = vtk.vtkTransform()
        T.RotateZ(value)
        gantry.SetUserTransform(T)
        collimator.SetUserTransform(T)  # set gantry angle
        collivision(collimator,couchtop)

    def slider_s(widget, event):
        """ Extends collimator """
        value = widget.GetRepresentation().GetValue()
        collimator.SetPosition(0,value,0)
        collivision(collimator,couchtop)

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
    