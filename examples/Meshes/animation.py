#!/usr/bin/env python

# have a look at this when manipulating the code:
# https://www.opengl.org/archives/resources/features/KilgardTechniques/oglpitfall/

# run with: python Meshes/animation.py


"""
Animated 3D vocal tract shapes and EMA points.

Help Menu:
h: print help menue to console
left mouse: rotate
right mouse: move around
mouse wheel: zoom
arrow up: increase speed
arrow down: decrease speed
space: pause
escape, q: quit animation
"""



import ctypes
from pathlib import Path

import pyglet
from pyglet.gl import *  # import von OpenGL symbolen
from pyglet.window import mouse
from pyglet.window import key
from pywavefront import visualization
import pywavefront
from tqdm import tqdm

WORD = 'apfelsine'


ema_path = Path(f'./Meshes/{WORD}/{WORD}-ema.txt')
mesh_path = f'./Meshes/{WORD}/{WORD}-meshes/blender_export/'
mesh_file_name_without_numeration_and_ending = WORD


rotation_x = 360.0 - 55.0
rotation_y = 5.0
zoom = -15.0
trans_x = 0.0
trans_y = 3.0

time = 0.0           # total amount of time animation is running
mesh_disp = 0.0      # amount of time current image is displayed
disp_time = 0.005      # how long should each image be displayed
disp_time_old = 0.005  # how long should each image be displayed before pause

num_img = 92  # amount of mesh files the currently displayed word consist of

mesh_counter = 0  # which mesh is currently displayed
meshes = []       # mesh-data
mesh_names = []   # names of mesh-data files
mesh = None       # currently displayed mesh data

emas = []              # stores all coordinates of ema at all timepoints
num_ema_coord = 0      # amount of ema coordinates in one timestep
num_emas = 0           # amount of emas to be displayes
emas_coord = []        # coordinates of emas at current timestep
index = 0              # index of x-value of currently displayed ema in emas[]
ema_point_size = 12.  # size of the visualized ema points


with open(ema_path, 'r') as ema_file:
    next(ema_file)
    for line in ema_file:
        single_elems = line.split()
        num_ema_coord = len(single_elems)-1  # amount of elems in one line
                                             # (first element is time and therefore not counted)
        for elem in single_elems[1:]:  # ignore first elem which tracks time
            number = float(elem)
            emas.append(number)

num_emas = num_ema_coord // 3


# store mesh names
for i in range(num_img - 1):
    mesh_names.append(mesh_path + mesh_file_name_without_numeration_and_ending + str(i) + '.obj')

print("Loading Mesh Data")

# store mesh data (takes quite a while)
for name in tqdm(mesh_names):
    meshes.append(pywavefront.Wavefront(name))

# select currently used mesh
mesh = meshes[mesh_counter]

window = pyglet.window.Window(1280, 720, resizable=True)

lightfv = ctypes.c_float * 4



@window.event
def on_resize(width, height):
    glViewport(0, 0, width, height)  # neccessary to have the scene printed in the whole window (maybe, it is better to use glOrtho)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60., float(width)/height, .1, 100.)
    glMatrixMode(GL_MODELVIEW)
    return True


@window.event
def on_draw():
    window.clear()
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()  #glLoadIdentity replaces the current matrix with the identity matrix

# Direction of axes:
#    +x: center to right
#    +y: center to top
#    +z: center to front

# Hint: If the scene is lit to bright, scale down GL_DIFFUSE of every light

    # Light source: outside the screen (front)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightfv(0.65, 0.65, 0.65, 1))
    glLightfv(GL_LIGHT0, GL_POSITION, lightfv(0.0, 0.0, 100.0, 0.0))

    # Light source: behind the screen (back)
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightfv(0.65, 0.65, 0.65, 1))
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightfv(1, 1, 1, 1))
    glLightfv(GL_LIGHT1, GL_POSITION, lightfv(0.0, 0.0, -100.0, 0.0))

    # Light source: right side of the screen
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lightfv(0.65, 0.65, 0.65, 1))
    glLightfv(GL_LIGHT2, GL_SPECULAR, lightfv(1, 1, 1, 1))
    glLightfv(GL_LIGHT2, GL_POSITION, lightfv(-100.0, 0.0, 0.0, 0.0))

    # Light source: left side of the screen
    glLightfv(GL_LIGHT3, GL_DIFFUSE, lightfv(0.65, 0.65, 0.65, 1))
    glLightfv(GL_LIGHT3, GL_SPECULAR, lightfv(1, 1, 1, 1))
    glLightfv(GL_LIGHT3, GL_POSITION, lightfv(100.0, 0.0, 0.0, 0.0))

    glEnable(GL_LIGHT0)
    glEnable(GL_LIGHT1)
    glEnable(GL_LIGHT2)
    glEnable(GL_LIGHT3)
    glEnable(GL_LIGHTING)

    # perspective and rotation influenced by mouse actions
    glTranslated(trans_x, trans_y, zoom)
    glRotatef(rotation_x, 0.0, 1.0, 0.0)
    glRotatef(rotation_y, 0.0, 0.0, -1.0)

    visualization.draw(mesh)

    # define color and size of Ema points
    glPointSize(ema_point_size)
    glColor3f(0.1, 0.3, 1.0)
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE)
    glColorMaterial(GL_FRONT, GL_DIFFUSE)
    glEnable(GL_COLOR_MATERIAL)  # if no material is enabled, ema points have default material value and are displayed greyish

    # draw all Ema Points
    glBegin(GL_POINTS)
    for ind in range(0, num_emas):
        glVertex3f(emas_coord[ind*3], emas_coord[ind*3 + 1], emas_coord[ind*3 + 2])
    glEnd()

    glDisable(GL_COLOR_MATERIAL)  # disabeling neccessary to not overwrite mesh material


@window.event
def on_mouse_scroll(x, y, scroll_x, scroll_y):
    global zoom
    zoom += 0.4 * scroll_y


@window.event
def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
    global rotation_x, rotation_y, trans_x, trans_y
    if buttons & mouse.LEFT:
        rotation_x += 1 * dx
        rotation_y += 1 * dy
    elif buttons & mouse.RIGHT:
        trans_x += 0.02 * dx
        trans_y += 0.02 * dy


@window.event
def on_key_press(symbol, modifiers):
    global disp_time, disp_time_old
    if symbol == key.UP:
        disp_time -= 0.001
        if disp_time < 0.001:
            disp_time = 0.001
            print("minimal display time per frame is 0.001.")
        print(f"speed up to {disp_time} sec/frame")
    if symbol == key.DOWN:
        if disp_time < 0.300:
            disp_time += 0.001
            print(f"speed down to {disp_time} sec/frame")
    if symbol == key.SPACE:
        if disp_time == 1000000:
            disp_time = disp_time_old
            print("ANIMATION RESUMED")
        else:
            disp_time_old = disp_time
            disp_time = 1000000
            print("ANIMATION PAUSED")
    if symbol == key.H: # print help in terminal
        print(__doc__)
    if symbol in (key.Q, key.ESCAPE):
        print("Animation closed")
        window.close()



def update(dt):
    """ Defines timesteps of the animation. """
    global rotation_x, rotation_y, time, mesh_disp, mesh_counter, mesh, \
           meshes, disp_time, num_ema_coord, num_emas, index, ema_coord
    time += dt  # keep in mind absoulte time (currently obsolete)
    mesh_disp += dt  # how long is current mesh displayed

    if mesh_disp > disp_time:  # is current mesh is displayed longer than it should?
        if mesh_counter == len(meshes) - 1:  # last mesh?
            print("  restart animation (press 'h' for help)")
            mesh_counter = 0  # start with first mesh again
        else:
            mesh_counter += 1
        mesh = meshes[mesh_counter]  # update currently used mesh
        mesh_disp = 0.0  # reset display time of this mesh

    emas_coord.clear()

    for ind in range(0, num_emas):
        index = mesh_counter * num_ema_coord + ind*3
        emas_coord.extend([emas[index], emas[index+1], emas[index+2]])


    if rotation_x > 720.0:
        rotation_x = 0.0
    if rotation_y > 720.0:
        rotation_y = 0.0


pyglet.clock.schedule(update)
if __name__ == '__main__':
    pyglet.app.run()
