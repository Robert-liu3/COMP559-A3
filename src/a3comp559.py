# YOUR NAME AND STUDENT NUMBER HERE
# Make any changes you need to this file (e.g., add UI elements to help you test)

import argparse
import numpy as np
import polyscope as ps
import polyscope.imgui as psim
from pynput.keyboard import Listener, Key
from integrators import *
from particle_system import *

# Initialize polyscope in 2D scene mode
ps.init()
ps.set_always_redraw(True)
ps.set_navigation_style("planar")
ps.set_ground_plane_mode("none")
ps.set_view_projection_mode("orthographic")

is_running = False
ui_int = 4
method = integration_methods[ui_int]

def main_display_loop():        
    global is_running, S, method, ui_int
    
    if psim.Button("Reset"): S.reset()
    do_step = psim.Button("Step")
    _, is_running	= psim.Checkbox("Run", is_running) 
    psim.TextUnformatted("Elapsed = " + str(S.elapsed))
    _, S.h = psim.SliderFloat("step size", S.h, v_min=0.001, v_max=0.1)
    _, S.stiffness = psim.SliderFloat("stiffness", S.stiffness, v_min=1, v_max=10000)
    _, S.damping = psim.SliderFloat("damping", S.damping, v_min=0, v_max=1)
    _, S.gravity = psim.SliderFloat("gravity y", S.gravity, v_min=-20, v_max=20)
    if psim.Button("Zero Gravity"): S.gravity = 0.0 

    _, ui_int = psim.InputInt("Method", ui_int,step=1 ) 
    ui_int = max(0, min(len(integration_methods)-1, ui_int)) 
    method = integration_methods[ui_int]
    psim.TextUnformatted("Current Example = " + S.name)
    psim.TextUnformatted("Current Solver = " + method.name)

    if is_running or do_step:
        S.advance_time( method )

parser = argparse.ArgumentParser(description='a3comp599')
parser.add_argument('--file', type=str, help='input json file')
args = parser.parse_args()
S = particle_system(args.file)
np.set_printoptions(precision=4) # helps interactive debugging?
ps.set_user_callback(main_display_loop)
ps.show()