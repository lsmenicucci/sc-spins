# -----------------------------------------------------------------------------
# Copyright (c) 2009-2016 Nicolas P. Rougier. All rights reserved.
# Distributed under the (new) BSD License.
# -----------------------------------------------------------------------------
import threading
from time import sleep
import numpy as np
from glumpy import app, gl, glm, gloo

vertex = """
uniform vec2 resolution;

attribute vec3 color;

attribute float rx;
attribute float ry;

attribute float sx;
attribute float sy;
attribute float sz;

varying vec3 v_color;

varying float v_rx;
varying float v_ry;

varying float v_sx;
varying float v_sy;
varying float v_sz;
void main()
{
    v_rx = rx;
    v_ry = ry;

    v_sx = sx;
    v_sy = sy;
    v_sz = sz;

    gl_PointSize = 20.0;
    gl_Position = vec4(vec2(rx, ry), 0.0, 1.0);

    v_color = color;
} """

fragment = """
#include "math/constants.glsl"
#include "arrows/arrows.glsl"
#include "antialias/antialias.glsl"

varying vec3 v_color;

varying float v_rx;
varying float v_ry;

varying float v_sx;
varying float v_sy;
varying float v_sz;
void main()
{ 
    const float M_PI = 3.14159265358979323846;
    const float size = 15.0;
    const float linewidth = 4.0;
    const float antialias =  1.0;
    const float body = 20;

    float intensity = sqrt(v_sx*v_sx + v_sy*v_sy + v_sz*v_sz)/0.5;
    
    vec2 texcoord = gl_FragCoord.xy;
    vec2 center = (2*gl_PointCoord.xy - 1)*size + texcoord.xy;

    texcoord -= center;
    float theta = M_PI - atan(v_sy, v_sx);
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    texcoord = vec2(cos_theta*texcoord.x - sin_theta*texcoord.y,
                    sin_theta*texcoord.x + cos_theta*texcoord.y);

    float d = arrow_stealth(texcoord, body, 0.25*body, linewidth, antialias);
    gl_FragColor = filled(d, linewidth, antialias, vec4(v_color, intensity));
} """

class InteractiveView(threading.Thread):
    def __init__(self, spintronics, record = False):
        threading.Thread.__init__(self)
        # Get the maximum length
        max_dim = np.max([
            np.max(spintronics.rx) - np.min(spintronics.rx), 
            np.max(spintronics.ry) - np.min(spintronics.ry)
            ])

        center_x = np.mean(spintronics.rx)
        center_y = np.mean(spintronics.ry)

        self.rx = (spintronics.rx - center_x)/max_dim * 1.8
        self.ry = (spintronics.ry - center_y)/max_dim * 1.8
        self.sx = spintronics.sx
        self.sy = spintronics.sy
        self.sz = spintronics.sz
        self.spin = spintronics.spin

        self.record = record

    def run(self):
        self.window = app.Window(512, 512, color=(1,1,1,1))
        points = gloo.Program(vertex, fragment)

        spin_types = set(self.spin)
        colors = [(0, 0, 0), (86/255, 180/255, 233/255)]
        color_map = dict( zip(spin_types, colors) )


        points["rx"] = self.rx
        points["ry"] = self.ry

        points["sx"] = self.sx
        points["sy"] = self.sy
        points["sz"] = self.sz

        points["color"] = [color_map[spin] for spin in self.spin]   

        @self.window.event
        def on_resize(width, height):
            points["resolution"] = width, height
            
        @self.window.event
        def on_draw(dt):
            points["sx"] = self.sx
            points["sy"] = self.sy
            points["sz"] = self.sz

            self.window.clear()
            points.draw(gl.GL_POINTS)

        # if (self.record == True):
        #     with recordWindow(self.window, 'test.mp4', fps = 60):
        #         app.run(framerate = 60, duration = 10.0)
        # else:
        app.run()
