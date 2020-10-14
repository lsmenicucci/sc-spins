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
attribute float rx;
attribute float ry;
attribute float sx;
attribute float sy;
varying float v_rx;
varying float v_ry;
varying float v_sx;
varying float v_sy;
void main()
{
    v_rx = rx;
    v_ry = ry;
    v_sx = sx;
    v_sy = sy;
    gl_PointSize = 15.0;
    gl_Position = vec4(vec2(rx, ry), 0.0, 1.0);
} """

fragment = """
#include "math/constants.glsl"
#include "arrows/arrows.glsl"
#include "antialias/antialias.glsl"

varying float v_rx;
varying float v_ry;
varying float v_sx;
varying float v_sy;
void main()
{ 
    const float M_PI = 3.14159265358979323846;
    const float size = 15.0;
    const float linewidth = 4.0;
    const float antialias =  1.0;
    const float body = 15;

    vec2 texcoord = gl_FragCoord.xy;
    vec2 center = (2*gl_PointCoord.xy - 1)*size + texcoord.xy;

    texcoord -= center;
    float theta = M_PI - atan(v_sy, v_sx);
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    texcoord = vec2(cos_theta*texcoord.x - sin_theta*texcoord.y,
                    sin_theta*texcoord.x + cos_theta*texcoord.y);

    float d = arrow_stealth(texcoord, body, 0.25*body, linewidth, antialias);
    gl_FragColor = filled(d, linewidth, antialias, vec4(0,0,0,1));
} """

class InteractiveView(threading.Thread):
    def __init__(self, spintronics):
        threading.Thread.__init__(self)
        # Get the maximum length
        max_dim = np.max([np.max(np.abs(spintronics.rx)), np.max(np.abs(spintronics.ry))])

        print(f"max dim: {max_dim}")

        self.rx = spintronics.rx/max_dim * 0.9
        self.ry = spintronics.ry/max_dim * 0.9
        self.sx = spintronics.sx
        self.sy = spintronics.sy

        print(f"max rx: {np.max(self.rx)}")

    def run(self):
        window = app.Window(512, 512, color=(1,1,1,1))
        points = gloo.Program(vertex, fragment)
        points["rx"] = self.rx
        points["ry"] = self.ry
        points["sx"] = self.sx
        points["sy"] = self.sy

        @window.event
        def on_resize(width, height):
            points["resolution"] = width, height
            
        @window.event
        def on_draw(dt):
            points["sx"] = self.sx
            points["sy"] = self.sy
            window.clear()
            points.draw(gl.GL_POINTS)

        app.run()
