# 2d animation for filling a boundary of a 2d shape .
from manim import *
import numpy as np
from numpy.typing import NDArray



class FillShader (ThreeDScene):

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)

        square = Square(color= BLUE, fill_opacity=0)

        self.play(FadeIn(axes, square))
        # self.add(axes)


        self.play(square.animate.set_fill(BLUE, opacity=1), run_time=2)
    
        self.stop_ambient_camera_rotation()
        



    