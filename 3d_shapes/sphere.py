from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
from manim import Sector, ValueTracker
from manim import Text
import numpy as np

class Basic (ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.add(axes)

        sphere = Sphere (radius=1, color=BLUE)
        self.add(sphere)

        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)


        # EITHER rotate the objects...
        # self.play(Rotate(VGroup(axes, sphere), angle=PI/4, axis=OUT, run_time=2))

        self.play(Rotate(VGroup(axes, sphere), angle=PI/4, axis=OUT, run_time=2))
        self.wait(2)
        # ...OR rotate the camera (more natural for 3D)
        self.begin_ambient_camera_rotation(rate=0.2)
        self.wait(2)
        self.stop_ambient_camera_rotation()

        self.play(FadeOut(sphere))


        group = VGroup(*[Dot().shift(RIGHT*i) for i in range(5)])
        group.arrange(RIGHT, buff=0.5)
        self.add(group)



