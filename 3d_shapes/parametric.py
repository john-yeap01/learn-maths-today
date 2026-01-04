from manim import *
import numpy as np

#----- Surfaces
class SurfacesAnimation(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        cylinder = Surface(
            lambda u, v: np.array([
                np.cos(TAU * v),
                np.sin(TAU * v),
                2 * (1 - u)
            ]),
            resolution=(6, 32)).fade(0.5) #Resolution of the surfaces

        paraboloid = Surface(
            lambda u, v: np.array([
                np.cos(v) * u,
                np.sin(v) * u,
                u**2
            ]),
            u_range=(0, 2),       # pick what you want
            v_range=(0, TAU),
            checkerboard_colors=[PURPLE_D, PURPLE_E],
            resolution=(10, 32),
        ).scale(2)

        para_hyp = Surface(
            lambda u, v: np.array([u, v, u**2 - v**2]),
            u_range=(-2, 2),
            v_range=(-2, 2),
            checkerboard_colors=[BLUE_D, BLUE_E],
            resolution=(15, 32),
        )

        cone = Surface(
            lambda u, v: np.array([u*np.cos(v), u*np.sin(v), u]),
            u_range=(0, 2),        # radius should usually be >= 0
            v_range=(0, TAU),
            checkerboard_colors=[GREEN_D, GREEN_E],
            resolution=(15, 32),
        )

        hip_one_side = Surface(
            lambda u, v: np.array([
                np.cosh(u)*np.cos(v),
                np.cosh(u)*np.sin(v),
                np.sinh(u)
            ]),
            u_range=(-2, 2),
            v_range=(0, TAU),
            checkerboard_colors=[YELLOW_D, YELLOW_E],
            resolution=(15, 32),
        )

        ellipsoid = Surface(
            lambda u, v: np.array([
                1*np.cos(u)*np.cos(v),
                2*np.cos(u)*np.sin(v),
                0.5*np.sin(u)
            ]),
            u_range=(-PI/2, PI/2),
            v_range=(0, TAU),
            checkerboard_colors=[TEAL_D, TEAL_E],
            resolution=(15, 32),
        ).scale(2)
        
        sphere = Surface(
            lambda u, v: np.array([
                1.5*np.cos(u)*np.cos(v),
                1.5*np.cos(u)*np.sin(v),
                1.5*np.sin(u)
            ]),
            u_range=(-PI/2, PI/2),
            v_range=(0, TAU),
            checkerboard_colors=[RED_D, RED_E],
            resolution=(15, 32),
        ).scale(2)


        self.set_camera_orientation(phi=75 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)


        self.add(axes)
        self.play(Write(sphere))
        self.wait()
        self.play(ReplacementTransform(sphere,ellipsoid))
        self.wait()
        self.play(ReplacementTransform(ellipsoid,cone))
        self.wait()
        self.play(ReplacementTransform(cone,hip_one_side))
        self.wait()
        self.play(ReplacementTransform(hip_one_side,para_hyp))
        self.wait()
        self.play(ReplacementTransform(para_hyp,paraboloid))
        self.wait()
        self.play(ReplacementTransform(paraboloid,cylinder))
        self.wait()
        self.play(FadeOut(cylinder))