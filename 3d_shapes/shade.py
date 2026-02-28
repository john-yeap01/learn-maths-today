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
        

class SquareColorChange(ThreeDScene):
    def construct(self):
        # Camera angled so XY plane is visible
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)

        axes = ThreeDAxes()
        self.add(axes)

        # Square on XY plane (normal is OUT)
        square = Square(side_length=2)
        square.set_fill(RED, opacity=0.8)
        square.set_stroke(RED)
        square.move_to(ORIGIN)

        self.add(square)
        self.wait(1)

        # Animate color change
        self.play(
            square.animate.set_fill(BLUE, opacity=0.8).set_stroke(BLUE),
            run_time=2
        )

        self.wait(1)

class CubeFaceSelection(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=65*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 1.5  # cube side length
        half = s / 2

        base_color = GREY_B
        select_color = YELLOW

        # --------------------------------------------------
        # Define 6 faces manually so we can animate them
        # --------------------------------------------------
        faces = [
            # +Z (top)
            Square(s).shift(OUT * half),

            # -Z (bottom)
            Square(s).shift(IN * half),

            # +X (right)
            Square(s).rotate(PI/2, axis=UP).shift(RIGHT * half),

            # -X (left)
            Square(s).rotate(PI/2, axis=UP).shift(LEFT * half),

            # +Y (back)
            Square(s).rotate(PI/2, axis=RIGHT).shift(UP * half),

            # -Y (front)
            Square(s).rotate(PI/2, axis=RIGHT).shift(DOWN * half),
        ]

        for f in faces:
            f.set_fill(base_color, 0.6)
            f.set_stroke(WHITE, 1)

        cube = VGroup(*faces)
        self.add(cube)

        self.wait(0.5)

        # --------------------------------------------------
        # Selection animation: one face at a time
        # --------------------------------------------------
        for face in faces:
            self.play(
                face.animate.set_fill(select_color, 0.9),
                run_time=0.6
            )
            self.wait(0.3)
            self.play(
                face.animate.set_fill(base_color, 0.6),
                run_time=0.4
            )

        self.wait(1)

    

class TriangularPrismFaceSelection(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=65*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        base_color   = GREY_B
        select_color = YELLOW

        # Prism parameters
        L = 2.0      # length of prism along Z
        half = L/2

        # Equilateral triangle in XY plane, centered
        r = 1.2  # "radius" from center to vertex
        angles = [90, 210, 330]  # degrees
        tri_xy = np.array([
            [r*np.cos(a*DEGREES), r*np.sin(a*DEGREES), 0.0]
            for a in angles
        ])

        # Front and back triangles (caps)
        tri_front_pts = [p + IN * half for p in tri_xy]   # z = -half
        tri_back_pts  = [p + OUT * half for p in tri_xy]  # z = +half

        cap_front = Polygon(*tri_front_pts)
        cap_back  = Polygon(*tri_back_pts)

        # 3 rectangular side faces (as quads)
        sides = []
        for i in range(3):
            j = (i + 1) % 3
            p0f, p1f = tri_front_pts[i], tri_front_pts[j]
            p0b, p1b = tri_back_pts[i],  tri_back_pts[j]
            side = Polygon(p0f, p1f, p1b, p0b)
            sides.append(side)

        faces = [cap_front, cap_back, *sides]

        # Style
        for f in faces:
            f.set_fill(base_color, 0.6)
            f.set_stroke(WHITE, 1)

        prism = VGroup(*faces)
        self.add(prism)
        self.wait(0.5)

        # Selection animation: one face at a time
        for face in faces:
            self.play(face.animate.set_fill(select_color, 0.9), run_time=0.6)
            self.wait(0.3)
            self.play(face.animate.set_fill(base_color, 0.6), run_time=0.4)

        self.wait(1)


class CylinderFaceSelection(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=65*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        base_color   = GREY_B
        select_color = YELLOW

        r = 1.2
        h = 2.0
        half = h / 2

        # --------------------------------------------------
        # Curved side (Manim cylinder)
        # --------------------------------------------------
        side = Cylinder(radius=r, height=h, direction=OUT)
        side.set_fill(base_color, 0.6).set_stroke(WHITE, 1)

        # --------------------------------------------------
        # Top & bottom caps (manual polygons)
        # --------------------------------------------------
        def make_cap(z, n=48):
            pts = [
                np.array([
                    r * np.cos(2 * PI * k / n),
                    r * np.sin(2 * PI * k / n),
                    z
                ])
                for k in range(n)
            ]
            return Polygon(*pts)

        top = make_cap(+half)
        bottom = make_cap(-half)

        for cap in (top, bottom):
            cap.set_fill(base_color, 0.6)
            cap.set_stroke(WHITE, 1)

        self.add(side, top, bottom)
        self.wait(0.5)

        # --------------------------------------------------
        # Selection animation
        # --------------------------------------------------
        parts = [side, top, bottom]

        for part in parts:
            self.play(part.animate.set_fill(select_color, 0.9), run_time=0.6)
            self.wait(0.3)
            self.play(part.animate.set_fill(base_color, 0.6), run_time=0.4)

        self.wait(1)