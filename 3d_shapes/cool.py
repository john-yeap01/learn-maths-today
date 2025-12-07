#contains cool geometric animatiosn for b roll


from manim import *

class CubeToSphere(ThreeDScene):
    def construct(self):
        # ---------- Camera: simple fixed 3D angle ----------
        self.set_camera_orientation(phi=65*DEGREES, theta=45*DEGREES)

        # ---------- 3D Objects ----------
        cube = Cube(
            side_length=2,
            fill_opacity=0.8,
            fill_color=BLUE,
            stroke_width=1,
            stroke_color=WHITE,
        )

        sphere = Sphere(
            radius=1.2,
            resolution=(32, 64),
            fill_opacity=0.8,
            fill_color=GREEN,
            stroke_width=0.8,
            stroke_color=WHITE,
        )
        sphere.move_to(cube.get_center())

        # ---------- 1. Show cube ----------
        self.play(FadeIn(cube, scale=0.5), run_time=1.5)

        # Optional small spin
        self.play(
            Rotate(cube, angle=PI/3, axis=UP),
            Rotate(cube, angle=PI/5, axis=RIGHT),
            run_time=2,
        )

        self.wait(0.5)

        # ---------- 2. Cube → Sphere transformation ----------
        self.play(
            Transform(cube, sphere),
            run_time=3.0,
            rate_func=smooth,
        )

        self.wait(1)

        # ---------- 3. Outro ----------
        self.play(FadeOut(cube), run_time=1.5)
        self.wait()




class SphereToPyramid(ThreeDScene):
    def construct(self):
        # Camera angle (fixed)
        self.set_camera_orientation(phi=65 * DEGREES, theta=45 * DEGREES)

        # ---------- Sphere ----------
        sphere = Sphere(
            radius=1.3,
            resolution=(24, 48),  # lower = faster
            fill_color=GREEN,
            fill_opacity=0.85,
            stroke_color=WHITE,
            stroke_width=0.6,
        )

        # ---------- Pyramid (square base) ----------
        vertex_coords = [
            [1, 1, 0],
            [1, -1, 0],
            [-1, -1, 0],
            [-1, 1, 0],
            [0, 0, 2],
        ]
        faces = [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
            [0, 1, 2, 3],  # base
        ]

        pyramid = Polyhedron(
            vertex_coords,
            faces,
            faces_config={
                "fill_color": RED,
                "fill_opacity": 0.9,
                "stroke_color": WHITE,
                "stroke_width": 0.8,
            },
        ).scale(0.9)

        pyramid.move_to(sphere.get_center())

        # ---------- Animation ----------

        # Show sphere
        self.play(FadeIn(sphere, scale=0.5), run_time=1.0)
        self.wait(0.2)

        # Crossfade: Sphere → Pyramid (no Transform)
        self.play(
            FadeOut(sphere),
            FadeIn(pyramid),
            run_time=1.5,
        )

        self.wait(0.5)

        # Fade out everything
        self.play(FadeOut(pyramid), run_time=1.0)
        self.wait()


class CubeSpherePyramid(ThreeDScene):
    def construct(self):
        # ---------- Camera ----------
        self.set_camera_orientation(phi=65*DEGREES, theta=45*DEGREES)

        # ---------- Shapes ----------
        # Cube
        cube = Cube(
            side_length=2,
            fill_opacity=0.8,
            fill_color=BLUE,
            stroke_width=1,
            stroke_color=WHITE,
        )

        # Sphere
        sphere = Sphere(
            radius=1.2,
            resolution=(32, 64),
            fill_opacity=0.8,
            fill_color=GREEN,
            stroke_width=0.8,
            stroke_color=WHITE,
        )
        sphere.move_to(cube.get_center())

        # Pyramid (square base pyramid)
        vertex_coords = [
            [1, 1, 0],     # Base corners
            [1, -1, 0],
            [-1, -1, 0],
            [-1, 1, 0],
            [0, 0, 2],     # Apex
        ]
        faces = [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
            [0, 1, 2, 3],  # Base
        ]

        pyramid = Polyhedron(
            vertex_coords,
            faces,
            faces_config={
                "fill_color": RED,
                "fill_opacity": 0.85,
                "stroke_color": WHITE,
                "stroke_width": 1,
            },
        ).scale(0.9)
        pyramid.move_to(cube.get_center())

        # ---------- 1. Show cube ----------
        self.play(FadeIn(cube, scale=0.5), run_time=1.5)
        self.play(
            Rotate(cube, angle=PI/3, axis=UP),
            Rotate(cube, angle=PI/5, axis=RIGHT),
            run_time=2,
        )
        self.wait(0.4)

        # ---------- 2. Cube → Sphere ----------
        self.play(
            Transform(cube, sphere),
            run_time=2.5,
            rate_func=smooth,
        )
        self.wait(0.4)

        # ---------- 3. Sphere → Pyramid ----------
        self.play(
            Transform(cube, pyramid),
            run_time=2.5,
            rate_func=smooth,
        )
        self.wait(1)

        # ---------- 4. Outro ----------
        self.play(FadeOut(cube), run_time=1.5)
        self.wait()



class TwoDToThreeD(ThreeDScene):
    def construct(self):
        # --- 0. Start the camera in a "2D-like" orientation ---
        # Camera looking straight at the XY-plane (flat view)
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES)

        # --- 1. 2D WORLD SETUP ---
        plane = NumberPlane(
            x_range=(-4, 4, 1),
            y_range=(-3, 3, 1),
            background_line_style={
                "stroke_color": GREY_D,
                "stroke_opacity": 0.8,
                "stroke_width": 1,
            },
        )

        square_2d = Square(side_length=2, color=BLUE, fill_opacity=0.4)
        circle_2d = Circle(radius=1, color=YELLOW, fill_opacity=0.4).shift(RIGHT * 2.5)

        label_2d = Text("2D View", font_size=36).to_edge(UP)

        self.add(plane)
        self.play(
            FadeIn(square_2d, shift=DOWN * 0.5),
            FadeIn(circle_2d, shift=UP * 0.5),
            FadeIn(label_2d),
            run_time=1.5,
        )
        self.wait(0.7)

        # Maybe a little 2D motion
        self.play(
            square_2d.animate.shift(LEFT * 1),
            circle_2d.animate.shift(LEFT * 1),
            run_time=1.2,
        )
        self.wait(0.5)

        # --- 2. PREP 3D OBJECTS (hidden at first) ---
        cube_3d = Cube(
            side_length=2,
            fill_color=BLUE_E,
            fill_opacity=0.9,
            stroke_width=0.8,
            stroke_color=WHITE,
        )
        cube_3d.move_to(square_2d.get_center())

        sphere_3d = Sphere(
            radius=1,
            resolution=(32, 64),
            fill_color=YELLOW_E,
            fill_opacity=0.9,
            stroke_width=0.6,
            stroke_color=WHITE,
        )
        sphere_3d.move_to(circle_2d.get_center())

        axes_3d = ThreeDAxes(
            x_range=(-4, 4, 1),
            y_range=(-3, 3, 1),
            z_range=(-3, 3, 1),
        )

        label_3d = Text("3D View", font_size=36).to_edge(UP).set_opacity(0)

        # --- 3. CAMERA TILT + SHAPE MORPH (the "dimension jump") ---
        # --- 3. CAMERA TILT + SHAPE MORPH (the "dimension jump") ---
        self.move_camera(
            phi=65 * DEGREES,
            theta=45 * DEGREES,
            run_time=3.0,
            added_anims=[
                FadeOut(label_2d),
                FadeIn(label_3d),
                Transform(square_2d, cube_3d),
                Transform(circle_2d, sphere_3d),
                FadeIn(axes_3d),
            ],
        )



        # At this point, `square_2d` is actually the cube, and `circle_2d` is the sphere.
        self.wait(0.8)

        # --- 4. Let the 3D world show off a bit ---
        self.begin_ambient_camera_rotation(rate=0.15)

        self.play(
            square_2d.animate.rotate(PI / 3, axis=OUT),
            circle_2d.animate.rotate(PI / 2, axis=UP),
            run_time=2.5,
        )
        self.wait(1.5)

        # --- 5. Outro ---
        self.play(
            FadeOut(square_2d),
            FadeOut(circle_2d),
            FadeOut(axes_3d),
            FadeOut(plane),
            FadeOut(label_3d),
            run_time=1.5,
        )
        self.wait()
