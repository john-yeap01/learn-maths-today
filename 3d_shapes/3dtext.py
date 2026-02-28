from manim import *
from modules import MyThreeDScene, make_fixed  #!!!!! THIS

# ---- Use modified scene and camera system ----
# with help from this gist
# https://gist.github.com/abul4fia/1419b181e8e3410ef78e6acc25c3df94



class Counter3D(MyThreeDScene):
    def construct(self):
        # Camera orientation
        self.set_camera_orientation(
            phi=65 * DEGREES,
            theta=45 * DEGREES
        )

        axes = ThreeDAxes()
        self.add(axes)

        value = ValueTracker(0)

        label = Text("Counter:", font_size=36)

        number = always_redraw(
            lambda: Text(
                str(int(value.get_value())),
                font_size=36
            ).next_to(label, RIGHT, buff=0.3)
        )

        counter = VGroup(label, number)
        counter.to_corner(UL)

        # Make fixed AFTER construction
        make_fixed(counter)

        self.add(counter)

        self.play(
            value.animate.set_value(100),
            run_time=5,
            rate_func=linear
        )

        self.wait()

class CubeWithHUD(MyThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=65 * DEGREES, theta=45 * DEGREES)
        self.add(ThreeDAxes())

        s = ValueTracker(1.0)

        # --- 3D cube (one object, scale updated) ---
        cube = Cube(side_length=1.0, fill_opacity=0.25, fill_color=BLUE, stroke_color=BLUE_E)
        cube.move_to(ORIGIN)
        self.add(cube)

        def cube_scale_updater(m: Mobject):
            # set overall scale so cube side length == s
            target = s.get_value()
            m.scale(target / m.get_width())  # since width == side length for axis-aligned cube

        # Better: avoid compounding scale by resetting from an initial copy
        cube_base = cube.copy()

        def cube_updater(m: Mobject):
            target = s.get_value()
            m.become(cube_base.copy().scale(target).move_to(ORIGIN))

        cube.add_updater(cube_updater)

        # --- HUD (same as working step) ---
        label_len = Text("Length:", font_size=36).to_corner(UL)
        label_vol = Text("Volume:", font_size=36).next_to(label_len, DOWN, aligned_edge=LEFT, buff=0.25)

        val_len = always_redraw(
            lambda: Text(f"{s.get_value():.2f}", font_size=36)
                    .next_to(label_len, RIGHT, buff=0.35)
                    .align_to(label_len, DOWN)
        )
        val_vol = always_redraw(
            lambda: Text(f"{(s.get_value()**3):.2f}", font_size=36)
                    .next_to(label_vol, RIGHT, buff=0.35)
                    .align_to(label_vol, DOWN)
        )

        make_fixed(label_len, label_vol, val_len, val_vol)
        self.add(label_len, val_len, label_vol, val_vol)

        # Camera motion (optional)
        self.begin_ambient_camera_rotation(rate=0.3)

        # Animate
        self.play(s.animate.set_value(3.0), run_time=8, rate_func=linear)
        self.wait(5)

        # Clean up updater if you want
        cube.remove_updater(cube_updater)

#---- Text on 3D
class Text3D1(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES,theta=-45*DEGREES)
        text3d=Text("This is a 3D text").scale(2)
        self.add(axes,text3d)

class Text3D2(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES,theta=-45*DEGREES)
        text3d=Text("This is a 3D text").scale(2).set_shade_in_3d(True) 
        text3d.rotate(PI/2,axis=RIGHT)
        self.add(axes,text3d)

# wARNING: this IS A VERY HEAVY ANIMATION TO RUN
class elipsoide(ThreeDScene):
	def construct(self):

		## x = a.cos(u)cos(v)
		## y = b.sen(u)cos(v)
		## z = c.sen(v)
		##np.array([x, y, z])
		##Sphere(radius = 2, color = , ....)


		ejes = ThreeDAxes(x_range = [-6, 6, 1],
			              y_range = [-6, 6, 1],
			              z_range = [-6, 6, 1],
			              x_length = 7,
			              y_length = 7,
			              z_length = 7).add_coordinates()

		a = ValueTracker(1)
		b = ValueTracker(1)
		c = ValueTracker(1)

		elipsoide = always_redraw(
			lambda: Surface(lambda u, v: ejes.c2p(*np.array([a.get_value()*np.cos(u)*np.cos(v), 
				            b.get_value()*np.sin(u)*np.cos(v), 
				            c.get_value()*np.sin(v)])),
				            u_range = [0, PI],
				            v_range = [0, 2*PI]))


		self.set_camera_orientation(phi = 65*DEGREES, theta = 60*DEGREES)
		self.play(Create(ejes), run_time = 2)
		self.wait(2)

		self.play(Create(elipsoide), run_time = 2)
		self.wait(2)
		self.play(a.animate.set_value(5), rate_func = there_and_back, run_time = 8) 
		self.wait(2)
		self.play(b.animate.set_value(5), rate_func = there_and_back, run_time = 8) 
		self.wait(2)
		self.play(c.animate.set_value(5), rate_func = there_and_back, run_time = 8) 
		self.wait(2)