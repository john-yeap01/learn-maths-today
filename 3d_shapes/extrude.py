from manim import *
import numpy as np
from numpy.typing import NDArray

class ExtrudePrism(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)

        axes = ThreeDAxes()
        self.add(axes)

        h = ValueTracker(0)  # extrusion height

        base = Square(side_length=2).set_fill(RED, 0.6).set_stroke(RED)

        # Top face follows base but lifted by OUT*h
        top = always_redraw(
            lambda: base.copy()
                    .shift(OUT * h.get_value())
                    .set_fill(RED, 0.6)
                    .set_stroke(RED)
        )

        # Side faces
        def make_side(i):
            return always_redraw(lambda: (
                lambda b, t: Polygon(
                    b[i], b[(i+1) % 4], t[(i+1) % 4], t[i]
                ).set_fill(RED, 0.3).set_stroke(RED)
            )(
                base.get_vertices(),
                [v + OUT * h.get_value() for v in base.get_vertices()]
            ))

        sides = VGroup(*[make_side(i) for i in range(4)])

        # --------------------------------------------------
        # Volume text (NO LaTeX)
        # Volume = base area × height = 4 × h
        # --------------------------------------------------
        volume_text = always_redraw(
            lambda: Text(
                f"Volume = 4 × {h.get_value():.2f} = {4*h.get_value():.2f}",
                font_size=28
            )
            .to_corner(UL)
        )

        area_text = always_redraw(
            lambda: Text("Area = 4", font_size=28)
                # Lay flat on XY plane
                .rotate(-PI/2, axis=RIGHT)
                # Un-mirror front/back
                .rotate(PI, axis=OUT)
                # FINAL FIX: correct left-right flip
                .rotate(PI, axis=LEFT)
        )
        
        height_text = Text(
            "Height = 2",
            font_size=28
        ).to_corner(UL).shift(DOWN * 0.8)
                
                
        

        self.add_fixed_in_frame_mobjects(volume_text)

        self.add(area_text)

        # IMPORTANT: do NOT also add them with self.add(...)
        self.add(base, top, sides)

        self.wait(4)


        # Animate extrusion
        self.play(h.animate.set_value(2), run_time=2, rate_func=smooth)
        # Fade in height text AFTER extrusion completes
        self.wait(3)

        self.add_fixed_in_frame_mobjects(height_text)

        self.play(FadeIn(height_text, shift=DOWN * 0.2), run_time=1)

        




# gets the points of the vertices for nice polygons and the outline for not so nice ones

def extract_points(obj):
    if hasattr(obj, "get_vertices"):
        vertices = obj.get_vertices()
        
    return vertices

# creates the side planes of the object to rise up along with the top

def build_sides(vertices, height):
    b = vertices
    hv = height.get_value()
    t = [p + OUT * hv for p in b]
    n = len(b)
    quads = []
    for i in range(n):
        j = (i + 1) % n
        q = Polygon(b[i], b[j], t[j], t[i])\
            .set_fill(BLUE, 0.35).set_stroke(WHITE, 1)
        quads.append(q)
    return VGroup(*quads)

class Generalised(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        axes = ThreeDAxes()
        self.add(axes)

        self.begin_ambient_camera_rotation(rate=0.2)

        base = Star(5, outer_radius=1, inner_radius=0.5)\
               .set_fill(BLUE, 0.6).set_stroke(WHITE, 1)

        vertices = extract_points(base)
        self.add(base)

        h = ValueTracker(0.01)

        # TOP FACE
        top = always_redraw(
            lambda: Polygon(*(p + OUT * h.get_value() for p in vertices))
                    .set_fill(BLUE, 0.6).set_stroke(WHITE, 1)
        )
        self.add(top)

        # SIDE FACES

        sides = always_redraw(lambda: build_sides(vertices, h))
        self.add(sides)

        self.play(h.animate.set_value(2.5), run_time=2)

        self.wait(2)
        self.stop_ambient_camera_rotation()



class HexPrismCrossSection(ThreeDScene):
    def construct(self):
        # --- Camera ---
        self.set_camera_orientation(phi=65 * DEGREES, theta=35 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.15)

        axes = ThreeDAxes()
        self.add(axes)

        # --- Prism parameters ---
        R = 2.0      # radius of hexagon (distance from center to vertex)
        H = 3.0      # total height of prism

        z_top =  H / 2
        z_bot = -H / 2

        # --- Hexagon vertices (in XY), then lift to z ---
        angles = np.linspace(0, TAU, 6, endpoint=False)
        xy = np.c_[R * np.cos(angles), R * np.sin(angles)]

        top_pts = [np.array([x, y, z_top]) for x, y in xy]
        bot_pts = [np.array([x, y, z_bot]) for x, y in xy]

        # --- Faces ---
        top_face = Polygon(*top_pts, color=BLUE, fill_opacity=0.35, stroke_width=2)
        bot_face = Polygon(*bot_pts, color=BLUE, fill_opacity=0.20, stroke_width=2)

        side_faces = VGroup()
        for i in range(6):
            j = (i + 1) % 6
            quad = Polygon(
                bot_pts[i], bot_pts[j], top_pts[j], top_pts[i],
                color=TEAL, fill_opacity=0.20, stroke_width=2
            )
            side_faces.add(quad)

        prism = VGroup(top_face, bot_face, side_faces)

        # --- Cross-sectional plane in the middle (z = 0) ---
        mid_pts = [np.array([x, y, 0.0]) for x, y in xy]
        cross_section1 = Polygon(
            *mid_pts,
            color=YELLOW,
            fill_color=YELLOW,
            fill_opacity=0.55,
            stroke_width=4
        )

        # --- Cross-sectional plane in the middle (z = 0) ---
        above_pts = [np.array([x, y, H/4]) for x, y in xy]
        cross_section2 = Polygon(
            *above_pts,
            color=GREEN,
            fill_color=GREEN,
            fill_opacity=0.55,
            stroke_width=4
        )

        # --- Cross-sectional plane in the middle (z = 0) ---
        below_pts = [np.array([x, y, -H/4]) for x, y in xy]
        cross_section3 = Polygon(
            *below_pts,
            color=BLUE,
            fill_color=BLUE,
            fill_opacity=0.55,
            stroke_width=4
        )

        # Optional: emphasize its boundary
        cross_outline1 = Polygon(*mid_pts, color=YELLOW, stroke_width=6)
        cross_outline2 = Polygon(*above_pts, color=YELLOW, stroke_width=6)
        cross_outline3 = Polygon(*below_pts, color=YELLOW, stroke_width=6)


        # --- Show ---
        self.play(FadeIn(prism), run_time=1.2)
        self.play(FadeIn(cross_section1), Create(cross_outline1), run_time=3.0)
        self.play(FadeIn(cross_section2), Create(cross_outline2), run_time=3.0)
        self.play(FadeIn(cross_section3), Create(cross_outline3), run_time=3.0)


        # Hold and rotate a bit
        self.wait(4)

        self.stop_ambient_camera_rotation()
        self.wait(0.5)
class TrapezoidPrismCrossSection_OnRectFace(ThreeDScene):
    def construct(self):
        # --- Camera ---
        self.set_camera_orientation(phi=65 * DEGREES, theta=35 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.15)

        axes = ThreeDAxes()
        self.add(axes)

        # --- Prism parameters ---
        H = 3.0
        z_top =  H / 2
        z_bot = -H / 2

        # --- Trapezoid vertices in XY (counterclockwise) ---
        xy = [
            np.array([-2.2, -1.2, 0.0]),
            np.array([ 2.2, -1.2, 0.0]),
            np.array([ 1.3,  1.2, 0.0]),
            np.array([-1.3,  1.2, 0.0]),
        ]

        top_pts = [p + np.array([0, 0, z_top]) for p in xy]
        bot_pts = [p + np.array([0, 0, z_bot]) for p in xy]

        # --- Faces ---
        top_face = Polygon(*top_pts, color=BLUE, fill_opacity=0.35, stroke_width=2)
        bot_face = Polygon(*bot_pts, color=BLUE, fill_opacity=0.20, stroke_width=2)

        side_faces = VGroup()
        n = len(xy)
        for i in range(n):
            j = (i + 1) % n
            quad = Polygon(
                bot_pts[i], bot_pts[j], top_pts[j], top_pts[i],
                color=TEAL, fill_opacity=0.20, stroke_width=2
            )
            side_faces.add(quad)

        prism = VGroup(top_face, bot_face, side_faces)

        # --- Cross sections (same trapezoid, different z heights) ---
        def trapezoid_at_z(z, color, fill_opacity=0.55, stroke_width=4):
            pts = [p + np.array([0, 0, z]) for p in xy]
            poly = Polygon(
                *pts,
                color=color,
                fill_color=color,
                fill_opacity=fill_opacity,
                stroke_width=stroke_width
            )
            outline = Polygon(*pts, color=YELLOW, stroke_width=6)
            return poly, outline

        cross1, outline1 = trapezoid_at_z(0.0, YELLOW)
        cross2, outline2 = trapezoid_at_z(H/4, GREEN)
        cross3, outline3 = trapezoid_at_z(-H/4, BLUE)

        # ------------------------------------------------------------
        # ORIENT IT SO A RECTANGULAR SIDE FACE IS ON THE GROUND PLANE
        # ------------------------------------------------------------
        whole = VGroup(prism, cross1, outline1, cross2, outline2, cross3, outline3)

        # Rotate 90° about the x-axis (RIGHT) so the prism lies on its side
        whole.rotate(90 * DEGREES, axis=RIGHT)

        # Shift up so the lowest points touch z = 0 (the "ground" plane)
        min_z = whole.get_all_points()[:, 2].min()
        whole.shift(OUT * (-min_z))

        # --- Animate ---
        self.play(FadeIn(prism), run_time=1.2)
        self.play(FadeIn(cross1), Create(outline1), run_time=2.0)
        self.play(FadeIn(cross2), Create(outline2), run_time=2.0)
        self.play(FadeIn(cross3), Create(outline3), run_time=2.0)

        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait(0.5)