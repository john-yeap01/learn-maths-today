from manim import *
import numpy as np

class ExtrudePrism(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)

        axes = ThreeDAxes()
        self.add(axes)

        h = ValueTracker(0.0)  # extrusion height

        base = Square(side_length=2).set_fill(RED, 0.6).set_stroke(RED)

        # Top face follows base but lifted by OUT*h
        top = always_redraw(
            lambda: base.copy()
                    .shift(OUT * h.get_value())
                    .set_fill(RED, 0.6).set_stroke(RED)
        )

        # Side faces are quads built from corresponding base/top edges
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

        # Show base, then animate the height to extrude  it out into a prism 
        self.add(base, top, sides)
        self.play(h.animate.set_value(2), run_time=2, rate_func=smooth)
        self.wait()





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
