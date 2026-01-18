from manim import *
import numpy as np

from utils import (
    make_square, make_rect, make_equilateral_triangle, make_wall_from_edge, make_lid_triangle_from_wall,
    Face,
    se3_about_edge,
    poly_from_V, line3d,
    vertex_labels, hinge_lines,
    make_point_hinge_updater, make_face_hinge_updater, make_child_hinge_updater, make_line_child_hinge_updater,
)


from utils.cylinder_net import map_flat_to_cylinder

# ----------------------------
# Demonstrate point about hinge
# ----------------------------
class HingePointDemo(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        p0 = np.array([-1.0, -1.0, -1.0])
        p1 = np.array([ 1.0,  1.0,  1.0])
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        x0 = np.array([0.5, 0.6, 0.0])
        dot = Dot3D(x0, color=RED)
        self.add(dot)

        theta = ValueTracker(0.0)
        dot.add_updater(make_point_hinge_updater(x0, p0, p1, theta, sign=+1))

        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)
        dot.clear_updaters()


# ----------------------------
# One face folding demo
# ----------------------------
class HingeFaceStepB(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        V, E = make_square(1.5, origin=(0.0, 0.0))
        V0 = V.copy()
        i0, i1 = E[0]
        p0, p1 = V0[int(i0)], V0[int(i1)]

        face = poly_from_V(V0, color=BLUE, fill_opacity=0.25, stroke_width=2)
        self.add(face, line3d(p0, p1, color=YELLOW))

        theta = ValueTracker(0.0)
        face.add_updater(make_face_hinge_updater(V0, p0, p1, theta, sign=+1))

        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        face.clear_updaters()
        self.wait(0.5)


# ----------------------------
# Four walls fold from center (your FiveSideCubeFourFaces)
# ----------------------------
class FiveSideCubeFourFaces(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 1.0
        VC, EC = make_square(s, (0, 0))
        VE, _  = make_square(s, (s, 0))
        VN, _  = make_square(s, (0, s))
        VW, _  = make_square(s, (-s, 0))
        VS, _  = make_square(s, (0, -s))

        centre = Face("centre", VC, color=BLUE)
        east   = Face("east",   VE, color=GREEN)
        north  = Face("north",  VN, color=RED)
        west   = Face("west",   VW, color=PURPLE)
        south  = Face("south",  VS, color=ORANGE)

        self.add(centre.poly, east.poly, north.poly, west.poly, south.poly)
        self.add(vertex_labels(VC))
        self.add(hinge_lines(VC, [EC[0], EC[1], EC[2], EC[3]], color=YELLOW))

        # hinges on centre
        e0, e1 = EC[1]  # [1,2]
        n0, n1 = EC[2]  # [2,3]
        w0, w1 = EC[3]  # [3,0]
        s0, s1 = EC[0]  # [0,1]
        pE0, pE1 = VC[int(e0)], VC[int(e1)]
        pN0, pN1 = VC[int(n0)], VC[int(n1)]
        pW0, pW1 = VC[int(w0)], VC[int(w1)]
        pS0, pS1 = VC[int(s0)], VC[int(s1)]

        theta = ValueTracker(0.0)

        east.poly.add_updater(make_face_hinge_updater(east.V0,  pE0, pE1, theta, sign=+1))
        north.poly.add_updater(make_face_hinge_updater(north.V0, pN0, pN1, theta, sign=+1))
        west.poly.add_updater(make_face_hinge_updater(west.V0,  pW0, pW1, theta, sign=+1))
        south.poly.add_updater(make_face_hinge_updater(south.V0, pS0, pS1, theta, sign=+1))

        self.play(theta.animate.set_value(-PI/2), run_time=3)
        self.play(theta.animate.set_value(0.0), run_time=3)

        for f in (east, north, west, south):
            f.poly.clear_updaters()


# ----------------------------
# Full cube net with lid (your FullCube)
# ----------------------------
class FullCube(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 2.0
        VC, EC = make_square(s, (0, 0))
        VE, _  = make_square(s, (s, 0))
        VN, _  = make_square(s, (0, s))
        VW, _  = make_square(s, (-s, 0))
        VS, _  = make_square(s, (0, -s))
        VExtra, _ = make_square(s, (2*s, 0))

        centre = Face("centre", VC, color=BLUE)
        east   = Face("east",   VE, color=GREEN)
        north  = Face("north",  VN, color=RED)
        west   = Face("west",   VW, color=PURPLE)
        south  = Face("south",  VS, color=ORANGE)
        extra  = Face("extra",  VExtra, color=WHITE)

        self.add(centre.poly, east.poly, north.poly, west.poly, south.poly, extra.poly)

        # centre hinges
        e0, e1 = EC[1]  # [1,2]
        n0, n1 = EC[2]  # [2,3]
        w0, w1 = EC[3]  # [3,0]
        s0, s1 = EC[0]  # [0,1]
        pE0, pE1 = VC[int(e0)], VC[int(e1)]
        pN0, pN1 = VC[int(n0)], VC[int(n1)]
        pW0, pW1 = VC[int(w0)], VC[int(w1)]
        pS0, pS1 = VC[int(s0)], VC[int(s1)]

        self.add(vertex_labels(VC))
        self.add(hinge_lines(VC, [EC[0], EC[1], EC[2], EC[3]], color=YELLOW))

        # child hinge on extra (local indices in extra face)
        extra_hinge_idx = (0, 3)

        theta = ValueTracker(0.0)
        phi   = ValueTracker(0.0)

        # walls
        east.poly.add_updater(make_face_hinge_updater(east.V0,  pE0, pE1, theta, sign=-1))
        north.poly.add_updater(make_face_hinge_updater(north.V0, pN0, pN1, theta, sign=+1))
        west.poly.add_updater(make_face_hinge_updater(west.V0,  pW0, pW1, theta, sign=+1))
        south.poly.add_updater(make_face_hinge_updater(south.V0, pS0, pS1, theta, sign=-1))

        # lid (extra) = child of east: parent hinge is pE0->pE1
        extra.poly.add_updater(
            make_child_hinge_updater(
                extra.V0,
                parent_p0=pE0, parent_p1=pE1,
                parent_tracker=theta, parent_sign=-1,
                child_hinge_idx=extra_hinge_idx,
                child_tracker=phi, child_sign=-1,
            )
        )

        # optional: visualize lid hinge line itself (tracks the same transforms)
        hinge_line = Line3D(VExtra[extra_hinge_idx[0]], VExtra[extra_hinge_idx[1]], color=BLUE)
        self.add(hinge_line)
        hinge_line.add_updater(
            make_line_child_hinge_updater(
                extra.V0,
                parent_p0=pE0, parent_p1=pE1,
                parent_tracker=theta, parent_sign=-1,
                line_idx=extra_hinge_idx,
                child_tracker=phi, child_sign=-1,
            )
        )

        self.play(theta.animate.set_value(PI/2), run_time=3)
        self.play(phi.animate.set_value(PI/2), run_time=3)

        for f in (east, north, west, south, extra):
            f.poly.clear_updaters()
        hinge_line.clear_updaters()


# ----------------------------
# Cuboid net fold (your CuboidNetFold)
# ----------------------------
class CuboidNetFold(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        l = 3.0
        w = 2.0
        h = 1.5

        V_base, _   = make_rect(l, w, origin=(0, 0))
        V_east, _   = make_rect(h, w, origin=(l, 0))
        V_west, _   = make_rect(h, w, origin=(-h, 0))
        V_north, _  = make_rect(l, h, origin=(0, w))
        V_south, _  = make_rect(l, h, origin=(0, -h))
        V_top, _    = make_rect(l, w, origin=(l + h, 0))

        base  = Face("base",  V_base,  color=BLUE)
        east  = Face("east",  V_east,  color=GREEN)
        west  = Face("west",  V_west,  color=PURPLE)
        north = Face("north", V_north, color=RED)
        south = Face("south", V_south, color=ORANGE)
        top   = Face("top",   V_top,   color=YELLOW)

        self.add(base.poly, east.poly, west.poly, north.poly, south.poly, top.poly)

        # base verts: 0:(0,0), 1:(l,0), 2:(l,w), 3:(0,w)
        pE0, pE1 = V_base[1], V_base[2]
        pW0, pW1 = V_base[0], V_base[3]
        pN0, pN1 = V_base[3], V_base[2]
        pS0, pS1 = V_base[0], V_base[1]

        theta = ValueTracker(0.0)
        phi   = ValueTracker(0.0)

        east.poly.add_updater(make_face_hinge_updater(east.V0,  pE0, pE1, theta, sign=-1))
        west.poly.add_updater(make_face_hinge_updater(west.V0,  pW0, pW1, theta, sign=+1))
        north.poly.add_updater(make_face_hinge_updater(north.V0, pN0, pN1, theta, sign=+1))
        south.poly.add_updater(make_face_hinge_updater(south.V0, pS0, pS1, theta, sign=-1))

        # top is child of east, hinge is its left edge (0,3) in its own vertices
        top_hinge_idx = (0, 3)
        top.poly.add_updater(
            make_child_hinge_updater(
                top.V0,
                parent_p0=pE0, parent_p1=pE1,
                parent_tracker=theta, parent_sign=-1,
                child_hinge_idx=top_hinge_idx,
                child_tracker=phi, child_sign=-1,
            )
        )

        self.wait(1)
        self.play(theta.animate.set_value(PI/2), run_time=3, rate_func=smooth)
        self.play(phi.animate.set_value(PI/2), run_time=2.5, rate_func=smooth)
        self.wait(1)
        self.play(phi.animate.set_value(0.0), run_time=2, rate_func=smooth)
        self.play(theta.animate.set_value(0.0), run_time=3, rate_func=smooth)

        for f in (east, west, north, south, top):
            f.poly.clear_updaters()
        self.wait(0.5)



# ----------------------------
# Triangular Prism (refactored)
# ----------------------------
class TriangularPrism(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.add(ThreeDAxes())

        s = 1.0  # triangle side
        h = 1.0  # wall "height" (net offset outward)

        # Base triangle + three walls + lid attached to wall0 outer edge
        V_base, E_base = make_equilateral_triangle(s)
        V_wall0, _ = make_wall_from_edge(V_base, E_base[0], h)  # AB
        V_wall1, _ = make_wall_from_edge(V_base, E_base[1], h)  # BC
        V_wall2, _ = make_wall_from_edge(V_base, E_base[2], h)  # CA
        V_lid, _, lid_hinge_idx = make_lid_triangle_from_wall(V_wall0, s)

        # Faces
        base  = Face("base",  V_base,  color=BLUE)
        wall0 = Face("wall0", V_wall0, color=GREEN)
        wall1 = Face("wall1", V_wall1, color=RED)
        wall2 = Face("wall2", V_wall2, color=PURPLE)
        lid   = Face("lid",   V_lid,   color=YELLOW)

        self.add(base.poly, wall0.poly, wall1.poly, wall2.poly, lid.poly)

        # Hinge lines are base edges (world-fixed)
        (a0, a1) = E_base[0]
        (b0, b1) = E_base[1]
        (c0, c1) = E_base[2]
        pA0, pA1 = V_base[int(a0)], V_base[int(a1)]
        pB0, pB1 = V_base[int(b0)], V_base[int(b1)]
        pC0, pC1 = V_base[int(c0)], V_base[int(c1)]

        # Optional labels
        self.add(vertex_labels(V_base))

        theta = ValueTracker(0.0)  # walls fold
        phi   = ValueTracker(0.0)  # lid folds relative to wall0

        # Walls fold about their corresponding base edges
        wall0.poly.add_updater(make_face_hinge_updater(wall0.V0, pA0, pA1, theta, sign=+1))
        wall1.poly.add_updater(make_face_hinge_updater(wall1.V0, pB0, pB1, theta, sign=+1))
        wall2.poly.add_updater(make_face_hinge_updater(wall2.V0, pC0, pC1, theta, sign=+1))

        # Lid is child of wall0:
        # parent hinge = base edge AB (pA0->pA1) with same theta
        # child hinge in lid local indices = lid_hinge_idx (0,1)
        lid.poly.add_updater(
            make_child_hinge_updater(
                lid.V0,
                parent_p0=pA0, parent_p1=pA1,
                parent_tracker=theta, parent_sign=+1,
                child_hinge_idx=lid_hinge_idx,
                child_tracker=phi, child_sign=+1,
            )
        )

        # Animate
        self.begin_ambient_camera_rotation(rate=0.5)
        self.play(theta.animate.set_value(-PI/2), run_time=3)
        self.play(phi.animate.set_value(-PI/2), run_time=3)

        # Cleanup
        for f in (wall0, wall1, wall2, lid):
            f.poly.clear_updaters()
        self.stop_ambient_camera_rotation()
        self.wait(0.5)


# ---------------------------------
# Cylinder with two lids (refactored)
# ---------------------------------
class CylinderWithTwoLids(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.add(ThreeDAxes())

        R = 1.0
        H = 2.0
        L = 2 * PI * R
        N = 40

        alpha = ValueTracker(0.0)  # roll side
        beta  = ValueTracker(0.0)  # fold lids

        panels = VGroup()
        dx = L / N

        # Flat sheet in xy-plane (z=0), centered in x
        for i in range(N):
            x0 = -L / 2 + i * dx
            x1 = -L / 2 + (i + 1) * dx

            v0 = np.array([x0, 0.0, 0.0])
            v1 = np.array([x1, 0.0, 0.0])
            v2 = np.array([x1, H, 0.0])
            v3 = np.array([x0, H, 0.0])

            poly = Polygon(
                v0, v1, v2, v3,
                color=BLUE,
                fill_color=BLUE,
                fill_opacity=0.7,
                stroke_width=0.5,
            )
            poly.base_vertices = np.array([v0, v1, v2, v3], dtype=np.float64)
            panels.add(poly)

        self.add(panels)

        # Axis reference
        axis_line = Line3D(start=[0, -0.5, 0], end=[0, H + 0.5, 0], color=YELLOW)
        self.add(axis_line)

        # Lids: represented by a VMobject polyline (smooth)
        disc_samples = 80
        angles = np.linspace(0, 2 * PI, disc_samples, endpoint=False)

        def build_lid(flat_center_y: float, cap_y: float):
            lid = VMobject(
                stroke_color=YELLOW,
                stroke_width=2,
                fill_color=YELLOW,
                fill_opacity=0.7,
            )

            base_pts = []
            cap_pts  = []
            for a in angles:
                # flat net circle
                x_flat = R * np.cos(a)
                y_flat = flat_center_y + R * np.sin(a)
                base_pts.append(np.array([x_flat, y_flat, 0.0], dtype=np.float64))

                # cap circle on cylinder (circle in plane y=cap_y)
                x_cap = R * np.cos(a)
                z_cap = R * np.sin(a)
                cap_pts.append(np.array([x_cap, cap_y, z_cap], dtype=np.float64))

            base_pts.append(base_pts[0])
            cap_pts.append(cap_pts[0])

            base_pts = np.array(base_pts, dtype=np.float64)
            cap_pts  = np.array(cap_pts, dtype=np.float64)

            lid.set_points_smoothly(base_pts)
            lid.base_points = base_pts
            lid.cap_points  = cap_pts
            return lid

        top_center_flat_y = H + R + 0.2
        bottom_center_flat_y = -R - 0.2

        top_lid = build_lid(flat_center_y=top_center_flat_y, cap_y=H)
        bottom_lid = build_lid(flat_center_y=bottom_center_flat_y, cap_y=0.0)

        self.add(top_lid, bottom_lid)

        # Updaters: panels roll to cylinder
        for poly in panels:
            base = poly.base_vertices.copy()

            def panel_updater(m, base_vertices=base):
                t = float(alpha.get_value())
                new_vertices = [map_flat_to_cylinder(v, t, R=R, L=L) for v in base_vertices]
                m.become(Polygon(
                    *new_vertices,
                    color=BLUE,
                    fill_color=BLUE,
                    fill_opacity=0.7,
                    stroke_width=0.5,
                ))
                return m

            poly.add_updater(panel_updater)

        # Updaters: lids lerp base->cap
        def lid_updater(m: VMobject):
            t = float(beta.get_value())
            pts = (1 - t) * m.base_points + t * m.cap_points
            m.set_points_smoothly(pts)
            return m

        top_lid.add_updater(lid_updater)
        bottom_lid.add_updater(lid_updater)

        # Animate
        self.wait(1)
        self.begin_ambient_camera_rotation(rate=0.2)
        self.play(alpha.animate.set_value(1.0), run_time=4, rate_func=smooth)
        self.play(beta.animate.set_value(1.0), run_time=3, rate_func=smooth)
        self.wait(2)
        self.stop_ambient_camera_rotation()

        # Cleanup
        for p in panels:
            p.clear_updaters()
        top_lid.clear_updaters()
        bottom_lid.clear_updaters()
        self.wait(0.5)


# ----------------------------
# L-shaped net fold (refactored)
# ----------------------------
class LShapedNetFold(ThreeDScene):
    """
    Net layout (side s, XY plane):

        [T-blue]
 [L-yellow][C1-green]
        [C2-red root]
        [C3-purple][R-white]

    - C2 is root (fixed)
    - C1 folds up from C2
    - T folds up from C1 (eta)
    - L hinges from C1 (phi)
    - C3 folds down from C2
    - R hinges from C3 (phi)
    """
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 2.0

        # Flat net squares
        VC2, _ = make_square(s, (0, 0))      # root
        VC1, _ = make_square(s, (0, s))      # above root
        VT,  _ = make_square(s, (0, 2*s))    # above C1
        VC3, _ = make_square(s, (0, -s))     # below root
        VL,  _ = make_square(s, (-s, s))     # left of C1
        VR,  _ = make_square(s, ( s, -s))    # right of C3

        C2 = Face("C2", VC2, color=RED)      # root
        C1 = Face("C1", VC1, color=GREEN)
        T  = Face("T",  VT,  color=BLUE)
        C3 = Face("C3", VC3, color=PURPLE)
        Lf = Face("L",  VL,  color=YELLOW)
        Rf = Face("R",  VR,  color=WHITE)

        self.add(T.poly, C1.poly, C2.poly, C3.poly, Lf.poly, Rf.poly)

        # Hinge lines on C2 (world-fixed)
        # For square ordering: 0:(0,0),1:(s,0),2:(s,s),3:(0,s)
        # C2–C1 hinge is top edge of C2: 3 -> 2
        pC2C1_0, pC2C1_1 = VC2[3], VC2[2]
        # C2–C3 hinge is bottom edge of C2: 0 -> 1
        pC2C3_0, pC2C3_1 = VC2[0], VC2[1]

        theta = ValueTracker(0.0)  # fold column (C1 up, C3 down)
        phi   = ValueTracker(0.0)  # fold side flaps (L and R)
        eta   = ValueTracker(0.0)  # fold top flap T about top edge of C1

        # C1: hinge about C2 top edge
        C1.poly.add_updater(make_face_hinge_updater(C1.V0, pC2C1_0, pC2C1_1, theta, sign=+1))

        # C3: hinge about C2 bottom edge, opposite direction
        C3.poly.add_updater(make_face_hinge_updater(C3.V0, pC2C3_0, pC2C3_1, theta, sign=-1))

        # T: child of C1
        # parent hinge = (pC2C1_0,pC2C1_1) with theta
        # child hinge = top edge of C1 in its own vertex indices: 3 -> 2
        T.poly.add_updater(
            make_child_hinge_updater(
                T.V0,
                parent_p0=pC2C1_0, parent_p1=pC2C1_1,
                parent_tracker=theta, parent_sign=+1,
                child_hinge_idx=(3, 2),
                child_tracker=eta, child_sign=+1,
            )
        )

        # L: child of C1, hinge about left edge of C1 (3 -> 0)
        Lf.poly.add_updater(
            make_child_hinge_updater(
                Lf.V0,
                parent_p0=pC2C1_0, parent_p1=pC2C1_1,
                parent_tracker=theta, parent_sign=+1,
                child_hinge_idx=(3, 0),
                child_tracker=phi, child_sign=+1,
            )
        )

        # R: child of C3, hinge about right edge of C3 (1 -> 2)
        Rf.poly.add_updater(
            make_child_hinge_updater(
                Rf.V0,
                parent_p0=pC2C3_0, parent_p1=pC2C3_1,
                parent_tracker=theta, parent_sign=-1,
                child_hinge_idx=(1, 2),
                child_tracker=phi, child_sign=+1,
            )
        )

        # Animate
        self.play(theta.animate.set_value(PI/2), run_time=5)
        self.play(phi.animate.set_value(-PI/2), run_time=5)
        self.play(eta.animate.set_value(PI/2), run_time=5)

        # Optional unfold
        self.play(phi.animate.set_value(0.0), run_time=3)
        self.play(theta.animate.set_value(0.0), run_time=3)
        self.play(eta.animate.set_value(0.0), run_time=3)

        # Cleanup
        for f in (C1, T, Lf, C3, Rf):
            f.poly.clear_updaters()

        self.wait(0.5)
