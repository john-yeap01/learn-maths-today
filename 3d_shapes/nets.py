from manim import *
import numpy as np
from typing import Tuple, Iterable
from numpy.typing import NDArray

# -----------------------------------------------------------------------------------------------
# HELPER FUNCTIONS 
# -----------------------------------------------------------------------------------------------
def unit_vec(v: Iterable[float]):
    x = np.asarray(v, dtype=float).reshape(3)

    if x.size != 3:
        raise ValueError("Expected 3 components")
    
    n = float(np.linalg.norm(x))
    if n == 0.0:
        raise ValueError("zero length vector")
    
    return (x/n).astype(np.float64, copy=False)

# returns matrix cross form of a vector
def _skew(v: Iterable[float]):
    vec = np.asarray(v, dtype=float)
    if vec.size != 3:
        raise ValueError("Expected 3 components")
    vec = vec.reshape(3)
    x_, y_, z_ = [float(c) for c in vec]
    skew = np.array([[0, -z_, y_], 
                     [z_, 0, -x_],
                     [-y_, x_, 0]], dtype=np.float64)
    return skew

# returns 3 by 3 generalised rotation matrix given an axis and an angle around that axis
def rodrigues(axis: Iterable[float], theta: float):
    k = unit_vec(axis)
    K = _skew(k)

    s = float(np.sin(theta))
    c = float(np.cos(theta))

    I = np.eye(3, dtype=np.float64)

    R = I + s*K + (1-c)*(K@K)

    return R

def make_square(side: float = 1.0, origin: Tuple[float, float] = (0.0, 0.0)) -> Tuple[np.ndarray, np.ndarray]:
    ox, oy = [float(c) for c in origin]
    s = float(side)
    V = np.array([
        [ox+0, oy+0, 0.0],
        [ox+s, oy+0, 0.0],
        [ox+s, oy+s, 0.0],
        [ox+0, oy+s, 0.0],
    ], dtype=np.float64)
    # edges as index pairs into V (counter-clockwise)
    E = np.array([[0,1],[1,2],[2,3],[3,0]], dtype=np.int32)
    return V, E

def se3_about_edge(p0: NDArray[np.float64],
                   p1: NDArray[np.float64],
                   theta: float):
    p0 = np.asarray(p0, dtype=np.float64).reshape(3)
    p1 = np.asarray(p1, dtype=np.float64).reshape(3)
    R = rodrigues(p1 - p0, theta)   # axis normalization handled inside rodrigues
    t = p0 - (R @ p0)               # keep the entire line fixed
    return R, t

# -----------------------------------------------------------------------------------------------

# Demonstrate how a point moves about a general hinge
class HingePointDemo(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        axes = ThreeDAxes()
        self.add(axes)

        # hinge line (yellow): p0 -> p1 (x-axis segment)
        p0 = np.array([-1.0, -1.0, -1.0])
        p1 = np.array([1.0, 1.0, 1.0])
        hinge_line = Line3D(p0, p1, color=YELLOW)
        a = Dot3D(p0, color=YELLOW); b = Dot3D(p1, color=YELLOW)

        # point to rotate (not on the line)
        x0 = np.array([0.5, 0.6, 0.0])
        dot = Dot3D(x0, color=RED)

        theta = ValueTracker(0)

        def update_dot(mob: Mobject):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            mob.move_to(R @ x0 + t)
            return mob

        dot.add_updater(update_dot)
        self.add(hinge_line, a, b, dot)

        self.play(theta.animate.set_value(PI/2), run_time=2)  # 0→+90°
        self.play(theta.animate.set_value(0.0), run_time=2)   # back to 0
        dot.clear_updaters()

class HingeFaceStepA(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        axes = ThreeDAxes()
        self.add(axes)

        V, E = make_square(1.5, origin=(0.0, 0.0))  # 1.5×1.5 at (0,0)

        # 1) draw the face from V (CCW). Slight fill so edges are visible.
        face = Polygon(*[tuple(p) for p in V], color=BLUE, fill_opacity=0.2, stroke_width=2)
        self.add(face)

        # 2) highlight the hinge edge [0,1] in yellow (and mark endpoints)
        i0, i1 = [int(x) for x in E[0]]          # first edge
        p0, p1 = V[i0], V[i1]
        hinge = Line3D(p0, p1, color=YELLOW)
        a = Dot3D(p0, color=YELLOW); b = Dot3D(p1, color=YELLOW)
        self.add(hinge, a, b)

        # 3) label vertices briefly to confirm ordering/orientation
        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0,0,0.02]))
            for i, p in enumerate(V)
        ])
        self.add(labels)

        self.wait(1)

class HingeFaceStepB(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        V, E = make_square(1.5, origin=(0.0, 0.0))
        V0 = V.copy()                              # keep an immutable copy V0
        i0, i1 = [int(x) for x in E[0]]
        p0, p1 = V0[i0], V0[i1]                    # hinge endpoints (world-fixed)

        # draw original face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        theta = ValueTracker(0.0)

        # updater: recompute rotated vertices and write back into the polygon
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            # Manim expects corners in order (no need to repeat first point for Polygon)
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        face.add_updater(update_face)

        # animate: 0 → +90° → 0 (watch that the yellow edge stays exactly in place)
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        face.clear_updaters()
        self.wait(0.5)

class HingeFaceStepC(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        V, E = make_square(1.5, origin=(0.0, 0.0))
        V0 = V.copy()
        i0, i1 = [int(x) for x in E[0]]
        p0, p1 = V0[i0], V0[i1]                    # hinge endpoints (world-fixed)

        # child geometry (we'll only show its hinge, no face)
        V1, E1 = make_square(1.5, origin=(0.0, 0.0))
        V1 = V1.copy()
        j0, j1 = [int(x) for x in E1[2]]          # child hinge indices 3->2
        j0, j1 = j1, j0                            # now direction is 3 -> 2
        q0, q1 = V1[j0], V1[j1]

        # (kept for minimal change; not strictly needed)
        u = unit_vec(q1 - q0)
        k = unit_vec(p1 - p0)

        # child hinge + dots (initial placement uses V1 directly)
        child_hinge = Line3D(V1[j0], V1[j1], color=GREEN)
        aG = Dot3D(V1[j0], color=GREEN)
        bG = Dot3D(V1[j1], color=GREEN)
        self.add(child_hinge, aG, bG)

        # draw original parent face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        theta = ValueTracker(0.0)

        # --- Updaters ---

        # parent folds about its hinge
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m
        face.add_updater(update_face)

        # child hinge rides along with the *same* parent transform
        def update_child_hinge(m: Line3D):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            P = (R @ V1[[j0, j1]].T).T + t      # uses broadcasting since the V1 points are rows 
            new_line = Line3D(P[0], P[1], color=GREEN)
            m.become(new_line)
            return m

        child_hinge.add_updater(update_child_hinge)

        def update_child_dot(m: Dot3D, idx: int):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            m.move_to((R @ V1[idx]) + t)
            return m
        aG.add_updater(lambda m: update_child_dot(m, j0))
        bG.add_updater(lambda m: update_child_dot(m, j1))

        # animate: 0 → +90° → 0
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        # cleanup
        face.clear_updaters()
        child_hinge.clear_updaters()
        aG.clear_updaters()
        bG.clear_updaters()
        self.wait(0.5)

class HingeFaceStepD(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        V, E = make_square(1.5, origin=(0.0, 0.0))
        V0 = V.copy()
        i0, i1 = [int(x) for x in E[0]]
        p0, p1 = V0[i0], V0[i1]                # hinge endpoints (world-fixed)

        # child geometry (we'll only show its hinge, no face)
        V1, E1 = make_square(1.5, origin=(0.0, 0.0))
        V1 = V1.copy()
        j0, j1 = [int(x) for x in E1[2]]       # child hinge indices 3->2
        j0, j1 = j1, j0                        # now direction is 3 -> 2
        q0, q1 = V1[j0], V1[j1]

        # child face
        child_face = Polygon(*[tuple(p) for p in V1], color=GREEN, fill_opacity=0.25, stroke_width=2)

        # rotate the child face by -90 first (behind the scenes)
        R1, t1 = se3_about_edge(q0, q1, -PI/2)
        V_rot1 = (R1 @ V1.T).T + t1
        V1 = V_rot1
        child_face.set_points_as_corners([tuple(p) for p in V_rot1])
        self.add(child_face)

        # (kept for minimal change; not strictly needed)
        u = unit_vec(q1 - q0)
        k = unit_vec(p1 - p0)

        # child hinge + dots (initial placement uses V1 directly)
        child_hinge = Line3D(V1[j0], V1[j1], color=GREEN)
        aG = Dot3D(V1[j0], color=GREEN)
        bG = Dot3D(V1[j1], color=GREEN)
        self.add(child_hinge, aG, bG)

        # draw original parent face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        # --- Rotation Updaters ---
        theta = ValueTracker(0.0)

        # parent folds about its hinge
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m
        face.add_updater(update_face)

        # child hinge rides along with the *same* parent transform
        def update_child_hinge(m: Line3D):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            P = (R @ V1[[j0, j1]].T).T + t      # uses broadcasting since the V1 points are rows 
            new_line = Line3D(P[0], P[1], color=GREEN)
            m.become(new_line)
            return m
        
        def update_child_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V1.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        child_hinge.add_updater(update_child_hinge)
        child_face.add_updater(update_child_face)

        def update_child_dot(m: Dot3D, idx: int):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            m.move_to((R @ V1[idx]) + t)
            return m
        aG.add_updater(lambda m: update_child_dot(m, j0))
        bG.add_updater(lambda m: update_child_dot(m, j1))

        # animate: 0 → +90° → 0
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        # cleanup
        face.clear_updaters()
        child_hinge.clear_updaters()
        aG.clear_updaters()
        bG.clear_updaters()
        self.wait(0.5)

# ----------------------------------------------------------------------
# Face class before we go on 
# ---------------------------------------------------------------------

class Face:
    def __init__(self, name: str, V0: NDArray[np.float64], hinge_edge: Tuple[int, int] | None = None,
                 color=BLUE):
        """
        V0: (4,3) base vertices in world coords (usually flat in z=0)
        hinge_edge: indices into V0 (i0, i1) for the edge this face hinges about
                    (None if this face is static / root)
        """
        self.name = name
        self.V0 = V0.copy() 
        self.hinge_edge = hinge_edge
        self.color = color

        # create the polygon in its base position
        self.poly = Polygon(*[tuple(p) for p in self.V0],
                            color=self.color,
                            fill_opacity=0.25,
                            stroke_width=2)




# ----------------------------------------------------------------------
# Creating a 
# ---------------------------------------------------------------------

# Creates five faces on the XY plane 
class FiveSideCube(ThreeDScene):

    def construct(self):
    
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        number_of_faces = 5


        # Counterclockwise faces with side s
        s = 1.0
        VC, EC = make_square(s, (0,0))
        center_face = Face(1,1)

        VE, EE = make_square(s, (s,0))
        VN, EN = make_square(s, (0,s))
        VW, EW = make_square(s, (-s,0))
        VS, ES = make_square(s, (0,-s))


        centre_square = Polygon(*VC)
        east_square = Polygon(*VE)
        north_square = Polygon(*VN)
        west_square = Polygon(*VW)
        south_square = Polygon(*VS)

        self.add(centre_square, east_square, west_square, south_square, north_square)
        for p in VC:
            print(p)

        self.wait(1)


        

class FiveSideCubeOneFace(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 1.0

        # --- Build the flat net geometry (same as your existing code) ---
        VC, EC = make_square(s, (0, 0))    # centre
        VE, EE = make_square(s, (s, 0))    # east

        # centre vertices (0,0,0), (s,0,0), (s,s,0), (0,s,0)
        # edges: [0,1], [1,2], [2,3], [3,0]
        # shared edge with east is (s,0) -> (s,s) = vertices 1 -> 2 on centre

        centre_face = Face("centre", VC, hinge_edge=None, color=BLUE)
        east_face   = Face("east",   VE, hinge_edge=(0, 3), color=GREEN)
        # NOTE: on the east face, the same physical edge is between vertices 0 -> 3 (s,0,0) -> (s,s,0)

        self.add(centre_face.poly, east_face.poly)

        # --- Visualise the shared hinge edge on the centre face ---
        # centre hinge = vertices 1 -> 2
        i0, i1 = EC[1]       # edge [1,2]
        p0, p1 = VC[i0], VC[i1]

        hinge_line = Line3D(p0, p1, color=YELLOW)
        a = Dot3D(p0, color=YELLOW)
        b = Dot3D(p1, color=YELLOW)
        self.add(hinge_line, a, b)

        # Optional: label centre vertices to confirm indices
        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0, 0, 0.02]))
            for i, p in enumerate(VC)
        ])
        self.add(labels)

        # --- Folding animation for the east face only ---
        theta = ValueTracker(0.0)

        def update_east(m: Polygon):
            # rotate east_face.V0 about the *world* hinge line of the centre
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ east_face.V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        east_face.poly.add_updater(update_east)

        # Animate: flat net (0) -> vertical (±90°) -> back
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        east_face.poly.clear_updaters()
        self.wait(0.5)




class FiveSideCubeFourFaces(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 1.0

        # --- Build flat net geometry ---
        VC, EC = make_square(s, (0, 0))      # centre
        VE, EE = make_square(s, (s, 0))      # east
        VN, EN = make_square(s, (0, s))      # north
        VW, EW = make_square(s, (-s, 0))     # west
        VS, ES = make_square(s, (0, -s))     # south

        # Faces
        centre_face = Face("centre", VC, hinge_edge=None, color=BLUE)
        east_face   = Face("east",   VE, hinge_edge=(0, 3), color=GREEN)
        north_face  = Face("north",  VN, hinge_edge=(0, 1), color=RED)
        west_face   = Face("west",   VW, hinge_edge=(1, 2), color=PURPLE)
        south_face  = Face("south",  VS, hinge_edge=(2, 3), color=ORANGE)

        self.add(
            centre_face.poly,
            east_face.poly,
            north_face.poly,
            west_face.poly,
            south_face.poly
        )

        # --- Hinge lines on the centre face ---
        # centre vertices:
        # 0: (0,0,0)
        # 1: (s,0,0)
        # 2: (s,s,0)
        # 3: (0,s,0)

        # East hinge: centre edge [1,2]
        e0, e1 = EC[1]  # [1,2]
        pE0, pE1 = VC[e0], VC[e1]

        # North hinge: centre edge [2,3]
        n0, n1 = EC[2]  # [2,3]
        pN0, pN1 = VC[n0], VC[n1]

        # West hinge: centre edge [3,0]
        w0, w1 = EC[3]  # [3,0]
        pW0, pW1 = VC[w0], VC[w1]

        # South hinge: centre edge [0,1]
        s0, s1 = EC[0]  # [0,1]
        pS0, pS1 = VC[s0], VC[s1]

        # Optional: draw all hinge lines for debugging
        hinges = VGroup(
            Line3D(pE0, pE1, color=YELLOW),
            Line3D(pN0, pN1, color=YELLOW),
            Line3D(pW0, pW1, color=YELLOW),
            Line3D(pS0, pS1, color=YELLOW),
        )
        self.add(hinges)

        # Optional: label centre vertices to confirm indices
        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0, 0, 0.02]))
            for i, p in enumerate(VC)
        ])
        self.add(labels)

        # --- One parameter to fold all side faces at once ---
        theta = ValueTracker(0.0)

        # East face folds around centre edge [1,2]
        def update_east(m: Polygon):
            R, t = se3_about_edge(pE0, pE1, theta.get_value())
            V_rot = (R @ east_face.V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # North face folds around centre edge [2,3]
        def update_north(m: Polygon):
            R, t = se3_about_edge(pN0, pN1, theta.get_value())
            V_rot = (R @ north_face.V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # West face folds around centre edge [3,0]
        def update_west(m: Polygon):
            R, t = se3_about_edge(pW0, pW1, theta.get_value())
            V_rot = (R @ west_face.V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # South face folds around centre edge [0,1]
        def update_south(m: Polygon):
            R, t = se3_about_edge(pS0, pS1, theta.get_value())
            V_rot = (R @ south_face.V0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        east_face.poly.add_updater(update_east)
        north_face.poly.add_updater(update_north)
        west_face.poly.add_updater(update_west)
        south_face.poly.add_updater(update_south)

        # Animate: all four fold up, then back down
        self.play(theta.animate.set_value(-PI/2), run_time=3)
        self.play(theta.animate.set_value(0.0),   run_time=3)

        # Cleanup
        for f in (east_face, north_face, west_face, south_face):
            f.poly.clear_updaters()

        self.wait(0.5)


class FullCube(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 2.0

        # --- Build flat net geometry ---
        VC, EC = make_square(s, (0, 0))      # centre
        VE, EE = make_square(s, (s, 0))      # east
        VN, EN = make_square(s, (0, s))      # north
        VW, EW = make_square(s, (-s, 0))     # west
        VS, ES = make_square(s, (0, -s))     # south
        VExtra, _ = make_square(s, (2*s, 0)) # extra on the far side of east

        # Base copies (never mutated)
        VE0    = VE.copy()
        VN0    = VN.copy()
        VW0    = VW.copy()
        VS0    = VS.copy()
        VExtra0 = VExtra.copy()

        # Faces
        centre_face = Face("centre", VC, hinge_edge=None, color=BLUE)
        east_face   = Face("east",   VE, hinge_edge=(0, 3), color=GREEN)
        north_face  = Face("north",  VN, hinge_edge=(0, 1), color=RED)
        west_face   = Face("west",   VW, hinge_edge=(1, 2), color=PURPLE)
        south_face  = Face("south",  VS, hinge_edge=(2, 3), color=ORANGE)
        extra_face  = Face("extra",  VExtra, hinge_edge=(0, 3), color=WHITE)

        self.add(
            centre_face.poly,
            east_face.poly,
            north_face.poly,
            west_face.poly,
            south_face.poly,
            extra_face.poly
        )

        # --- Hinge lines on the centre face ---
        # centre vertices:
        # 0: (0,0,0)
        # 1: (s,0,0)
        # 2: (s,s,0)
        # 3: (0,s,0)

        # East hinge: centre edge [1,2]
        e0, e1 = EC[1]  # [1,2]
        pE0, pE1 = VC[e0], VC[e1]

        # North hinge: centre edge [2,3]
        n0, n1 = EC[2]  # [2,3]
        pN0, pN1 = VC[n0], VC[n1]

        # West hinge: centre edge [3,0]
        w0, w1 = EC[3]  # [3,0]
        pW0, pW1 = VC[w0], VC[w1]

        # South hinge: centre edge [0,1]
        s0, s1 = EC[0]  # [0,1]
        pS0, pS1 = VC[s0], VC[s1]

        # Shared edge E–Extra in flat net:
        # East:  VE[1] -> VE[2] = (2s,0,0) -> (2s,s,0)
        # Extra: VExtra[0] -> VExtra[3]   = (2s,0,0) -> (2s,s,0)
        extra_hinge_idx = (0, 3)    # edge in the frame of the extra

        # Optional hinge visuals
        hinges = VGroup(
            Line3D(pE0, pE1, color=YELLOW),
            Line3D(pN0, pN1, color=YELLOW),
            Line3D(pW0, pW1, color=YELLOW),
            Line3D(pS0, pS1, color=YELLOW),
        )
        self.add(hinges)

        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0, 0, 0.02]))
            for i, p in enumerate(VC)
        ])
        self.add(labels)

        # --- Two angles: walls vs lid ---
        theta = ValueTracker(0.0)  # walls fold from centre
        phi   = ValueTracker(0.0)  # lid folds from east

        # Helper: parent transform for east wall (centre edge [1,2])
        def parent_transform(angle: float):
            return se3_about_edge(pE0, pE1, angle)

        # East face folds around centre edge [1,2]
        def update_east(m: Polygon):
            R, t = parent_transform(theta.get_value())
            V_rot = (R @ VE0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # North face folds around centre edge [2,3]
        def update_north(m: Polygon):
            R, t = se3_about_edge(pN0, pN1, theta.get_value())
            V_rot = (R @ VN0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # West face folds around centre edge [3,0]
        def update_west(m: Polygon):
            R, t = se3_about_edge(pW0, pW1, theta.get_value())
            V_rot = (R @ VW0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # South face folds around centre edge [0,1]
        def update_south(m: Polygon):
            R, t = se3_about_edge(pS0, pS1, theta.get_value())
            V_rot = (R @ VS0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # Extra face = child of east
        def update_extra(m: Polygon):
            # 1) Apply the same transform as east (parent)
            R_p, t_p = parent_transform(theta.get_value())
            Vp = (R_p @ VExtra0.T).T + t_p

            # 2) Get the hinge edge E–Extra in world coords after parent transform
            h0, h1 = extra_hinge_idx
            H0, H1 = Vp[h0], Vp[h1]         # H0 - H1 defines the moved hinge line

            # 3) Rotate around that hinge by phi
            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            V_final = (R_c @ Vp.T).T + t_c

            m.set_points_as_corners([tuple(p) for p in V_final])          # set the new vertices of the polygon
            return m

        east_face.poly.add_updater(update_east)
        north_face.poly.add_updater(update_north)
        west_face.poly.add_updater(update_west)
        south_face.poly.add_updater(update_south)
        extra_face.poly.add_updater(update_extra)

        # Also visualise the child hinge line from the extra face reference frame
        extra_hinge_line = Line3D(VExtra0[extra_hinge_idx[0]], VExtra0[extra_hinge_idx[1]], color=BLUE)
        self.add(extra_hinge_line)

        def update_extra_hinge_line(m: Line3D):

            # Rotate about the centre-east face edge
            R_p, t_p = parent_transform(theta.get_value())
            Vp = (R_p @ VExtra0.T).T + t_p
            h0, h1 = extra_hinge_idx
            H0, H1 = Vp[h0], Vp[h1]

            # same child rotation as extra face -- technically dont need this
            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            H0_f = R_c @ H0 + t_c
            H1_f = R_c @ H1 + t_c

            m.become(Line3D(H0_f, H1_f, color=BLUE))            # set the new points of the line
            return m

        extra_hinge_line.add_updater(update_extra_hinge_line)

        # --- Animate: walls up, then lid down ---
        # Step 1: fold all walls up (extra follows east, still coplanar with it)
        self.play(theta.animate.set_value(-PI/2), run_time=10)

        # Step 2: now fold the lid (extra) relative to east
        self.play(phi.animate.set_value(-PI/2), run_time=10)

        # Cleanup
        for f in (east_face, north_face, west_face, south_face, extra_face):
            f.poly.clear_updaters()
        extra_hinge_line.clear_updaters()

        self.wait(0.5)



class TriangularPrism(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.add(ThreeDAxes())

        s = 1.0       # triangle side length
        h = 1.0       # prism height (length of rectangular walls)

        # --- Helper: build base equilateral triangle in xy-plane ---
        # Vertices: A(0,0,0), B(s,0,0), C(s/2, sqrt(3)/2 * s, 0)
        def make_equilateral_triangle(side: float):
            A = np.array([0.0, 0.0, 0.0])
            B = np.array([side, 0.0, 0.0])
            C = np.array([0.5 * side, np.sqrt(3) / 2 * side, 0.0])
            V = np.array([A, B, C])
            E = [(0, 1), (1, 2), (2, 0)]  # AB, BC, CA
            return V, E

        # --- Helper: make rectangle (wall) from one edge of triangle ---
        # For CCW triangle, interior is on the left of each directed edge.
        # Outward normal is "right" side (dy, -dx) / |edge|.
        def make_wall_from_edge(V_tri: np.ndarray, edge, height: float):
            i, j = edge
            P0 = V_tri[i][:2]
            P1 = V_tri[j][:2]
            d = P1 - P0
            dx, dy = d
            length = np.hypot(dx, dy)
            if length == 0:
                raise ValueError("Degenerate edge for wall")

            # outward unit_vec normal (right of edge)
            u_out = np.array([dy, -dx]) / length
            Q0 = P0 + u_out * height
            Q1 = P1 + u_out * height

            # Embed in 3D, z = 0
            P0_3 = np.array([P0[0], P0[1], 0.0])
            P1_3 = np.array([P1[0], P1[1], 0.0])
            Q1_3 = np.array([Q1[0], Q1[1], 0.0])
            Q0_3 = np.array([Q0[0], Q0[1], 0.0])

            # Order: edge (P0 -> P1), then outer edge (Q1 -> Q0)
            V_wall = np.array([P0_3, P1_3, Q1_3, Q0_3])
            E_wall = [(0, 1), (1, 2), (2, 3), (3, 0)]
            return V_wall, E_wall

        # --- Helper: make lid triangle attached to the outer edge of wall 0 ---
        # Wall0 vertices: [P0, P1, Q1, Q0], outer edge is (2,3) = Q1->Q0.
        # We build an equilateral triangle using base Q0-Q1 lying in xy-plane.
        def make_lid_triangle_from_wall(V_wall0: np.ndarray, side: float):
            Q1 = V_wall0[2][:2]
            Q0 = V_wall0[3][:2]
            L0 = Q0
            L1 = Q1
            d = L1 - L0
            dx, dy = d
            length = np.hypot(dx, dy)
            if abs(length - side) > 1e-6:
                side_eff = length
            else:
                side_eff = side

            # Interior of lid on the "outside" of the wall.
            # Take interior normal to be "left" of edge L0->L1
            u_int = -np.array([-dy, dx]) / length       # <--- THIS CHOOSES THE SIDE

            h_tri = np.sqrt(3) / 2 * side_eff
            L2 = (L0 + L1) / 2 + u_int * h_tri         # apex

            L0_3 = np.array([L0[0], L0[1], 0.0])
            L1_3 = np.array([L1[0], L1[1], 0.0])
            L2_3 = np.array([L2[0], L2[1], 0.0])

            V_lid = np.array([L0_3, L1_3, L2_3])
            E_lid = [(0, 1), (1, 2), (2, 0)]
            lid_hinge = (0, 1)
            return V_lid, E_lid, lid_hinge


        # --- Build flat net geometry ---
        V_base, E_base = make_equilateral_triangle(s)

        # 3 walls, each from one edge of the base triangle
        V_wall0, E_wall0 = make_wall_from_edge(V_base, E_base[0], h)  # edge AB
        V_wall1, E_wall1 = make_wall_from_edge(V_base, E_base[1], h)  # edge BC
        V_wall2, E_wall2 = make_wall_from_edge(V_base, E_base[2], h)  # edge CA

        # Lid: triangle attached to the outer edge of wall0
        V_lid, E_lid, lid_hinge_idx = make_lid_triangle_from_wall(V_wall0, s)

        # Base copies (never mutated)
        V_base0   = V_base.copy()
        V_wall0_0 = V_wall0.copy()
        V_wall1_0 = V_wall1.copy()
        V_wall2_0 = V_wall2.copy()
        V_lid0    = V_lid.copy()

        # Faces
        base_face  = Face("base",  V_base,   hinge_edge=None,      color=BLUE)
        wall0_face = Face("wall0", V_wall0,  hinge_edge=(0, 1),    color=GREEN)
        wall1_face = Face("wall1", V_wall1,  hinge_edge=(0, 1),    color=RED)
        wall2_face = Face("wall2", V_wall2,  hinge_edge=(0, 1),    color=PURPLE)
        lid_face   = Face("lid",   V_lid,    hinge_edge=lid_hinge_idx, color=YELLOW)

        self.add(
            base_face.poly,
            wall0_face.poly,
            wall1_face.poly,
            wall2_face.poly,
            lid_face.poly
        )

        # --- Hinge lines on the base triangle ---
        # Base edges: E_base[0]=(0,1), E_base[1]=(1,2), E_base[2]=(2,0)
        (a0, a1) = E_base[0]
        (b0, b1) = E_base[1]
        (c0, c1) = E_base[2]

        pA0, pA1 = V_base0[a0], V_base0[a1]
        pB0, pB1 = V_base0[b0], V_base0[b1]
        pC0, pC1 = V_base0[c0], V_base0[c1]

        hinges = VGroup(
            Line3D(pA0, pA1, color=WHITE),
            Line3D(pB0, pB1, color=WHITE),
            Line3D(pC0, pC1, color=WHITE),
        )
        self.add(hinges)

        # Optional labels on base vertices
        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0, 0, 0.02]))
            for i, p in enumerate(V_base0)
        ])
        self.add(labels)

        # --- Two angles: walls vs lid ---
        theta = ValueTracker(0.0)   # 3 walls fold from base
        phi   = ValueTracker(0.0)   # lid folds from wall0

        # Parent transform for wall0: around base edge E_base[0] (A-B)
        def parent_transform(angle: float):
            return se3_about_edge(pA0, pA1, angle)

        # --- Updaters for walls ---
        def update_wall0(m: Polygon):
            R, t = parent_transform(theta.get_value())
            V_rot = (R @ V_wall0_0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        def update_wall1(m: Polygon):
            R, t = se3_about_edge(pB0, pB1, theta.get_value())
            V_rot = (R @ V_wall1_0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        def update_wall2(m: Polygon):
            R, t = se3_about_edge(pC0, pC1, theta.get_value())
            V_rot = (R @ V_wall2_0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # --- Lid = child of wall0 ---
        def update_lid(m: Polygon):
            # 1) Apply same transform as wall0 (parent)
            R_p, t_p = parent_transform(theta.get_value())
            Vp = (R_p @ V_lid0.T).T + t_p

            # 2) Hinge edge in lid frame after parent transform
            h0, h1 = lid_hinge_idx  # should be (0,1)
            H0, H1 = Vp[h0], Vp[h1]

            # 3) Rotate around that hinge by phi
            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            V_final = (R_c @ Vp.T).T + t_c

            m.set_points_as_corners([tuple(p) for p in V_final])
            return m

        wall0_face.poly.add_updater(update_wall0)
        wall1_face.poly.add_updater(update_wall1)
        wall2_face.poly.add_updater(update_wall2)
        lid_face.poly.add_updater(update_lid)

        # Visualise the lid hinge line
        lid_hinge_line = Line3D(V_lid0[lid_hinge_idx[0]], V_lid0[lid_hinge_idx[1]], color=YELLOW)
        self.add(lid_hinge_line)

        def update_lid_hinge_line(m: Line3D):
            # Same transforms as lid
            R_p, t_p = parent_transform(theta.get_value())
            Vp = (R_p @ V_lid0.T).T + t_p
            h0, h1 = lid_hinge_idx
            H0, H1 = Vp[h0], Vp[h1]

            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            H0_f = R_c @ H0 + t_c
            H1_f = R_c @ H1 + t_c

            m.become(Line3D(H0_f, H1_f, color=YELLOW))
            return m

        lid_hinge_line.add_updater(update_lid_hinge_line)

        # --- Animate: walls up, then lid down ---
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.5)

        # Step 1: fold all walls up from the base
        self.play(theta.animate.set_value(-PI / 2), run_time=3)

        # Step 2: fold the lid relative to wall0
        self.play(phi.animate.set_value(-PI / 2), run_time=3)

        # Cleanup
        for f in (wall0_face, wall1_face, wall2_face, lid_face):
            f.poly.clear_updaters()
        lid_hinge_line.clear_updaters()

        self.wait(0.5)

        self.stop_ambient_camera_rotation()





from manim import *
import numpy as np

class CylinderWithTwoLids(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.add(ThreeDAxes())

        R = 1.0          # cylinder radius
        H = 2.0          # cylinder height
        L = 2 * PI * R   # rectangle width = circumference
        N = 40           # number of vertical panels

        alpha = ValueTracker(0.0)  # roll side
        beta  = ValueTracker(0.0)  # fold lids

        panels = VGroup()
        dx = L / N

        # --- Flat sheet in xy-plane, centred in x ---
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
            poly.base_vertices = np.array([v0, v1, v2, v3])
            panels.add(poly)

        self.add(panels)

        # Cylinder axis (for reference)
        axis_line = Line3D(
            start=[0, -0.5, 0],
            end=[0, H + 0.5, 0],
            color=YELLOW,
        )
        self.add(axis_line)

        # --- Lids: start flat in xy-plane, above and below the sheet ---
        disc_samples = 80
        angles = np.linspace(0, 2 * PI, disc_samples, endpoint=False)

        # TOP lid (flat position)
        top_center_flat = np.array([0.0, H + R + 0.2, 0.0])
        top_lid = VMobject(
            stroke_color=YELLOW,
            stroke_width=2,
            fill_color=YELLOW,
            fill_opacity=0.7,
        )

        top_base_pts = []
        top_cap_pts  = []

        for a in angles:
            # flat net pos
            x_flat = top_center_flat[0] + R * np.cos(a)
            y_flat = top_center_flat[1] + R * np.sin(a)
            top_base_pts.append(np.array([x_flat, y_flat, 0.0]))

            # cap pos: circle in plane y = H
            x_cap = R * np.cos(a)
            y_cap = H
            z_cap = R * np.sin(a)
            top_cap_pts.append(np.array([x_cap, y_cap, z_cap]))

        top_base_pts.append(top_base_pts[0])
        top_cap_pts.append(top_cap_pts[0])

        top_base_pts = np.array(top_base_pts)
        top_cap_pts  = np.array(top_cap_pts)

        top_lid.set_points_smoothly(top_base_pts)
        top_lid.base_points = top_base_pts
        top_lid.cap_points  = top_cap_pts

        self.add(top_lid)

        # BOTTOM lid (flat position)
        bottom_center_flat = np.array([0.0, -R - 0.2, 0.0])
        bottom_lid = VMobject(
            stroke_color=YELLOW,
            stroke_width=2,
            fill_color=YELLOW,
            fill_opacity=0.7,
        )

        bottom_base_pts = []
        bottom_cap_pts  = []

        for a in angles:
            # flat net pos
            x_flat = bottom_center_flat[0] + R * np.cos(a)
            y_flat = bottom_center_flat[1] + R * np.sin(a)
            bottom_base_pts.append(np.array([x_flat, y_flat, 0.0]))

            # cap pos: circle in plane y = 0
            x_cap = R * np.cos(a)
            y_cap = 0.0
            z_cap = R * np.sin(a)
            bottom_cap_pts.append(np.array([x_cap, y_cap, z_cap]))

        bottom_base_pts.append(bottom_base_pts[0])
        bottom_cap_pts.append(bottom_cap_pts[0])

        bottom_base_pts = np.array(bottom_base_pts)
        bottom_cap_pts  = np.array(bottom_cap_pts)

        bottom_lid.set_points_smoothly(bottom_base_pts)
        bottom_lid.base_points = bottom_base_pts
        bottom_lid.cap_points  = bottom_cap_pts

        self.add(bottom_lid)

        # --- Mapping: flat sheet -> symmetric rolled cylinder ---
        def flat_to_cylinder_coord(p_flat):
            x, y, z = p_flat
            u = x / (L / 2)                # u in [-1, 1]
            theta_final = PI * (1.0 + u)   # u=-1->0, u=0->π, u=1->2π
            Xc = R * np.cos(theta_final)
            Zc = R * np.sin(theta_final)
            return np.array([Xc, y, Zc])

        def map_flat_to_cylinder(p_flat, t):
            p_cyl = flat_to_cylinder_coord(p_flat)
            return (1 - t) * p_flat + t * p_cyl

        # --- Updaters: side panels ---
        for poly in panels:
            base = poly.base_vertices.copy()

            def updater(m, base_vertices=base):
                t = alpha.get_value()
                new_vertices = [map_flat_to_cylinder(v, t) for v in base_vertices]
                m.become(Polygon(
                    *new_vertices,
                    color=BLUE,
                    fill_color=BLUE,
                    fill_opacity=0.7,
                    stroke_width=0.5,
                ))
                return m

            poly.add_updater(updater)

        # --- Updaters: lids (both use the same beta) ---
        def lid_updater_top(m: VMobject):
            t = beta.get_value()
            pts = (1 - t) * m.base_points + t * m.cap_points
            m.set_points_smoothly(pts)
            return m

        def lid_updater_bottom(m: VMobject):
            t = beta.get_value()
            pts = (1 - t) * m.base_points + t * m.cap_points
            m.set_points_smoothly(pts)
            return m

        top_lid.add_updater(lid_updater_top)
        bottom_lid.add_updater(lid_updater_bottom)

        # --- Animation ---
        self.wait(1)  # flat net with two lids

        self.begin_ambient_camera_rotation(rate=0.2)

        # 1) roll side into cylinder
        self.play(alpha.animate.set_value(1.0), run_time=4, rate_func=smooth)

        # 2) fold both lids into top and bottom caps
        self.play(beta.animate.set_value(1.0), run_time=3, rate_func=smooth)

        self.wait(2)
        self.stop_ambient_camera_rotation()

        for p in panels:
            p.clear_updaters()
        top_lid.clear_updaters()
        bottom_lid.clear_updaters()
        self.wait(0.5)


class LShapedNetFold(ThreeDScene):
    """
    Net layout (side length s, in XY plane):

        [T-blue]
 [L-yellow][C1-green]
        [C2-red root]
        [C3-purple][R-white]

    - C2 is the root (fixed)
    - C1 folds up from C2
    - T folds up from C1
    - L hinges from C1 (NOT from T)
    - C3 folds down from C2
    - R hinges from C3
    """
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        s = 2.0

        # --- Squares in the flat net ---
        VC2, _ = make_square(s, (0, 0))      # root
        VC1, _ = make_square(s, (0, s))      # above root
        VT,  _ = make_square(s, (0, 2*s))    # above C1
        VC3, _ = make_square(s, (0, -s))     # below root

        # Yellow flap now attaches to the LEFT of C1 (not left of T)
        VL,  _ = make_square(s, (-s, s))

        # Right flap attaches to the RIGHT of C3
        VR,  _ = make_square(s, ( s, -s))

        # Immutable copies
        VC2_0 = VC2.copy()
        VC1_0 = VC1.copy()
        VT0   = VT.copy()
        VC3_0 = VC3.copy()
        VL0   = VL.copy()
        VR0   = VR.copy()

        # Polys
        poly_T  = Polygon(*[tuple(p) for p in VT0],   color=BLUE,   fill_opacity=0.25)
        poly_C1 = Polygon(*[tuple(p) for p in VC1_0], color=GREEN,  fill_opacity=0.25)
        poly_C2 = Polygon(*[tuple(p) for p in VC2_0], color=RED,    fill_opacity=0.25)  # root
        poly_C3 = Polygon(*[tuple(p) for p in VC3_0], color=PURPLE, fill_opacity=0.25)
        poly_L  = Polygon(*[tuple(p) for p in VL0],   color=YELLOW, fill_opacity=0.25)
        poly_R  = Polygon(*[tuple(p) for p in VR0],   color=WHITE,  fill_opacity=0.25)

        self.add(poly_T, poly_C1, poly_C2, poly_C3, poly_L, poly_R)

        # --- Hinge lines (flat reference) ---
        # C2–C1 hinge: top edge of C2 = 3 -> 2
        pC2C1_0, pC2C1_1 = VC2_0[3], VC2_0[2]

        # C2–C3 hinge: bottom edge of C2 = 0 -> 1
        pC2C3_0, pC2C3_1 = VC2_0[0], VC2_0[1]

        theta = ValueTracker(0.0)  # fold column
        phi   = ValueTracker(0.0)  # fold flaps
        eta = ValueTracker(0.0)

        # --- Updaters ---

        # C1 folds about C2–C1
        def update_C1(m: Polygon):
            R, t = se3_about_edge(pC2C1_0, pC2C1_1, theta.get_value())
            V_rot = (R @ VC1_0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # T is child of C1: follow theta, then hinge about C1–T top edge by phi
        def update_T(m: Polygon):
            R_p, t_p = se3_about_edge(pC2C1_0, pC2C1_1, theta.get_value())
            VC1_p = (R_p @ VC1_0.T).T + t_p
            VT_p  = (R_p @ VT0.T).T   + t_p

            # hinge line is top edge of C1: 3 -> 2 (in world after parent)
            H0, H1 = VC1_p[3], VC1_p[2]
            R_c, t_c = se3_about_edge(H0, H1, eta.get_value())
            V_final = (R_c @ VT_p.T).T + t_c

            m.set_points_as_corners([tuple(p) for p in V_final])
            return m

        # L is child of C1: follow theta, then hinge about LEFT edge of C1 by phi
        def update_L(m: Polygon):
            R_p, t_p = se3_about_edge(pC2C1_0, pC2C1_1, theta.get_value())
            VC1_p = (R_p @ VC1_0.T).T + t_p
            VL_p  = (R_p @ VL0.T).T   + t_p

            # left edge of C1: 3 -> 0
            H0, H1 = VC1_p[3], VC1_p[0]
            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            V_final = (R_c @ VL_p.T).T + t_c

            m.set_points_as_corners([tuple(p) for p in V_final])
            return m

        # C3 folds about C2–C3 (opposite direction)
        def update_C3(m: Polygon):
            R, t = se3_about_edge(pC2C3_0, pC2C3_1, -theta.get_value())
            V_rot = (R @ VC3_0.T).T + t
            m.set_points_as_corners([tuple(p) for p in V_rot])
            return m

        # R is child of C3: follow -theta, then hinge about RIGHT edge of C3 by phi
        def update_R(m: Polygon):
            R_p, t_p = se3_about_edge(pC2C3_0, pC2C3_1, -theta.get_value())
            VC3_p = (R_p @ VC3_0.T).T + t_p
            VR_p  = (R_p @ VR0.T).T   + t_p

            # right edge of C3: 1 -> 2
            H0, H1 = VC3_p[1], VC3_p[2]
            R_c, t_c = se3_about_edge(H0, H1, phi.get_value())
            V_final = (R_c @ VR_p.T).T + t_c

            m.set_points_as_corners([tuple(p) for p in V_final])
            return m

        poly_C1.add_updater(update_C1)
        poly_T.add_updater(update_T)
        poly_L.add_updater(update_L)
        poly_C3.add_updater(update_C3)
        poly_R.add_updater(update_R)

        # --- Animate ---
        self.play(theta.animate.set_value(PI/2), run_time=5)
        self.play(phi.animate.set_value(-PI/2),   run_time=5)
        self.play(eta.animate.set_value(PI/2), run_time=5)

        # optional unfold
        self.play(phi.animate.set_value(0.0), run_time=3)
        self.play(theta.animate.set_value(0.0), run_time=3)
        self.play(eta.animate.set_value(0.0), run_time=3)


        for p in (poly_C1, poly_T, poly_L, poly_C3, poly_R):
            p.clear_updaters()

        self.wait(0.5)



# ----------------------------
# Rect helper (like make_square)
# ----------------------------
def make_rect(width: float, height: float, origin=(0.0, 0.0)):
    ox, oy = origin
    w = float(width)
    h = float(height)
    V = np.array([
        [ox + 0, oy + 0, 0.0],
        [ox + w, oy + 0, 0.0],
        [ox + w, oy + h, 0.0],
        [ox + 0, oy + h, 0.0],
    ], dtype=np.float64)
    E = np.array([[0,1],[1,2],[2,3],[3,0]], dtype=np.int32)
    return V, E

# ----------------------------
# Your hinge math (assumed already in file):
# - unit_vec
# - _skew
# - rodrigues
# - se3_about_edge
# ----------------------------

class CuboidNetFold(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # Cuboid dimensions
        l = 3.0   # length (x direction on base)
        w = 2.0   # width  (y direction on base)
        h = 1.5   # height (wall height)

        # ----------------------------
        # Build flat net in XY plane
        # Layout:
        #   north (l×h)
        # west (h×w)  base (l×w)  east (h×w)  top (l×w)
        #   south (l×h)
        # ----------------------------

        V_base, E_base   = make_rect(l, w, origin=(0, 0))
        V_east, _        = make_rect(h, w, origin=(l, 0))
        V_west, _        = make_rect(h, w, origin=(-h, 0))
        V_north, _       = make_rect(l, h, origin=(0, w))
        V_south, _       = make_rect(l, h, origin=(0, -h))
        V_top, _         = make_rect(l, w, origin=(l + h, 0))  # attached to outer edge of east

        # Immutable copies
        V_base0  = V_base.copy()
        V_east0  = V_east.copy()
        V_west0  = V_west.copy()
        V_north0 = V_north.copy()
        V_south0 = V_south.copy()
        V_top0   = V_top.copy()

        # Polys
        base  = Polygon(*map(tuple, V_base0),  color=BLUE,   fill_opacity=0.25)
        east  = Polygon(*map(tuple, V_east0),  color=GREEN,  fill_opacity=0.25)
        west  = Polygon(*map(tuple, V_west0),  color=PURPLE, fill_opacity=0.25)
        north = Polygon(*map(tuple, V_north0), color=RED,    fill_opacity=0.25)
        south = Polygon(*map(tuple, V_south0), color=ORANGE, fill_opacity=0.25)
        top   = Polygon(*map(tuple, V_top0),   color=YELLOW, fill_opacity=0.25)

        self.add(base, east, west, north, south, top)

        # ----------------------------
        # Hinge lines on BASE (world-fixed)
        # base verts:
        # 0:(0,0), 1:(l,0), 2:(l,w), 3:(0,w)
        # ----------------------------
        pE0, pE1 = V_base0[1], V_base0[2]  # right edge: 1->2 (east wall)
        pW0, pW1 = V_base0[0], V_base0[3]  # left edge: 0->3 (west wall)
        pN0, pN1 = V_base0[3], V_base0[2]  # top edge:  3->2 (north wall)
        pS0, pS1 = V_base0[0], V_base0[1]  # bottom:    0->1 (south wall)

        # trackers
        theta = ValueTracker(0.0)  # walls up/down
        phi   = ValueTracker(0.0)  # top lid close/open (child of east)

        # ----------------------------
        # Updaters for walls
        # Choose signs so they fold "up" (positive z) nicely.
        # If any face folds the wrong way, flip the sign on its angle.
        # ----------------------------
        def upd_east(m: Polygon):
            R, t = se3_about_edge(pE0, pE1, -theta.get_value())
            V = (R @ V_east0.T).T + t
            m.set_points_as_corners(list(map(tuple, V)))
            return m

        def upd_west(m: Polygon):
            R, t = se3_about_edge(pW0, pW1, theta.get_value())
            V = (R @ V_west0.T).T + t
            m.set_points_as_corners(list(map(tuple, V)))
            return m

        def upd_north(m: Polygon):
            R, t = se3_about_edge(pN0, pN1, theta.get_value())
            V = (R @ V_north0.T).T + t
            m.set_points_as_corners(list(map(tuple, V)))
            return m

        def upd_south(m: Polygon):
            R, t = se3_about_edge(pS0, pS1, -theta.get_value())
            V = (R @ V_south0.T).T + t
            m.set_points_as_corners(list(map(tuple, V)))
            return m

        # ----------------------------
        # Top face is a CHILD of east wall:
        # 1) move with east wall (parent)
        # 2) then rotate about the east wall's OUTER edge (x = l+h)
        # In top-face local vertices, the hinge is its LEFT edge: 0 -> 3
        # ----------------------------
        top_hinge_idx = (0, 3)

        def upd_top(m: Polygon):
            # parent = east wall transform
            R_p, t_p = se3_about_edge(pE0, pE1, -theta.get_value())
            Vp = (R_p @ V_top0.T).T + t_p

            # hinge line in world coords after parent
            h0, h1 = top_hinge_idx
            H0, H1 = Vp[h0], Vp[h1]

            # child rotation (closing lid)
            R_c, t_c = se3_about_edge(H0, H1, -phi.get_value())
            V_final = (R_c @ Vp.T).T + t_c

            m.set_points_as_corners(list(map(tuple, V_final)))
            return m

        east.add_updater(upd_east)
        west.add_updater(upd_west)
        north.add_updater(upd_north)
        south.add_updater(upd_south)
        top.add_updater(upd_top)

        # ----------------------------
        # Animate: net -> fold -> unfold
        # ----------------------------
        self.wait(1)

        # Fold walls up
        self.play(theta.animate.set_value(PI/2), run_time=3, rate_func=smooth)

        # Close the top
        self.play(phi.animate.set_value(PI/2), run_time=2.5, rate_func=smooth)

        self.wait(1)

        # Unfold (open top first, then walls down)
        self.play(phi.animate.set_value(0.0), run_time=2, rate_func=smooth)
        self.play(theta.animate.set_value(0.0), run_time=3, rate_func=smooth)

        # cleanup
        for m in (east, west, north, south, top):
            m.clear_updaters()
        self.wait(0.5)