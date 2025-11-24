from manim import *
import numpy as np
from typing import Tuple, Iterable
from numpy.typing import NDArray

# -----------------------------------------------------------------------------------------------
# HELPER FUNCTIONS 
# -----------------------------------------------------------------------------------------------
def unit(v: Iterable[float]):
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
    k = unit(axis)
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
        u = unit(q1 - q0)
        k = unit(p1 - p0)

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
        u = unit(q1 - q0)
        k = unit(p1 - p0)

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

        s = 1.0

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
        self.play(theta.animate.set_value(-PI/2), run_time=3)

        # Step 2: now fold the lid (extra) relative to east
        self.play(phi.animate.set_value(-PI/2), run_time=3)

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

            # outward unit normal (right of edge)
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