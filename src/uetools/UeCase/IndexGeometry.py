from __future__ import annotations
from uetools import Case
import uetools
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Tuple, Union, Sequence, Optional, Iterable, Type
from numpy.typing import NDArray
from matplotlib.pyplot import Figure, Axes
import networkx as nx


Handle = str
NeighborVal = Union[str, List[str], None]
IndexRange = Tuple[int, int]
PatchIndexBounds = Tuple[IndexRange, IndexRange]


def _as_list(v: NeighborVal) -> List[str]:
    """Normalize neighbor value to a list of strings."""
    if v is None:
        return []
    if isinstance(v, list):
        return [x for x in v if x is not None]
    return [v]


@dataclass(frozen=True)
class Walls:
    W: str = "__WALL_W__"
    E: str = "__WALL_E__"
    S: str = "__WALL_S__"
    N: str = "__WALL_N__"


class TopoPatch:
    """
    Contains topological information about patches

    """

    def __init__(
        self,
        patch_info: Dict,
        location: Tuple[int,int],
        region_bounds: Tuple[Tuple[int], Tuple[int]],
        rm: NDArray,
        zm: NDArray
    ) -> None:
        self.info = patch_info
        self.location = location
        self.region_bounds = region_bounds
        self.rm = rm[
                        region_bounds[0][0]:region_bounds[0][1]+1, 
                        region_bounds[1][0]:region_bounds[1][1]+1
        ] 
        self.zm = zm[
                        region_bounds[0][0]:region_bounds[0][1]+1, 
                        region_bounds[1][0]:region_bounds[1][1]+1
        ] 
        return


    def plot_patch_geo(
        self,
        ax: Optional[Type[Axes]] = None,
        color: Optional[Type[Axes]] = 'k',
    ) -> Type[Figure]:
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()
        
        (nx, ny, _) = self.rm.shape
        for ix in range(nx):
            for iy in range(ny):
                ax.plot(
                        [self.rm[ix, iy, r] for r in [1, 2, 4, 3, 1]],
                        [self.zm[ix, iy, r] for r in [1, 2, 4, 3, 1]],
                        '-', color=color
                )

        return ax.get_figure()

    def plot_patch(self, ax: Type[Axes]) -> Type[Figure]:
        from random import random
        from matplotlib.patches import Rectangle
        x0, x1 = self.region_bounds[0]
        y0, y1 = self.region_bounds[1]

        width = x1 - x0 + 1
        height = y1 - y0 + 1 

        color = (random(), random(), random())

        rect = Rectangle(
            (x0-.5, y0-.5),
            width,
            height,
            facecolor=color,
            edgecolor=color,
            alpha=0.4,
            linewidth=1,
        )

        ax.add_patch(rect)
        return ax.get_figure()

class TopoGeo:
    """
    Contains topological geometry information about grid
    """
    def __init__(
        self, 
        case: Type[Case],
    ) -> None:
        """
        Creates topological grid information 
        """
        from yaml import safe_load_all
        import bisect

        self.patches = {}
        # Get index and physical dimensions from case
        self.rm = case.get('rm')
        self.zm = case.get('zm')
        self.geometry = case.get('geometry')[0].decode('UTF-8').strip()
        var_store = [   'ixlb', 
                        'ixrb', 
                        'ixpt1', 
                        'ixpt2', 
                        'iysptrx1',
                        'iysptrx2',
                        'nx',
                        'ny'
        ]
        [self.__setattr__(var, case.get(var)) for var in var_store]
        self.ixpt1 += 1
        self.ixpt2 += 1
        self.iyb = (0, self.ny)
        self.xpoints = len(self.ixlb)
        # Create x/y segmentation map
        self.ix = sorted((
            (val, name, i)
            for name, arr in {
                name: self.__getattribute__(name) for name in [
                        "ixlb",
                        "ixrb",
                        "ixpt1",
                        "ixpt2",
                ]
            }.items()
            for i, val in enumerate(arr)
        ))
        # TODO: remove this
        self.iysptrx1 += 1
        self.iysptrx2 += 1
        self.iy = sorted((
            (val, name, i)
            for name, arr in {
                name: self.__getattribute__(name) for name in [
                        "iyb",
                        "iysptrx1",
                ]
            }.items()
            for i, val in enumerate(arr)
        ))
        ix_indices = [int(v) for v, _, _ in self.ix]
        ix_cuts = bisect.bisect_left(ix_indices, self.ixrb[0]+1)
        self.region_xbounds = [ix_indices[:ix_cuts], ix_indices[ix_cuts:]]
        for i in range(self.xpoints):
            self.region_xbounds[i][-1] += 2
        self.region_ybounds = [int(v) for v, _, _ in self.iy]
        self.region_ybounds[-1] += 1
        


        # Read all geometry specifiers and create a dict from them
        yaml_file = "{}/{}".format( 
                        uetools.__path__[0], 
                        "yamls/geometries.yaml"
        )
        # Detect whether SF+/- and flip accordingly
        if "snowflake" in self.geometry.lower():
            sfmap = {105: 75, 135: 45, 165: 15}
            sfangle = int(self.geometry.replace('snowflake',''))
            if sfangle in [105, 135, 165]:
                self.geometry_info = geometry_decks[f"snowflake{sfmap[sfangle]}"]
                # Flip E/W neighbors for each patch to recreate SF- counterpart
                for name, patch in self.geometry_info['patches']:
                    for var in ['neighbors', 'boundaries']:
                        try:
                            patch[var] = swap_e_w(patch[var])
                        except:
                            pass
            


        geometry_decks = {}
        with open(yaml_file, 'r') as f:
            for data in safe_load_all(f):
                geometry_decks[data['name']] = data
        # Identify the correct topology deck and store it to self
        try:
            self.geometry_info = geometry_decks[self.geometry]
        except:
            self.geometry_info = None
            print(  f"Warning! Geometry {self.geometry} not specified "
                    "in specifier deck. \nPlotting and interpolation "
                    "routines will not work as expected.")
            return
        # Create patch objects containing all patch information
        if self.xpoints <= 2:
            for ixpt in range(1,1+self.xpoints):
                try:
                    patchlist = [patch for patch in  self.geometry_info[f'subregion{ixpt}']]
                except:
                    patchlist = list(self.geometry_info['patches'].keys())
                patchdict = {a: self.geometry_info['patches'][a] for a in patchlist}
                # Order segments according to indices
                order = TopologicalPacker2D(patchdict).arrange()
                bounds = self.patch_bounds_from_neighbors(
                    order,
                    patchdict,
                    self.region_xbounds[ixpt-1],
                    self.region_ybounds
                )
                # Create and assign patch objects, appending to self.patches
                self.patches.update(
                    {
                        name: TopoPatch(
                                self.geometry_info['patches'][name], 
                                location,
                                bounds[name],
                                self.rm,
                                self.zm
                        ) for name, location in order.items()
                    }
                )
        else:
            raise NotImplementedError("xpoints>2 not implemented")
        return


    def swap_e_w(d: dict[str, object]) -> dict[str, object]:
        out = dict(d)  # shallow copy

        e = out.pop("E", None)
        w = out.pop("W", None)

        if e is not None or w is not None:
            out["E"] = w
            out["W"] = e

        return out

    def plot_geo(
        self,
        plot_patches: Optional[bool] = False,
    ) -> Type[Figure]:
        from matplotlib.pyplot import subplots
        from random import random

        f, ax = subplots()
        for key, patch in self.patches.items():
            color = (random(), random(), random())
            if plot_patches:
                patch.plot_patch_geo()
            patch.plot_patch_geo(ax=ax, color=color)

        return ax.get_figure()

    def plot_patchmap(self) -> Type[Figure]:
        from matplotlib.pyplot import subplots
        f, ax = subplots()
        ax.plot([-0.5,-0.5], [-.5, self.ny+.5], 'k-')
        ax.plot([self.ixrb.max()+1.5, self.ixrb.max()+1.5], [-.5, self.ny+.5], 'k-')
        ax.plot([-.5, self.ixrb.max()+1.5], [-.5, -.5], 'k-')
        ax.plot([-.5, self.ixrb.max()+1.5], [self.ny+.5, self.ny+.5], 'k-')

        if self.xpoints > 1:
            ax.plot([self.ixlb[1]-.5, self.ixlb[1]-.5], [-.5, self.ny+.5], 'k-')
    
        for key, patch in self.patches.items():
            patch.plot_patch(ax=ax)

        return ax.get_figure()


    def _neighbors_to_list(self, entry: Any) -> list[str]:
        """
        Neighbor entry may be:
          - missing / None
          - a string: "patch_name"
          - a dict of strings: {something: "patch_name", ...}
          - a list/tuple of strings (optional support)
        Returns a list of neighbor patch names.
        """
        if entry is None:
            return []

        if isinstance(entry, str):
            return [entry]

        if isinstance(entry, dict):
            out: list[str] = []
            for v in entry.values():
                if isinstance(v, str):
                    out.append(v)
                else:
                    raise TypeError(f"Neighbor dict values must be strings, got {type(v)}")
            return out

        if isinstance(entry, (list, tuple)):
            if all(isinstance(v, str) for v in entry):
                return list(entry)
            raise TypeError("Neighbor list/tuple must contain only strings")

        raise TypeError(f"Unsupported neighbor entry type: {type(entry)}")


    def patch_bounds_from_neighbors(
        self,
        patch_starts: Mapping[str, tuple[int, int]],
        patch_info: Mapping[str, Mapping[str, Any]],
        x_pts: Sequence[int],
        y_pts: Sequence[int],
    ) -> dict[str, tuple[tuple[int, int], tuple[int, int]]]:
        """
        Compute per-patch bounds as values:
            {name: ((x0, x1), (y0, y1))}

        Assumptions:
          - patch_starts[name] = (ix_start, iy_start) uses 1-based segment indices
          - x_pts and y_pts are segmentation *points* (edges), so nx = len(x_pts)-1 segments, ny = len(y_pts)-1
          - patch_info[name]["neighbors"] is a dict with keys in {"N","S","E","W"} (case-insensitive),
            each value either a string neighbor name or a dict of strings (multiple neighbors along that edge).
          - The x-extent ends at the boundary defined by east neighbor start ix (or domain boundary if none).
          - The y-extent ends at the boundary defined by north neighbor start iy (or domain boundary if none).
        """
        nx = len(x_pts) - 1
        ny = len(y_pts) - 1
        if nx <= 0 or ny <= 0:
            raise ValueError("x_pts and y_pts must each have at least two points (>=1 segment).")

        out: dict[str, tuple[tuple[int, int], tuple[int, int]]] = {}

        for name, (ix0, iy0) in patch_starts.items():
            if not (1 <= ix0 <= nx):
                raise IndexError(f"{name}: ix_start={ix0} out of range 1..{nx}")
            if not (1 <= iy0 <= ny):
                raise IndexError(f"{name}: iy_start={iy0} out of range 1..{ny}")

            info = patch_info.get(name, {})
            neigh = info.get("neighbors", {}) or {}

            # allow lowercase keys
            def get_dir(d: str) -> Any:
                return neigh.get(d) if d in neigh else neigh.get(d.lower()) if d.lower() in neigh else neigh.get(d.upper())

            east = self._neighbors_to_list(get_dir("E"))
            north = self._neighbors_to_list(get_dir("N"))

            # Infer ix_end from east neighbors (or right boundary)
            if east:
                east_ix = []
                for nb in east:
                    if nb not in patch_starts:
                        raise KeyError(f"{name}: east neighbor {nb!r} not found in patch_starts")
                    east_ix.append(patch_starts[nb][0])
                ix1 = min(east_ix)

                if ix1 <= ix0:
                    raise ValueError(f"{name}: inferred ix_end={ix1} is not > ix_start={ix0} (east neighbors: {east})")
            else:
                ix1 = nx + 1  # right boundary (in 1-based boundary-index space)

            # Infer iy_end from north neighbors (or top boundary)
            if north:
                north_iy = []
                for nb in north:
                    if nb not in patch_starts:
                        raise KeyError(f"{name}: north neighbor {nb!r} not found in patch_starts")
                    north_iy.append(patch_starts[nb][1])
                iy1 = min(north_iy)

                if iy1 <= iy0:
                    raise ValueError(f"{name}: inferred iy_end={iy1} is not > iy_start={iy0} (north neighbors: {north})")
            else:
                iy1 = ny + 1  # top boundary

            # Convert boundary indices to coordinate values
            x0, x1 = x_pts[ix0 - 1], x_pts[ix1 - 1]
            y0, y1 = y_pts[iy0 - 1], y_pts[iy1 - 1]

            out[name] = ((x0, x1-1), (y0, y1-1))

        return out


class TopologicalPacker2D:
    """
    Assign integer (x, y) coordinates from N/S/E/W internal neighbor constraints.

    NEW behavior (requested):
      - You DO NOT need to specify wall connections in the input.
      - The class automatically assigns which patches touch each wall by detecting
        "sources" and "sinks" in the constraint graphs:
          * X graph (Gx): sources -> West wall, sinks -> East wall
          * Y graph (Gy): sources -> South wall, sinks -> North wall

    Conventions:
      - E neighbor: A -> B in Gx  (B is to the right of A)
      - W neighbor: B -> A in Gx  (B is to the left of A)
      - N neighbor: A -> B in Gy  (B is above A)
      - S neighbor: B -> A in Gy  (B is below A)

    If cycles exist, SCC condensation is used for wall attachment and ranking.
    """

    def __init__(
        self,
        patches: Mapping[Handle, Mapping[str, Any]],
        *,
        neighbor_key: str = "neighbors",
        walls: Walls = Walls(),
        boundary_aliases: Optional[Mapping[str, str]] = None,
        auto_attach_walls: bool = True,
    ) -> None:
        """
        patches:
          dict like {handle: {"neighbors": {"N": ..., "S": ..., "E": ..., "W": ...}, ...}, ...}
          Only internal patch handles are required. Wall tokens are optional.

        boundary_aliases:
          Optional mapping of any boundary tokens you might still use to canonical walls.
          (Not required for auto wall attachment, but kept for backwards compatibility.)

        auto_attach_walls:
          If True, automatically connect inferred boundary patches to walls.
        """
        self.patches = patches
        self.neighbor_key = neighbor_key
        self.walls = walls
        self.auto_attach_walls = auto_attach_walls

        default_aliases = {
            # If you still use these tokens, they can be canonicalized.
            "BOUNDARY_W": walls.W,
            "BOUNDARY_E": walls.E,
            "BOUNDARY_S": walls.S,
            "BOUNDARY_N": walls.N,
            walls.W: walls.W,
            walls.E: walls.E,
            walls.S: walls.S,
            walls.N: walls.N,
        }
        if boundary_aliases:
            default_aliases.update(boundary_aliases)
        self.boundary_aliases = default_aliases

        self.Gx: nx.DiGraph = nx.DiGraph()
        self.Gy: nx.DiGraph = nx.DiGraph()

    # ---------- Public API ----------

    def arrange(self) -> Dict[Handle, Tuple[int, int]]:
        """Return integer (x, y) coords for each patch handle."""
        self.build_graphs()
        rx = self.rank(self.Gx)
        ry = self.rank(self.Gy)
        return {h: (rx[h], ry[h]) for h in self.patches.keys()}

    def build_graphs(self) -> Tuple[nx.DiGraph, nx.DiGraph]:
        """(Re)build and store constraint graphs Gx and Gy from internal couplings."""
        self.Gx = nx.DiGraph()
        self.Gy = nx.DiGraph()

        # Add patch nodes
        for h in self.patches.keys():
            self.Gx.add_node(h)
            self.Gy.add_node(h)

        # Add internal neighbor constraints
        for a, pd in self.patches.items():
            nbd = (pd.get(self.neighbor_key) or {}) if pd else {}

            # East: a -> b in Gx
            for b in _as_list(nbd.get("E")):
                b = self._canon_if_boundary(b)
                self.Gx.add_edge(a, b)

            # West: b -> a in Gx
            for b in _as_list(nbd.get("W")):
                b = self._canon_if_boundary(b)
                self.Gx.add_edge(b, a)

            # North: a -> b in Gy
            for b in _as_list(nbd.get("N")):
                b = self._canon_if_boundary(b)
                self.Gy.add_edge(a, b)

            # South: b -> a in Gy
            for b in _as_list(nbd.get("S")):
                b = self._canon_if_boundary(b)
                self.Gy.add_edge(b, a)

        # Add wall nodes (always present)
        for w in (self.walls.W, self.walls.E, self.walls.S, self.walls.N):
            self.Gx.add_node(w)
            self.Gy.add_node(w)

        # Automatically attach walls based on constraint sources/sinks
        if self.auto_attach_walls:
            self._auto_attach_walls(self.Gx, west=self.walls.W, east=self.walls.E)
            self._auto_attach_walls(self.Gy, west=self.walls.S, east=self.walls.N)

        return self.Gx, self.Gy

    # ---------- Ranking / SCC handling ----------

    def rank(self, G: nx.DiGraph) -> Dict[Handle, int]:
        """
        If DAG: longest-path layering.
        If cyclic: condense SCCs -> DAG, layer that, then expand.
        """
        if nx.is_directed_acyclic_graph(G):
            return self._longest_path_layering(G)

        sccs = list(nx.strongly_connected_components(G))
        comp_index: Dict[Handle, int] = {}
        for i, comp in enumerate(sccs):
            for n in comp:
                comp_index[n] = i

        CG = nx.DiGraph()
        CG.add_nodes_from(range(len(sccs)))
        for u, v in G.edges():
            cu, cv = comp_index[u], comp_index[v]
            if cu != cv:
                CG.add_edge(cu, cv)

        comp_rank = self._longest_path_layering(CG)

        out: Dict[Handle, int] = {}
        for node, ci in comp_index.items():
            out[node] = comp_rank[ci]
        return out

    @staticmethod
    def _longest_path_layering(G: nx.DiGraph) -> Dict[Any, int]:
        """Layering for a DAG: rank[v] = max(rank[u] + 1 for u->v), starting at 0."""
        r: Dict[Any, int] = {n: 0 for n in G.nodes()}
        for v in nx.topological_sort(G):
            preds = list(G.predecessors(v))
            if preds:
                r[v] = max(r[p] + 1 for p in preds)
        return r

    # ---------- Auto wall attachment ----------

    def _auto_attach_walls(self, G: nx.DiGraph, *, west: str, east: str) -> None:
        """
        Attach wall nodes automatically:
          - Connect west -> sources
          - Connect sinks -> east

        For cyclic graphs, we compute SCC condensation and use its sources/sinks,
        then connect walls to all members of those SCCs.

        Notes:
          - We only consider patch nodes when determining sources/sinks.
          - Walls themselves are ignored in source/sink detection.
        """
        patch_nodes = set(self.patches.keys())

        # Nothing to do if no patches
        if not patch_nodes:
            return

        # Helper: add edges west->node and node->east for each node in iterable
        def attach_sources_sinks(sources: Iterable[Handle], sinks: Iterable[Handle]) -> None:
            for n in sources:
                if n in patch_nodes:
                    G.add_edge(west, n)
            for n in sinks:
                if n in patch_nodes:
                    G.add_edge(n, east)

        # If it's already a DAG (among all nodes), we can just compute sources/sinks directly on patches.
        # But because walls are present as nodes, filter degrees on patch-induced subgraph.
        H = G.subgraph(patch_nodes).copy()

        if nx.is_directed_acyclic_graph(H):
            sources = [n for n in H.nodes() if H.in_degree(n) == 0]
            sinks = [n for n in H.nodes() if H.out_degree(n) == 0]
            attach_sources_sinks(sources, sinks)
            return

        # Cyclic: compute SCCs on patch subgraph, condense to DAG
        sccs = list(nx.strongly_connected_components(H))
        comp_index: Dict[Handle, int] = {}
        for i, comp in enumerate(sccs):
            for n in comp:
                comp_index[n] = i

        CG = nx.DiGraph()
        CG.add_nodes_from(range(len(sccs)))
        for u, v in H.edges():
            cu, cv = comp_index[u], comp_index[v]
            if cu != cv:
                CG.add_edge(cu, cv)

        source_comps = [c for c in CG.nodes() if CG.in_degree(c) == 0]
        sink_comps = [c for c in CG.nodes() if CG.out_degree(c) == 0]

        sources: List[Handle] = []
        sinks: List[Handle] = []
        for ci in source_comps:
            sources.extend(list(sccs[ci]))
        for ci in sink_comps:
            sinks.extend(list(sccs[ci]))

        attach_sources_sinks(sources, sinks)

    # ---------- Boundary canonicalization (optional) ----------

    def _canon_if_boundary(self, h: str) -> str:
        """
        If h matches a known boundary alias (optional), map it to the canonical wall node.
        Otherwise return h unchanged.

        With auto wall attachment, you typically won't use this at all.
        """
        mapped = self.boundary_aliases.get(h, None)
        return mapped if mapped is not None else h

