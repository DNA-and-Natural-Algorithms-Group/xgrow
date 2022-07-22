from __future__ import annotations
from dataclasses import dataclass, field
from collections import UserList
import itertools
from io import StringIO
from typing import (
    Any,
    IO,
    List,
    Literal,
    Mapping,
    Optional,
    Sequence,
    TextIO,
    Dict,
    Tuple,
    TypeVar,
    Union,
    cast,
    overload,
)
from warnings import warn
from copy import deepcopy

__all__ = ["XgrowArgs", "Bond", "Tile", "Glue", "TileSet", "InitState"]

import yaml

# Classes for controlling flow/block style of various yaml elements. These are
# used in dump()
class blockseq(dict):
    pass


def blockseq_rep(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:seq", data, flow_style=False)


class blockmap(dict):
    pass


def blockmap_rep(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:map", data, flow_style=False)


class flowmap(dict):
    pass


def flowmap_rep(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:map", data, flow_style=True)


yaml.add_representer(blockseq, blockseq_rep)
yaml.add_representer(flowmap, flowmap_rep)
yaml.add_representer(blockmap, blockmap_rep)


@dataclass(init=False)
class XgrowArgs:
    """Xgrow arguments."""

    block: Optional[int] = None
    "display block size, in pixels"
    size: Optional[int] = None
    "field side length (power-of-two) [default 256]"
    rand: Optional[str] = None
    "random number seed"
    k: Optional[float] = None
    "hybridization rate constant (/sec)"
    Gmc: Optional[float] = None
    "initiation free energy  (units kT)"
    Gse: Optional[float] = None
    "interaction free energy per binding"
    Gas: Optional[float] = None
    Gam: Optional[float] = None
    Gae: Optional[float] = None
    Gfc: Optional[float] = None
    T: Optional[float] = None
    blast_rate_alpha: Optional[float] = None
    blast_rate_beta: Optional[float] = None
    blast_rate_gamma: Optional[float] = None
    seed: Optional[str] = None
    update_rate: Optional[int] = None
    tracefile: Optional[str] = None
    untiltiles: Optional[Sequence[Union[str, int]]] = None
    "Run until the named or numbered tiles are present."
    tmax: Optional[float] = None
    emax: Optional[int] = None
    smax: Optional[int] = None
    smin: Optional[int] = None
    untiltilescount: Optional[str] = None
    clean_cycles: Optional[int] = None
    error_radius: Optional[float] = None
    datafile: Optional[str] = None
    arrayfile: Optional[str] = None
    exportfile: Optional[str] = None
    importfile: Optional[str] = None
    min_strength: Optional[float] = None
    window: Optional[bool] = None
    doubletiles: Optional[Sequence[Tuple[str | int, str | int]]] = tuple()
    vdoubletiles: Optional[Sequence[Tuple[str | int, str | int]]] = tuple()
    pause: Optional[bool] = None
    wander: Optional[bool] = None
    periodic: Optional[bool] = None
    fission: Optional[Literal["off", "on", "chunk"]] = None
    font: Optional[str] = None

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> XgrowArgs:
        return cls(**d)

    def to_dict(self) -> Dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v is not None}

    def __init__(self, **kwargs: Dict[str, Any]):
        for k, v in kwargs.items():
            if k in self.__dataclass_fields__:  # type: ignore
                setattr(self, k, v)
            else:
                warn(f"Ignoring {k}={v}")

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            + ", ".join(
                f"{k}={repr(getattr(self, k))}"
                for k in self.__dataclass_fields__  # type: ignore
                if getattr(self, k) not in [None, tuple()]
            )
            + ")"
        )

    def update(self, other: Dict[str, Any] | XgrowArgs) -> None:
        if isinstance(other, dict):
            self.__dict__.update(other)
        else:
            self.__dict__.update(other.__dict__)

    def merge(self, other: Dict[str, Any] | XgrowArgs) -> XgrowArgs:
        d = deepcopy(self)
        d.update(other)
        return d

    def __bool__(self):
        return len(self.__dict__) > 0 and all(
            not bool(v) for v in self.__dict__.values()
        )


@dataclass
class Tile:
    "An xgrow tile, either single or double."
    edges: List[str | int]
    """A lists of either string bond names, or integer bond references.  0 is a null bond.
    Edges proceed clockwise, starting from the North (or Northwest for horizontal double tiles)."""
    name: Optional[str] = None
    "The tile name."
    shape: Literal["S", "H", "V"] = "S"
    "Whether the tile is a Single, Horizontal Double Tile, or Vertical Double Tile"
    stoic: Optional[float] = None
    "Stoichiometric ratio of the tile to base tile concentration (default 1.0)"
    color: Optional[str] = None
    "Color of the tile, as an X11 color name or reference"

    def _xgstring(self, tilenum: int | None = None) -> str:
        assert self.shape == "S"
        ts = "{ " + " ".join(str(x) for x in self.edges) + " }"
        if self.stoic is not None:
            ts += f"[{self.stoic}]"
        if self.color is not None:
            ts += f"({self.color})"
        if self.name is not None:
            ts += f"<{self.name}>"
        ts += "   % "
        if tilenum:
            ts += f" (tile #{tilenum})\n"
        return ts

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> Tile:
        return cls(**d)

    def to_dict(self) -> Dict[str, str | int]:
        return {k: v for k, v in self.__dict__.items() if v is not None}

    def ident(self) -> str:
        if self.name:
            return self.name
        else:
            raise ValueError


@dataclass
class Bond:
    "An xgrow bond type"
    name: Optional[str] = None
    strength: int = 1

    @property
    def xgname(self) -> str:
        raise NotImplementedError

    @classmethod
    def from_dict(cls, d: Mapping[str, str | int]) -> Bond:
        return cls(**d) # type: ignore

    def to_dict(self) -> Dict[str, str | int]:
        return {k: v for k, v in self.__dict__.items() if v is not None}


@dataclass
class Glue:
    "An xgrow glue declaration, setting the interaction strength between two bond types."
    bond1: int | str
    bond2: int | str
    strength: float

    @classmethod
    def from_dict(cls, d: Sequence[int | str | float]) -> Glue:
        return cls(*d)  # type: ignore

    def to_dict(self) -> Tuple[int | str, int | str, float]:
        return (self.bond1, self.bond2, self.strength)


def _updatebonds(
    edges: List[str | int],
    bi: int,
    bm: int,
    bondnames: Dict[str, int],
    extrabonds: List[Bond],
) -> Tuple[int, int]:
    for e in edges:
        if isinstance(e, int):
            bm = max(e, bm)
        if isinstance(e, str) and e != "0" and e not in bondnames:
            bondnames[e] = bi
            extrabonds.append(Bond(e))
            bm = max(bi, bm)
            bi += 1
    return bi, bm


class InitState(UserList):  # [Tuple[int, int, Union[str, int]]]
    """An initial state for xgrow.  This consists of a list of (x, y, tile_name_or_int) tuples."""

    def to_importfile(
        self, size: int, tilenums: Mapping[str, int], stream: Optional[IO[str]] = None
    ) -> Optional[str]:
        # fixme: check connectivity?
        if stream is None:
            stream = StringIO()

        stream.write("flake{1}={...\n")
        stream.write("[ 0 0 0 0 0 0 0 0 0 0 ],...\n")
        stream.write("[ 1 1 1 ],...\n")

        asnums = {(x, y): _get_or_int(tilenums, t) for x, y, t in self.data}

        padlen = max(len(str(t)) for t in asnums.values())

        for x in range(0, size):
            if x == 0:
                stream.write("[")
            else:
                stream.write(" ")
            for y in range(0, size):
                stream.write(" " + str(asnums.get((x, y), 0)).rjust(padlen))
            stream.write("; ...\n")
        stream.write("] };\n")

        stream.flush()

        if isinstance(stream, StringIO):
            return stream.getvalue()
        else:
            return None


@dataclass
class TileSet:
    """An xgrow tileset."""

    tiles: List[Tile]
    bonds: List[Bond] = field(default_factory=list)
    glues: List[Glue] = field(default_factory=list)
    xgrowargs: XgrowArgs = field(default_factory=XgrowArgs)
    initstate: Optional[InitState] = None

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> TileSet:
        if "xgrowargs" in d.keys():
            xga = XgrowArgs.from_dict(d["xgrowargs"])
        else:
            xga = XgrowArgs()
        return cls(
            tiles=[Tile.from_dict(x) for x in d["tiles"]],
            bonds=[Bond.from_dict(x) for x in d.get("bonds", [])],
            glues=[Glue.from_dict(x) for x in d.get("glues", [])],
            xgrowargs=xga,
            initstate=d.get("initstate", None),
        )

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {}
        for k, v in self.__dict__.items():
            if v not in [[], None]:
                if isinstance(v, InitState):
                    d[k] = [x for x in v]
                elif isinstance(v, Sequence):
                    d[k] = [x.to_dict() for x in v]  # type: ignore
                else:
                    d[k] = v.to_dict()

        return d

    @classmethod
    def from_yaml(cls, file_or_handle) -> TileSet:
        import yaml

        if isinstance(file_or_handle, str):
            try:
                file_or_handle = open(file_or_handle, "r")
            except FileNotFoundError:
                pass
        return cls.from_dict(yaml.safe_load(file_or_handle))

    def to_yaml(self, file_or_stream=None) -> None:
        import yaml

        tileset = self.to_dict()

        # Make tile and bond sections flow-style
        tileset["tiles"] = [flowmap(x) for x in tileset["tiles"]]
        tileset["bonds"] = [flowmap(x) for x in tileset["bonds"]]

        # If xgrowargs is there, make it block-style
        if "xgrowargs" in tileset.keys():
            tileset["xgrowargs"] = blockmap(tileset["xgrowargs"])
        return yaml.dump(tileset, file_or_stream)

    @overload
    def to_xgrow(
        self,
        stream: TextIO,
        extraparams: Dict[str, Any] = {},
        *,
        return_tilenums: Literal[True],
        double_tiles_at_end: bool = True,
    ) -> Tuple[None, Dict[str, int]]:
        ...

    @overload
    def to_xgrow(
        self,
        stream: None = None,
        extraparams: Dict[str, Any] = {},
        *,
        return_tilenums: Literal[True],
        double_tiles_at_end: bool = True,
    ) -> Tuple[str, Dict[str, int]]:
        ...

    @overload
    def to_xgrow(
        self,
        stream: None = None,
        extraparams: Dict[str, Any] = {},
        return_tilenums: Literal[False] = False,
        double_tiles_at_end: bool = True,
    ) -> str:
        ...

    @overload
    def to_xgrow(
        self,
        stream: TextIO,
        extraparams: Dict[str, Any] = {},
        return_tilenums: Literal[False] = False,
        double_tiles_at_end: bool = True,
    ) -> None:
        ...

    def to_xgrow(
        self,
        stream: Optional[TextIO] = None,
        extraparams: Dict[str, Any] = {},
        return_tilenums: bool = False,
        double_tiles_at_end: bool = True,
    ) -> None | str | Tuple[str, Dict[str, int]] | Tuple[None, Dict[str, int]]:
        if stream is None:
            stream = StringIO()
            writing = False
        else:
            writing = True

        tstrings: list[str] = []
        ti = 1
        tile_to_i: Dict[str, int] = dict()
        bondnames: Dict[str, int] = dict()
        extrabonds: list[Bond] = []
        bi = 1
        bm = 0
        hdoubles: list[Tuple[int, int]] = []
        vdoubles: list[Tuple[int, int]] = []
        doubletilestrings: list[str] = []
        double_i = len(self.tiles) + 1

        for b in self.bonds:
            if b.name is None:
                if len(bondnames):
                    raise ValueError
                else:
                    continue
            bondnames[b.name] = bi
            bi += 1
            bm += 1

        for tile in self.tiles:
            e = tile.edges
            if tile.shape == "S":
                tstrings.append(tile._xgstring(ti))
                if tile.name:
                    tile_to_i[tile.name] = ti
                bi, bm = _updatebonds(e, bi, bm, bondnames, extrabonds)
                ti += 1
            else:
                if tile.shape == "H":
                    t1 = Tile(
                        [e[0], f"{tile.name}_db", e[4], e[5]],
                        tile.name,
                        "S",
                        tile.stoic,
                        tile.color,
                    )
                    t2 = Tile(
                        [e[1], e[2], e[3], f"{tile.name}_db"],
                        _def(tile.name) + "_right",
                        "S",
                        tile.stoic,
                        tile.color,
                    )
                    if double_tiles_at_end:
                        hdoubles.append((ti, double_i))

                    else:
                        hdoubles.append((ti, ti + 1))
                elif tile.shape == "V":
                    t1 = Tile(
                        [e[0], e[1], f"{tile.name}_db", e[5]],
                        tile.name,
                        "S",
                        tile.stoic,
                        tile.color,
                    )
                    t2 = Tile(
                        [f"{tile.name}_db", e[2], e[3], e[4]],
                        _def(tile.name) + "_bottom",
                        "S",
                        tile.stoic,
                        tile.color,
                    )
                    if double_tiles_at_end:
                        vdoubles.append((ti, double_i))
                    else:
                        vdoubles.append((ti, ti + 1))
                else:
                    raise ValueError
                tstrings.append(t1._xgstring(ti))
                if t1.name:
                    tile_to_i[t1.name] = ti
                bi, bm = _updatebonds(t1.edges, bi, bm, bondnames, extrabonds)
                ti += 1
                if double_tiles_at_end:
                    doubletilestrings.append(t2._xgstring(double_i))
                    if t2.name:
                        tile_to_i[t2.name] = double_i
                    bi, bm = _updatebonds(t2.edges, bi, bm, bondnames, extrabonds)
                    double_i += 1
                else:
                    tstrings.append(t2._xgstring(ti))
                    if t2.name:
                        tile_to_i[t2.name] = ti
                    bi, bm = _updatebonds(t2.edges, bi, bm, bondnames, extrabonds)
                    ti += 1

        stream.write(f"num tile types={double_i-1}\n" f"num binding types={bm}\n")

        if len(bondnames):
            stream.write(
                "binding type names={ " + " ".join(b for b in bondnames) + " }\n"
            )

        stream.write("tile edges={\n")

        for t in itertools.chain(tstrings, doubletilestrings):
            stream.write(t)

        stream.write("}\n")

        stream.write("binding strengths={ ")

        stream.write(" ".join(str(b.strength) for b in self.bonds + extrabonds))

        if (nb := bm - len(self.bonds) - len(extrabonds)) > 0:
            stream.write(" " + " ".join("1" for _ in range(0, nb)))

        stream.write(" }\n")

        for glue in self.glues:
            n1 = _get_or_int(bondnames, glue.bond1)
            n2 = _get_or_int(bondnames, glue.bond2)
            stream.write(f"g({n1},{n2})={glue.strength:g}\n")

        xgrowargs = self.xgrowargs.merge(extraparams)
        for f in xgrowargs.__dataclass_fields__:  # type: ignore
            v = getattr(xgrowargs, f)
            if v is None:
                continue
            if f == "doubletiles":
                for i1, i2 in v:
                    i1 = _get_or_int(tile_to_i, i1)
                    i2 = _get_or_int(tile_to_i, i2)
                    stream.write(f"doubletile={i1},{i2}\n")
            elif f == "vdoubletiles":
                for i1, i2 in v:
                    i1 = _get_or_int(tile_to_i, i1)
                    i2 = _get_or_int(tile_to_i, i2)
                    stream.write(f"vdoubletile={i1},{i2}\n")
            elif isinstance(v, (list, tuple)):
                v = ",".join(str(x) for x in v)
                stream.write(f"{f}={v}\n")
            else:
                stream.write(f"{f}={v}\n")

        for i1, i2 in hdoubles:
            stream.write(f"doubletile={i1},{i2}\n")
        for i1, i2 in vdoubles:
            stream.write(f"vdoubletile={i1},{i2}\n")

        if writing:
            if return_tilenums:
                return None, tile_to_i
            else:
                return None
        else:
            if return_tilenums:
                return cast(StringIO, stream).getvalue(), tile_to_i
            else:
                return cast(StringIO, stream).getvalue()


def _get_or_int(d: Mapping[str, int], v: int | str):
    if isinstance(v, int):
        return v
    else:
        return d[v]


T = TypeVar("T")


def _def(k: T | None) -> T | str:
    if k is None:
        return ""
    else:
        return k
