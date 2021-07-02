from __future__ import annotations
from dataclasses import dataclass, field
from io import StringIO
from typing import Any, List, Literal, Optional, Sequence, TextIO, Tuple, cast, overload
from warnings import warn

__all__ = ["XgrowArgs", "Bond", "Tile", "Glue", "TileSet"]


@dataclass(init=False)
class XgrowArgs:
    block: Optional[int] = None
    size: Optional[int] = None
    rand: Optional[str] = None
    k: Optional[float] = None
    Gmc: Optional[float] = None
    Gse: Optional[float] = None
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
    untiltiles: Optional[str] = None
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
    window: Optional[bool] = None

    @classmethod
    def from_dict(cls, d: dict) -> XgrowArgs:
        return cls(**d)

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v}

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if k in self.__dataclass_fields__:
                setattr(self, k, v)
            else:
                warn(f"Ignoring {k}={v}")

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            + ", ".join(
                f"{k}={repr(getattr(self, k))}"
                for k in self.__dataclass_fields__
                if getattr(self, k) not in [None, tuple()]
            )
            + ")"
        )

    def __bool__(self):
        return len(self.__dict__) > 0 and all(
            not bool(v) for v in self.__dict__.values()
        )


@dataclass
class Tile:
    edges: List[str | int]
    name: Optional[str] = None
    shape: Literal["S", "H", "V"] = "S"
    stoic: Optional[float] = None
    color: Optional[str] = None

    def _xgstring(self, tilenum=None) -> str:
        assert self.shape == "S"
        ts = "{ " + " ".join(str(x) for x in self.edges) + " }"
        if self.stoic is not None:
            ts += f"[{self.stoic}]"
        if self.color is not None:
            ts += f"({self.color})"
        ts += "   % "
        if self.name is not None:
            ts += self.name
        if tilenum:
            ts += f" (tile #{tilenum})\n"
        return ts

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v}

    @classmethod
    def from_dict(cls, d: dict) -> Tile:
        return cls(**d)


@dataclass
class Bond:
    name: Optional[str] = None
    strength: int = 1

    @property
    def xgname(self) -> str:
        raise NotImplementedError

    @classmethod
    def from_dict(cls, d: dict) -> Bond:
        return cls(**d)

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v}


@dataclass
class Glue:
    bond1: int | str
    bond2: int | str
    strength: float

    @classmethod
    def from_dict(cls, d: list) -> Glue:
        return cls(*d)

    def to_dict(self) -> tuple[int | str, int | str, float]:
        return (self.bond1, self.bond2, self.strength)


def _updatebonds(
    edges: List[str | int],
    bi: int,
    bm: int,
    bondnames: dict[str, int],
    extrabonds: List[Bond],
) -> tuple[int, int]:
    for e in edges:
        if isinstance(e, int):
            bm = max(e, bm)
        if isinstance(e, str) and e != "0" and e not in bondnames:
            bondnames[e] = bi
            extrabonds.append(Bond(e))
            bm = max(bi, bm)
            bi += 1
    return bi, bm


@dataclass
class TileSet:
    tiles: List[Tile]
    bonds: List[Bond] = field(default_factory=list)
    glues: List[Glue] = field(default_factory=list)
    xgrowargs: XgrowArgs = XgrowArgs()

    @classmethod
    def from_dict(cls, d: dict) -> TileSet:
        if "xgrowargs" in d.keys():
            xga = XgrowArgs.from_dict(d["xgrowargs"])
        else:
            xga = XgrowArgs()
        return cls(
            tiles=[Tile.from_dict(x) for x in d["tiles"]],
            bonds=[Bond.from_dict(x) for x in d.get("bonds", [])],
            glues=[Glue.from_dict(x) for x in d.get("glues", [])],
            xgrowargs=xga,
        )

    def to_dict(self) -> dict[str, Any]:
        d = {}
        for k, v in self.__dict__.items():
            if v:
                if isinstance(v, Sequence):
                    d[k] = [x.to_dict() for x in v]
                else:
                    d[k] = v.to_dict()
        return d

    @overload
    def to_xgrow(self, stream: None = None) -> str:
        ...

    @overload
    def to_xgrow(self, stream: TextIO) -> None:
        ...

    def to_xgrow(self, stream: Optional[TextIO] = None):
        if stream is None:
            stream = StringIO()
            writing = False
        else:
            writing = True

        tstrings = []
        ti = 1
        tile_to_i = dict()
        bondnames = dict()
        extrabonds = []
        bi = 1
        bm = 1
        hdoubles = []
        vdoubles = []

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
                    vdoubles.append((ti, ti + 1))
                else:
                    raise ValueError
                tstrings.append(t1._xgstring(ti))
                if t1.name:
                    tile_to_i[t1.name] = ti
                bi, bm = _updatebonds(t1.edges, bi, bm, bondnames, extrabonds)
                ti += 1
                tstrings.append(t2._xgstring(ti))
                if t2.name:
                    tile_to_i[t2.name] = ti
                bi, bm = _updatebonds(t2.edges, bi, bm, bondnames, extrabonds)
                ti += 1

        stream.write(f"num tile types={ti-1}\n" f"num binding types={bm}\n")  # fixme

        if len(bondnames):
            stream.write(
                "binding type names={ " + " ".join(b for b in bondnames) + " }\n"
            )

        stream.write("tile edges={\n")

        for tile in tstrings:
            stream.write(tile)

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

        for f in self.xgrowargs.__dataclass_fields__:
            v = getattr(self.xgrowargs, f)
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
            else:
                stream.write(f"{f}={v}\n")

        for i1, i2 in hdoubles:
            stream.write(f"doubletile={i1},{i2}\n")
        for i1, i2 in vdoubles:
            stream.write(f"vdoubletile={i1},{i2}\n")

        if writing:
            return None
        else:
            return cast(StringIO, stream).getvalue()


def _get_or_int(d: dict[str, int], v: int | str):
    if isinstance(v, int):
        return v
    else:
        return d[v]


def _def(k):
    if k is None:
        return ""
    else:
        return k
