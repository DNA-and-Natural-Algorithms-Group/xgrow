# Notes on the STXG input format

STXG files are in [YAML](yaml.org) format. Their top level should be a mapping containing at least a `tiles` key, which will contain a sequence of tiles.  It can also contain `bonds`, `glues`, and `xgrowargs`. (The difference between `bonds` and `glues` is historical, sorry.  For the purposes of Xgrow, "bond" refers to a type of glue, while "glue" refers to a specific binding strength between two potentially non-identical bond types.)

## Tiles

The `tiles` key is *required*, as it contains the tiles for the tileset.  It should be a sequence of mappings, containing the following keys, of which only `edges` is required:

- (*required*) `edges`: a sequence of 4 bond types, referenced by either name or number (names thus can't start with a digit). The number `0` is special, and refers to a null bond. Other than `0`, I recommend against using numbered glues. Using a name that is not included in the `bonds` key will automatically create a new strength-1 bond.
- `stoic`: the ratio of the tile's concentration to the standard tile concentration.  This can be zero if, eg, you want to use a seed tile that will never be able to attach to a flake.
- `color`: the color of the tile within Xgrow. This can be a X11 color name, or a "#RRGGBB"-hex-format color.
- `name`: a tile name. This has no effect within Xgrow, but is included in the comments of the (old-format) input file that is generated. It is also used for double tiles.

## Bonds

The `bonds` key is *optional*: if all bonds are strength-1 and have no special properties, it is sufficient to let Xgrow create them automatically.  If you want to set options for them, however, `bonds` contains a sequence of mappings with the following keys:

- (*required*) `name`: the name of the bond.  This is required, so that the bonds in `tiles` can refer to it. It *is* used by Xgrow, and thus is somewhat restricted in what it can contain, but letters, numbers and underscores should be fine, as long as the first character is not a number.
- `strength`: the strength of the bond. You may want to set this to `0` if you want only non-diagonal bonds.

## Glues

`glues` (*optional*) is used to set non-diagonal bond interactions.  It should be a sequence of sequences of format `[bond_name_1, bond_name_2, strength]`.

## Xgrowargs

`xgrowargs` (*required* at the moment, but shouldn't be!) is a mapping that includes arguments to xgrow, as key-value pairs.  Arguments that would normally have no value use `true` and `false` to set and unset them.

One difference between this and Xgrow itself is with the `doubletiles` and `vdoubletiles` arguments.  These should be sequences of sequences of the format `[tile_name_1, tile_name_2]`, rather than needing to use the tile number.

## An example

This example is a modification of sierpinski.stxg in the examples folder in order to include more options:

```yaml
bonds:
- { name: B, strength: 2 }

tiles:
- { name: corner, edges: [ B, 0, 0, B ], color: red, stoic: 0  }
- { name: topboundary, edges: [B, 0, B, v1], color: magenta, stoic: 0.5 }
- { name: leftboundary, edges: [v1, B, 0, B], color: purple, stoic: 0.5 }
- { edges: [ v0, v0, v0, v0 ], color: blue3 }
- { edges: [ v0, v1, v1, v0 ], color: green }
- { edges: [ v1, v0, v1, v1 ], color: yellow }
- { edges: [ v1, v1, v0, v1 ], color: tan }

xgrowargs:
    block: 8
    size: 64
    Gse: 8.03
    Gmc: 16
    update_rate: 5000
    pause: true
```
