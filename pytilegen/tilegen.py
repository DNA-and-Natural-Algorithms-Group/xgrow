#!/bin/env python

class Tileset:
	"""A generator for xgrow tilesets. This takes various forms of inputs, and
	outputs (with self.generate_tileset) an xgrow-compatible tileset. One
	needs, essentially, to add tiles and glues, but tiles can be named, edges
	can be named and can have arbitrary interactions, and so on. Note that ' is
	considered a special character in bond names, so that the complementarize
	option can, when False, use bonds without complements using the same data
	as otherwise. I suggest using the csv module for reading data from
	spreadsheets; it seems to be a very simple way of working with the data.
	""" 
	
	def __init__(s, complementarize=False, options=[]):
		# Various options
		s.options = options

        s.tiles = []
        s.dtiles = []
        s.ends = []
        s.glues = {}
        s.compl = complementarize

    def add_tiles_from_list( s, tilelist ):
        for t in tilelist:
            if not s.compl:
                t[0]=t[0].rstrip("'")
                t[1]=t[1].rstrip("'")
                t[2]=t[2].rstrip("'")
                t[3]=t[3].rstrip("'")
            
            # Add ends if they aren't already added. We'll sort them later.
            for e in t[0:4]: 
                if e not in s.ends and e != '0': s.ends += [e]
            
            opts = {}
            if len(t)>4 and t[4]: opts['conc'] = t[4]
            if len(t)>5 and t[5]: opts['color'] = t[5]
            if len(t)>6 and t[6]: opts['name'] = t[6]

            s.tiles += [( (t[0],t[1],t[2],t[3]), opts )]
    
    def add_double_tiles_from_list( s,tilelist ):
        "Format of list is [tl tr r br bl l color concentration]"
        for t in tilelist:
            if not s.compl:
                t[0]=t[0].rstrip("'")
                t[1]=t[1].rstrip("'")
                t[2]=t[2].rstrip("'")
                t[3]=t[3].rstrip("'")
                t[4]=t[4].rstrip("'")
                t[5]=t[5].rstrip("'")
            
            for e in t[0:6]: 
                if e not in s.ends and e != '0': s.ends += [e]

            opts = {}
            if len(t)>6 and t[6]: opts['conc'] = t[6]
            if len(t)>7 and t[7]: opts['color'] = t[7]
            if len(t)>8 and t[8]: opts['name'] = t[8]

            s.dtiles += [( (t[0],t[1],t[2],t[3],t[4],t[5]), opts )]

    def glues_from_arrays( s,ends,array ):
        error("not done!")

    def glues_from_list( s, gluelist ):
        for g in gluelist:
            if not s.compl:
                g[0]=g[0].rstrip("'")
                g[1]=g[1].rstrip("'")
            s.glues[(g[0],g[1])] = g[2]

    
    def glues_from_seqlist( s, seqlist ):
        error("not done!")

    def generate_tileset(s):
        ts = ""

        # We should sort ends here
        def endcmp(a,b):
            p = (a[-1]=="'")-(b[-1]=="'")
            if p: return p
            a=a.rstrip("'"); b=b.rstrip("'")
            p = b.isdigit() - a.isdigit()
            if p: return p
            if a.isdigit() and b.isdigit(): return cmp(int(a),int(b))
            return cmp(a,b)
        
        s.ends.sort(cmp=endcmp)

        out_ends = s.ends
        out_tiles = s.tiles

        # We need to deal with double tiles first.
        dtbn = 1
        dttn = len(s.tiles)
        dpairs = []
        for dt in s.dtiles:
            bn = "dtb"+str(dtbn)
            
            dpairs += [(dttn,dttn+1)]
            out_ends += [bn]
            out_tiles += [( (dt[0][0], bn, dt[0][4], dt[0][5]), dt[1] )]
            out_tiles += [( (dt[0][1], dt[0][2], dt[0][3], bn), dt[1] )]
            dttn += 2
        
        ts += "num tile types=%d\n" % len(out_tiles)
        ts += "num binding types=%d\n" % len(out_ends)
        ts += "tile edges={\n"
        out_ends = ['0'] + out_ends
        for t in out_tiles:
            tline = "{%d %d %d %d}" % tuple([ out_ends.index(x) for x in t[0] ])
            if 'conc' in t[1]: tline+="[%s]" % t[1]['conc']
            if 'color' in t[1]: tline+="(%s)" % t[1]['color']
            if 'name' in t[1]: tline+="\t\t\t%% %s" % t[1]['name']
            ts += tline + "\n"
        
        ts += "}\n"
        # We set all strengths with glues.
        ts += "binding strengths={" + ( "0 "*(len(out_ends)-1) ).rstrip() + "}\n"
        
        for key in s.glues:
            ts += "g(%d,%d)=%s\n" % ( out_ends.index(key[0]), out_ends.index(key[1]), s.glues[key] )
        
        ts += "\n"
        
        # Add the double tiles
        for pair in dpairs:
            ts += "doubletile=%d,%d\n" % tuple(x+1 for x in pair)
        
        # Add the options
        for option in s.options:
            ts += option+"\n"
        
        return ts
